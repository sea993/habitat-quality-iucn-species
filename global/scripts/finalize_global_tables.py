"""Verify completed outputs and summarize observed results, without editing manuscripts."""
import json
import sqlite3
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from global_common import ROOT, OUT, CFG
TABLE=OUT/'tables'

def main():
    if not (OUT/'progress/statistics_complete.json').exists():
        raise RuntimeError('Global statistics have not finished')
    cfg=CFG
    expected=len(pd.read_csv(ROOT/'tables/global_species_screened.csv'))
    db=sqlite3.connect(OUT/'species/metrics.sqlite')
    counts=dict(db.execute('SELECT status, COUNT(*) FROM completed GROUP BY status'))
    if sum(counts.values())!=expected:
        raise ValueError('Incomplete species processing')
    invalid=db.execute('''SELECT COUNT(*) FROM metrics WHERE mean_hq NOT BETWEEN 0 AND 1
        OR q05_hq NOT BETWEEN 0 AND 1 OR low_quality_pct NOT BETWEEN 0 AND 100.0000001
        OR low_quality_pct_0239 NOT BETWEEN 0 AND 100.0000001
        OR valid_area_km2<=0 OR intersecting_cells<=0''').fetchone()[0]
    if invalid:
        raise ValueError(f'{invalid} invalid metric rows')
    incomplete=db.execute('''SELECT COUNT(*) FROM
        (SELECT species_id,scope FROM metrics GROUP BY species_id,scope
        HAVING COUNT(DISTINCT year)!=21 OR MIN(year)!=2000 OR MAX(year)!=2020
        OR MAX(valid_area_km2)-MIN(valid_area_km2)>0.000001
        OR MAX(intersecting_cells)!=MIN(intersecting_cells))''').fetchone()[0]
    if incomplete:
        raise ValueError(f'{incomplete} incomplete or spatially inconsistent trajectories')
    analyzed=db.execute("SELECT COUNT(DISTINCT species_id) FROM metrics WHERE scope='Global'").fetchone()[0]
    if analyzed!=counts.get('analyzed',0):
        raise ValueError('Global species count differs from extraction status')
    sample=pd.read_sql_query('''SELECT taxon,category,COUNT(DISTINCT species_id) AS n_species
        FROM metrics WHERE scope='Global' GROUP BY taxon,category ORDER BY taxon,category''',db)
    sample.to_csv(TABLE/'final_sample_composition.csv',index=False,encoding='utf-8-sig')
    hfp_unequal=db.execute('''SELECT COUNT(*) FROM hfp a JOIN hfp b ON a.species_id=b.species_id
        WHERE a.year=2000 AND b.year=2020 AND ABS(a.valid_area_km2-b.valid_area_km2)>0.000001''').fetchone()[0]
    if hfp_unequal or not (OUT/'progress/hfp_common_support_complete.json').exists():
        raise ValueError('HFP endpoints do not yet use verified common support')
    db.close()
    effects=pd.read_csv(TABLE/'group_effects_ci.csv')
    annual=pd.read_csv(TABLE/'group_annual_trajectories.csv')
    if not np.isfinite(effects.estimate).all():
        raise ValueError('Nonfinite effect estimate')
    relevant=effects[effects.n_species>1]
    if not np.isfinite(relevant[['ci95_low','ci95_high']]).all().all():
        raise ValueError('Missing confidence interval for a multi-species group')
    if (relevant.ci95_low>relevant.ci95_high).any():
        raise ValueError('Reversed confidence interval')
    if effects.duplicated(['scope','group','effect']).any():
        raise ValueError('Duplicate group effect')
    if annual.duplicated(['scope','group','year','metric']).any():
        raise ValueError('Duplicate group-year metric')
    background=pd.read_csv(TABLE/'background_annual.csv')
    glob=background[background.scope=='Global'].sort_values('year')
    if list(glob.year)!=cfg['years']:
        raise ValueError('Missing global background year')
    threshold=json.loads((TABLE/'global_low_quality_threshold.json').read_text())
    bridge=pd.read_csv(TABLE/'china_matched_species_comparison.csv')
    comparisons=[]
    for metric,old,new,unit in [
        ('Mean HQ','Q_mean_OLS_slope','mean_hq_slope','HQ/year'),
        ('5th-percentile HQ','Q_5pct_OLS_slope','q05_hq_slope','HQ/year'),
        ('Low-quality area, same 0.239 threshold','R_low_OLS_slope','low_quality_pct_0239_slope','percentage points/year')]:
        values=bridge[[old,new]].dropna()
        x,y=values[old].to_numpy(),values[new].to_numpy()
        rng=np.random.default_rng(cfg['random_seed'])
        delta=y-x
        boot=np.array([delta[rng.integers(0,len(delta),len(delta))].mean() for _ in range(cfg['bootstrap_replicates'])]) if len(delta) else np.array([])
        ci=np.quantile(boot,[.025,.975]) if len(boot) else [np.nan,np.nan]
        comparisons.append({'metric':metric,'n_matched_species':len(values),'original_mean_slope':x.mean() if len(x) else np.nan,
            'global_flow_mean_slope':y.mean() if len(y) else np.nan,'mean_paired_slope_difference':delta.mean() if len(delta) else np.nan,
            'ci95_low':ci[0],'ci95_high':ci[1],'same_slope_direction_fraction':np.mean(np.sign(x)==np.sign(y)) if len(x) else np.nan,
            'spearman_rho':spearmanr(x,y).statistic if len(x)>2 else np.nan,'units':unit})
    pd.DataFrame(comparisons).to_csv(TABLE/'china_three_metric_bridge.csv',index=False,encoding='utf-8-sig')
