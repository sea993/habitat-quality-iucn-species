from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats
import scipy

import argparse
ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Recompute pairwise bootstrap intervals and annual trend tests.")
parser.add_argument('--output-dir', type=Path, default=ROOT.parents[2] / 'review-output/china')
out = parser.parse_args().output_dir
out.mkdir(parents=True, exist_ok=True)
d=pd.read_csv(ROOT/'species.csv').sort_values('id_no')
ref=pd.read_csv(ROOT/'MannWhitney_Tests.csv')
metrics=dict(zip(ref.metric.unique(),['Q_mean_period_mean','Q_5pct_period_mean','R_low_period_mean','ANCI_2000','ANCI_2020','ANCI_change']))
assert len(d)==1820 and d.id_no.is_unique and len(ref)==96
B=10000
for i,r in ref.iterrows():
    x=d.loc[d[r.grouping]==r.group1,metrics[r.metric]].to_numpy()
    y=d.loc[d[r.grouping]==r.group2,metrics[r.metric]].to_numpy()
    assert len(x)==r.n1 and len(y)==r.n2
    u=stats.mannwhitneyu(x,y,method='asymptotic',alternative='two-sided',use_continuity=True)
    assert np.isclose(u.statistic,r.statistic,rtol=0,atol=1e-8)
    assert np.isclose(u.pvalue,r.p_raw,rtol=1e-8,atol=1e-15)
    assert np.isclose(2*u.statistic/(len(x)*len(y))-1,r.effect,rtol=0,atol=1e-12)
    seed=20261110+i
    rng=np.random.default_rng(seed); boot=np.empty(B)
    for start in range(0,B,200):
        n=min(200,B-start)
        xx=x[rng.integers(0,len(x),(n,len(x)))]; yy=y[rng.integers(0,len(y),(n,len(y)))]
        ranks=stats.rankdata(np.concatenate([xx,yy],axis=1),axis=1,method='average')
        uu=ranks[:,:len(x)].sum(axis=1)-len(x)*(len(x)+1)/2
        boot[start:start+n]=2*uu/(len(x)*len(y))-1
    lo,hi=np.quantile(boot,[.025,.975],method='linear')
    assert -1<=lo<=hi<=1
    ref.loc[i,['ci_low','ci_high']]=lo,hi
    ref.loc[i,'bootstrap_seed']=seed
    if (i+1)%16==0: print(f'S10B: {i+1}/96 completed',flush=True)
for _,g in ref.groupby(['metric','grouping']):
    assert np.allclose(stats.false_discovery_control(g.p_raw.to_numpy()),g.p_BH,rtol=1e-7,atol=1e-15)
ref.to_csv(out/'MannWhitney_recomputed.csv',index=False)

annual=pd.read_csv(ROOT/'annual_groups.csv'); national=pd.read_csv(ROOT/'annual_national.csv')
tr=pd.read_csv(ROOT/'annual_trend_source.csv'); series=[]
for i,r in tr.iterrows():
    if r.group=='China':
        z=national.sort_values('year'); y=z['mean'].to_numpy()
    else:
        typ={'All IUCN Red List species':'All','IUCN Red List category':'Category','Taxonomic group':'Taxon'}[r.grouping]
        group='All' if r.group=='All species' else r.group
        z=annual[(annual.group_type==typ)&(annual.group==group)].sort_values('year')
        y=z[r.metric_code+'_mean'].to_numpy()
    x=z.year.to_numpy()
    assert np.array_equal(x,np.arange(2000,2021)) and len(y)==21
    fit=stats.linregress(x,y)
    s=sum(np.sign(y[j:]-y[j-1]).sum() for j in range(1,len(y)))
    _,ties=np.unique(y,return_counts=True)
    var=(21*20*47-np.sum(ties*(ties-1)*(2*ties+5)))/18
    p=2*stats.norm.sf(abs((s-np.sign(s))/np.sqrt(var)))
    assert np.isclose(fit.slope*10,r.slope_per_decade,rtol=1e-9,atol=1e-12)
    assert np.isclose(s/210,r.tau,rtol=0,atol=1e-12) and np.isclose(p,r.mk_p,rtol=1e-8,atol=1e-15)
    margin=stats.t.ppf(.975,19)*fit.stderr*10
    sen=stats.theilslopes(y,x,alpha=.95)
    tr.loc[i,['ols_ci_low','ols_ci_high','sen_slope_per_decade','sen_ci_low','sen_ci_high']]=[
        fit.slope*10-margin,fit.slope*10+margin,sen.slope*10,sen.low_slope*10,sen.high_slope*10]
    series.extend(dict(row=int(i),year=int(xx),value=float(yy)) for xx,yy in zip(x,y))
tr.to_csv(out/'Annual_Trend_recomputed.csv',index=False)
pd.DataFrame(series).to_csv(out/'annual_trend_series.csv',index=False)
qa={'original_U_effect_p_BH_reproduced':True,'original_slope_tau_MK_p_reproduced':True,
    'S10B_rows':len(ref),'S13_rows':len(tr),'bootstrap_replicates':B,'seed_start':20261110,'seed_end':20261110+95,
    'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,
    'S10B_method':'Group-stratified species bootstrap, percentile CI, independent resampling of two groups, unadjusted pointwise intervals',
    'S13_method':'OLS t interval df=19; Sen slope interval from scipy.stats.theilslopes; no serial-correlation correction'}
(out/'statistics_integrity.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa,indent=2),flush=True)
