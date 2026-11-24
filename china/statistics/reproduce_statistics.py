"""Recompute China rank tests, pairwise comparisons, scale diagnostics and annual trends."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
FIELDS = {
    'Mean habitat quality, 2000-2020': 'Q_mean_period_mean',
    '5th-percentile habitat quality, 2000-2020': 'Q_5pct_period_mean',
    'Proportion of low-quality habitat, 2000-2020': 'R_low_period_mean',
    'ANCI in 2000': 'ANCI_2000', 'ANCI in 2020': 'ANCI_2020',
    'ANCI change, 2020 minus 2000': 'ANCI_change',
}


def report(out, name, result, keys, fields):
    reference = pd.read_csv(ROOT / (name+'.csv')).set_index(keys, verify_integrity=True)
    table = pd.DataFrame(result).set_index(keys, verify_integrity=True).sort_index()
    assert table.index.equals(reference.sort_index().index)
    np.testing.assert_allclose(table[fields], reference.loc[table.index, fields], rtol=2e-8, atol=1e-12)
    table.to_csv(out / (name+'.csv'))
    print(f'{name}: {len(table)} rows verified')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT.parents[2] / 'review-output/china')
    out = parser.parse_args().output_dir
    out.mkdir(parents=True, exist_ok=True)
    species = pd.read_csv(ROOT / 'species.csv').sort_values('id_no')
    roster = pd.read_csv(ROOT.parent / 'parameters/species_screening.csv')
    assert len(roster) == 1830 and len(species) == 1820 and species.id_no.is_unique
    assert set(species.id_no) == set(roster.loc[roster.included, 'id_no'])
    rows = []
    for label, field in FIELDS.items():
        for group in ['category', 'class']:
            subsets = [frame[field] for _, frame in species.groupby(group)]
            h, p = stats.kruskal(*subsets)
            n, k = len(species), len(subsets)
            rows.append(dict(metric=label, grouping=group, n=n,H=h,df=k-1,p_raw=p,
                             rank_eta_squared=max(0,(h-k+1)/(n-k))))
    report(out, 'KruskalWallis_Tests', rows, ['metric','grouping'], ['n','H','df','p_raw','rank_eta_squared'])
    for name in ['Welch_Tests','MannWhitney_Tests']:
        keys = ['metric','grouping','group1','group2']
        rows = []
        for r in pd.read_csv(ROOT / (name+'.csv'))[keys].itertuples(index=False):
            x = species.loc[species[r.grouping] == r.group1, FIELDS[r.metric]]
            y = species.loc[species[r.grouping] == r.group2, FIELDS[r.metric]]
            result = dict(r._asdict(),n1=len(x),n2=len(y))
            if name == 'Welch_Tests':
                v, w = x.var()/len(x), y.var()/len(y)
                df = (v+w)**2/(v*v/(len(x)-1)+w*w/(len(y)-1))
                effect = x.mean()-y.mean()
                t = effect / np.sqrt(v+w)
                radius = stats.t.ppf(.975,df)*np.sqrt(v+w)
                result.update(statistic=t,df=df,p_raw=2*stats.t.sf(abs(t),df),effect=effect,
                              ci_low=effect-radius,ci_high=effect+radius)
            else:
                u,p = stats.mannwhitneyu(x,y,method='asymptotic',use_continuity=True)
                result.update(statistic=u,p_raw=p,effect=2*u/(len(x)*len(y))-1)
            rows.append(result)
        data = pd.DataFrame(rows)
        if name == 'MannWhitney_Tests':
            data['p_BH'] = data.groupby(['metric','grouping']).p_raw.transform(stats.false_discovery_control)
        report(out,name,data,keys,[f for f in data if f not in keys])
    rows = []
    annual = pd.read_csv(ROOT/'annual_groups.csv')
    national = pd.read_csv(ROOT/'annual_national.csv')
    for r in pd.read_csv(ROOT/'annual_trend_source.csv').itertuples():
        if r.group=='China':
            series=national.set_index('year')['mean']
        else:
            code={'All IUCN Red List species':'All','IUCN Red List category':'Category','Taxonomic group':'Taxon'}[r.grouping]
            group='All' if r.group=='All species' else r.group
            series=annual[(annual.group_type==code)&(annual.group==group)].set_index('year')[r.metric_code+'_mean']
        series=series.sort_index()
        x=series.index.to_numpy()
        y=series.to_numpy()
        assert np.array_equal(x,np.arange(2000,2021))
        ols=stats.linregress(x,y)
        sen=stats.theilslopes(y,x)
        sign=sum(np.sign(y[j:]-y[j-1]).sum() for j in range(1,len(y)))
        counts=np.unique(y,return_counts=True)[1]
        variance=(21*20*47-np.sum(counts*(counts-1)*(2*counts+5)))/18
        margin=stats.t.ppf(.975,19)*ols.stderr*10
        rows.append(dict(metric=r.metric,grouping=r.grouping,group=r.group,tau=sign/210,
                         mk_p=2*stats.norm.sf(abs((sign-np.sign(sign))/np.sqrt(variance))),
                         slope_per_decade=10*ols.slope,ols_ci_low=10*ols.slope-margin,
                         ols_ci_high=10*ols.slope+margin,sen_slope_per_decade=10*sen.slope,
                         sen_ci_low=10*sen.low_slope,sen_ci_high=10*sen.high_slope))
    report(out,'Annual_Trend_Tests',rows,['metric','grouping','group'],list(rows[0])[3:])
    scale=pd.read_csv(ROOT/'multiscale_model_metrics.csv')
    scale['split']=scale.iloc[:,0].ffill().map({'训练集':'train','测试集':'test'})
    scale=scale.rename(columns={'Scale':'scale','Morans_I':'moran_i','P_value':'p_value'})
    report(out,'Multiscale_Evaluation',scale,['scale','split'],['moran_i','p_value'])
    stored=pd.read_csv(ROOT.parent/'figure_data/GWRF_annual_metrics_as_stored.csv')
    diag=stored.rename(columns={'Year':'year','MoranI_Train':'moran_i_train','p_val_Train':'p_value_train',
                                'MoranI_Test':'moran_i_test','p_val_Test':'p_value_test'})
    residual=pd.read_csv(ROOT.parent/'shap/Annual_Residual_Tests.csv').set_index('year')
    np.testing.assert_allclose(diag.set_index('year').loc[residual.index,residual.columns],residual)
    print('Annual_Residual_Tests: 21 archived annual diagnostics verified')


if __name__ == '__main__':
    main()
