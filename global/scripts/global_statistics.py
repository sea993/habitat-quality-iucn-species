"""Species-equal summaries and 2,000 whole-trajectory bootstrap replicates."""
import csv
import json
import sqlite3
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from global_common import *

METRICS=['mean_hq','q05_hq','low_quality_pct','low_quality_pct_0239']
CHINA_ORIGINAL = ROOT.parent / 'china/statistics/species_habitat_summaries.csv'
CHINA_IDS = ROOT.parent / 'china/statistics/species.csv'

def ols_slopes(values,years):
    x=np.array(years,dtype=float,copy=True);x-=x.mean()
    values=np.asarray(values,dtype=float)
    centered=values-values[:,0:1,...]
    return np.tensordot(centered,x,axes=([1],[0]))/np.dot(x,x)

def bootstrap_means(values,replicates,seed):
    """Sample species, keeping all years and all metrics coupled within a draw."""
    values=np.asarray(values,dtype=float)
    rng=np.random.default_rng(seed);n=len(values)
    results=np.empty((replicates,)+values.shape[1:])
    for i in range(replicates):
        selection=rng.integers(0,n,n)
        results[i]=values[selection].mean(axis=0)
    return results

def intervals(values):return np.quantile(values,[.025,.975],axis=0)

def calculate():
    setup()
    if (OUT/'progress/statistics_complete.json').exists():return
    if not (OUT/'progress/species_complete.json').exists():raise RuntimeError('Species extraction incomplete')
    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'statistics_and_confidence_intervals'})
    db=sqlite3.connect(OUT/'species/metrics.sqlite')
    df=pd.read_sql_query('SELECT * FROM metrics',db)
    background=pd.read_csv(OUT/'tables/background_annual.csv')
    years=CFG['years'];x=np.array(years,dtype=float);x-=x.mean();denom=x@x
    original_ids=set(pd.read_csv(CHINA_IDS).id_no.astype(int))
    if len(original_ids)!=1820:raise ValueError('Original China species roster is not the manuscript 1,820-species sample')
    effects=[];trajectories=[];species_slopes=[];counter=0
    for scope in sorted(df.scope.unique()):
        sub=df[df.scope==scope].copy()
        if scope=='China':sub=sub[sub.species_id.isin(original_ids)]
        if sub.empty:continue
        info=sub.drop_duplicates('species_id').set_index('species_id')
        ids=np.array(sorted(info.index))
        data=np.stack([sub.pivot(index='species_id',columns='year',values=m).reindex(index=ids,columns=years).to_numpy() for m in METRICS],axis=2)
        if not np.isfinite(data).all():raise ValueError('Incomplete species-year trajectory')
        slopes=ols_slopes(data,years)
        for i,key in enumerate(ids):
            species_slopes.append({'species_id':int(key),'scope':scope,'taxon':info.loc[key,'taxon'],'category':info.loc[key,'category'],
                **{m+'_slope':float(slopes[i,j]) for j,m in enumerate(METRICS)},
                **{m+'_change_2000_2020':float(data[i,-1,j]-data[i,0,j]) for j,m in enumerate(METRICS)},
                'q05_minus_mean_slope':float(slopes[i,1]-slopes[i,0]),
                'q05_declines_faster':bool(slopes[i,1]<min(0,slopes[i,0]))})
        labels={'All species':np.ones(len(ids),bool),'Threatened species':info.loc[ids,'category'].isin(['VU','EN','CR']).to_numpy()}
        for cat in ['LC','NT','VU','EN','CR']:labels['Category:'+cat]=(info.loc[ids,'category']==cat).to_numpy()
        for taxon in sorted(info.taxon.unique()):labels['Taxon:'+taxon]=(info.loc[ids,'taxon']==taxon).to_numpy()
        bg=background[background.scope==scope].set_index('year').reindex(years)
        bg_slope=float((bg.mean_hq.to_numpy()-bg.mean_hq.iloc[0])@x/denom) if len(bg.dropna())==len(years) else np.nan
        for label,use in labels.items():
            n=int(use.sum())
            if not n:continue
            counter+=1;sample=data[use];point=sample.mean(axis=0)
            boot=bootstrap_means(sample,CFG['bootstrap_replicates'],CFG['random_seed']+counter)
            ci=intervals(boot);bs=ols_slopes(boot,years);sl=slopes[use].mean(axis=0);sci=intervals(bs)
            annual_change=point-point[0]
            annual_change_ci=intervals(boot-boot[:,0:1,:])
            change=point[-1]-point[0];bchange=boot[:,-1]-boot[:,0];dci=intervals(bchange)
            for j,m in enumerate(METRICS):
                for yi,year in enumerate(years):
                    trajectories.append({'scope':scope,'group':label,'n_species':n,'year':year,'metric':m,
                        'estimate':point[yi,j],'ci95_low':ci[0,yi,j] if n>1 else None,'ci95_high':ci[1,yi,j] if n>1 else None,
                        'change_from_2000':annual_change[yi,j],
                        'change_ci95_low':annual_change_ci[0,yi,j] if n>1 else None,
                        'change_ci95_high':annual_change_ci[1,yi,j] if n>1 else None})
                effects.append({'scope':scope,'group':label,'n_species':n,'effect':m+'_annual_OLS_slope','estimate':sl[j],
                    'ci95_low':sci[0,j] if n>1 else None,'ci95_high':sci[1,j] if n>1 else None,'units':'HQ/year' if j<2 else 'percentage points/year'})
                effects.append({'scope':scope,'group':label,'n_species':n,'effect':m+'_change_2000_2020','estimate':change[j],
                    'ci95_low':dci[0,j] if n>1 else None,'ci95_high':dci[1,j] if n>1 else None,'units':'HQ' if j<2 else 'percentage points'})
            difference=bs[:,1]-bs[:,0];dci=intervals(difference)
            effects.append({'scope':scope,'group':label,'n_species':n,'effect':'q05_minus_mean_paired_slope','estimate':sl[1]-sl[0],
                'ci95_low':dci[0] if n>1 else None,'ci95_high':dci[1] if n>1 else None,'units':'HQ/year'})
            faster=(slopes[use,1]<np.minimum(0,slopes[use,0])).astype(float)
            pci=intervals(bootstrap_means(faster,CFG['bootstrap_replicates'],CFG['random_seed']+counter))
            effects.append({'scope':scope,'group':label,'n_species':n,'effect':'fraction_species_q05_declines_faster','estimate':faster.mean(),
                'ci95_low':pci[0] if n>1 else None,'ci95_high':pci[1] if n>1 else None,'units':'proportion'})
            if np.isfinite(bg_slope):
                for j,m in enumerate(METRICS[:2]):
                    diff=bs[:,j]-bg_slope;ci_diff=intervals(diff)
                    effects.append({'scope':scope,'group':label,'n_species':n,'effect':m+'_slope_minus_background_mean_slope','estimate':sl[j]-bg_slope,
                        'ci95_low':ci_diff[0] if n>1 else None,'ci95_high':ci_diff[1] if n>1 else None,'units':'HQ/year'})
                    if bg_slope< -1e-6 and sl[j]<0:
                        ratio_ci=intervals(bs[:,j]/bg_slope)
                        effects.append({'scope':scope,'group':label,'n_species':n,'effect':m+'_decline_ratio_to_background_mean','estimate':sl[j]/bg_slope,
                            'ci95_low':ratio_ci[0] if n>1 else None,'ci95_high':ratio_ci[1] if n>1 else None,'units':'ratio'})
                level_difference=boot[:,:,0].mean(axis=1)-bg.mean_hq.mean();lci=intervals(level_difference)
                effects.append({'scope':scope,'group':label,'n_species':n,'effect':'period_mean_hq_minus_background','estimate':point[:,0].mean()-bg.mean_hq.mean(),
                    'ci95_low':lci[0] if n>1 else None,'ci95_high':lci[1] if n>1 else None,'units':'HQ'})
            event('statistics_group_complete',scope=scope,group=label,n_species=n)
            atomic_json(OUT/'progress/current.json',{'status':'running','stage':'statistics_and_confidence_intervals',
                'completed_groups':counter,'last_scope':scope,'last_group':label})
    write_csv(OUT/'tables/group_annual_trajectories.csv',trajectories)
    write_csv(OUT/'tables/group_effects_ci.csv',effects)
    write_csv(OUT/'tables/species_slopes.csv',species_slopes)
    slopes_df=pd.DataFrame(species_slopes)
    original=pd.read_csv(CHINA_ORIGINAL);original.id_no=original.id_no.astype(int)
    original=original[original.id_no.isin(original_ids)]
    china=slopes_df[slopes_df.scope=='China'].merge(original,left_on='species_id',right_on='id_no',how='inner',suffixes=('_globalflow','_original'))
    china.to_csv(OUT/'tables/china_matched_species_comparison.csv',index=False,encoding='utf-8-sig')
    atomic_json(OUT/'progress/core_statistics_complete.json',{'status':'complete','effect_rows':len(effects),
        'bootstrap_replicates':CFG['bootstrap_replicates'],'hfp_auxiliary':'pending common-support calculation'})
    from global_hfp_support import harmonize
    harmonize()
    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'hfp_auxiliary_statistics'})
    hfp=pd.read_sql_query('SELECT * FROM hfp',db).pivot(index='species_id',columns='year',values='mean_hfp')
    joined=slopes_df[slopes_df.scope=='Global'].set_index('species_id').join(hfp)
    joined['hfp_change']=joined[2020]-joined[2000]
    joined.to_csv(OUT/'tables/hfp_hq_species_comparison.csv',encoding='utf-8-sig')
    auxiliary=[]
    for m in ['mean_hq_change_2000_2020','q05_hq_change_2000_2020']:
        pair=joined[['hfp_change',m]].dropna()
        result=spearmanr(pair.hfp_change,pair[m])
        rng=np.random.default_rng(CFG['random_seed']);values=pair.to_numpy();boot=[]
        if len(pair)>2:
            for i in range(CFG['bootstrap_replicates']):
                sample=values[rng.integers(0,len(values),len(values))]
                boot.append(spearmanr(sample[:,0],sample[:,1]).statistic)
        finite=np.array(boot)[np.isfinite(boot)];ci=np.quantile(finite,[.025,.975]) if len(finite) else [np.nan,np.nan]
        auxiliary.append({'comparison':'HFP change vs '+m,'n_species':len(pair),'spearman_rho':result.statistic,
            'ci95_low':ci[0],'ci95_high':ci[1],'p_value':result.pvalue,'interpretation':'association, not independent causal attribution'})
    write_csv(OUT/'tables/hfp_auxiliary_correlations.csv',auxiliary)
    atomic_json(OUT/'progress/statistics_complete.json',{'status':'complete','bootstrap_replicates':CFG['bootstrap_replicates'],
        'seed':CFG['random_seed'],'resampling_unit':'species; full 21-year trajectory and paired metrics kept together',
        'confidence_interval':'percentile 95%; species sampling uncertainty, not model/input uncertainty',
        'china_original_sample':len(original_ids),'china_matched_species':len(china),'effect_rows':len(effects),
        'ratio_background_threshold':'background mean-HQ slope < -0.000001 HQ/year and species-group slope < 0',
        'q05_faster_definition':'q05 slope < mean slope and q05 slope < 0'})
    db.close()

if __name__=='__main__':calculate()
