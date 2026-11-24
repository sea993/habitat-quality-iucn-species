"""Common 21-year support, pooled weighted Jenks breaks and backgrounds."""
import json
from functools import lru_cache
import numpy as np
from global_common import *

def jenks_histogram(counts,classes=5):
    """Weighted Fisher-Jenks optimum for a fixed fine histogram, no sampling.

    Every global land-cell/year contributes once. Value bin width <= 0.0001;
    returned boundaries are upper bin edges and are reused for all years.
    """
    counts=np.asarray(counts,dtype=np.float64)
    occupied=np.flatnonzero(counts)
    if len(occupied)<classes:raise ValueError('Fewer distinct HQ bins than classes')
    values=(occupied+.5)/len(counts);weights=counts[occupied];n=len(values)
    sw=np.r_[0,np.cumsum(weights)];sx=np.r_[0,np.cumsum(weights*values)];sx2=np.r_[0,np.cumsum(weights*values**2)]
    old=np.full(n+1,np.inf);old[0]=0
    trace=np.zeros((classes+1,n+1),np.int32)
    for k in range(1,classes+1):
        new=np.full(n+1,np.inf)
        for end in range(k,n+1):
            starts=np.arange(k-1,end)
            w=sw[end]-sw[starts];x=sx[end]-sx[starts]
            cost=np.maximum(0,sx2[end]-sx2[starts]-x*x/w)
            candidates=old[starts]+cost
            best=int(np.argmin(candidates));new[end]=candidates[best];trace[k,end]=starts[best]
        old=new
    ends=[];end=n
    for k in range(classes,0,-1):
        ends.append(end);end=int(trace[k,end])
    ends.reverse()
    return [float((occupied[e-1]+1)/len(counts)) for e in ends]

@lru_cache(maxsize=1)
def realm_names():
    names=json.loads((OUT/'inputs/realm_names.json').read_text())
    return {int(k):v for k,v in names.items() if v not in ('N/A','NA','')}

def clean_realms(a):
    valid=realm_names()
    out=a.copy();out[~np.isin(out,list(valid))]=0
    return out

def build_background():
    setup()
    if (OUT/'progress/background_complete.json').exists():return
    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'common_support_and_threshold'})
    for year in CFG['years']:
        if not (OUT/'progress'/f'year_{year}.json').exists():raise RuntimeError(f'Annual raster incomplete: {year}')
    years=CFG['years'];names=realm_names();scopes=['Global','China']+[f'Realm:{v}' for v in names.values()]
    bins=CFG['jenks_histogram_bins'];hist=np.zeros((len(scopes),len(years),bins),np.uint64)
    sums=np.zeros((len(scopes),len(years)),np.float64)
    counts=np.zeros(len(scopes),np.uint64);original_low=np.zeros_like(sums)
    land=gdal.Open(str(OUT/'inputs/land.tif'));realm=gdal.Open(str(OUT/'inputs/realm.tif'));china=gdal.Open(str(OUT/'inputs/china.tif'))
    rasters=[gdal.Open(str(OUT/'rasters'/f'HQ_{year}.tif')) for year in years]
    common=new_raster(OUT/'inputs/common_valid_land.partial.tif')
    processed=0;initial_land=0
    for y in range(0,HEIGHT,512):
        for x in range(0,WIDTH,512):
            nx,ny=min(512,WIDTH-x),min(512,HEIGHT-y)
            mask=land.ReadAsArray(x,y,nx,ny)>0
            if not mask.any():continue
            initial_land+=int(mask.sum())
            arrays=[d.ReadAsArray(x,y,nx,ny) for d in rasters]
            for a in arrays:mask &= np.isfinite(a)&(a>=0)&(a<=1)
            common.GetRasterBand(1).WriteArray(mask.astype(np.uint8),x,y)
            r=clean_realms(realm.ReadAsArray(x,y,nx,ny));c=china.ReadAsArray(x,y,nx,ny)>0
            selections=[mask,mask&c]+[mask&(r==key) for key in names]
            for group,use in enumerate(selections):
                n=int(use.sum());counts[group]+=n
                if not n:continue
                for yi,a in enumerate(arrays):
                    vals=a[use];indices=np.minimum((vals.astype(float)*bins).astype(np.int32),bins-1)
                    hist[group,yi]+=np.bincount(indices,minlength=bins).astype(np.uint64)
                    sums[group,yi]+=vals.sum(dtype=np.float64)
                    original_low[group,yi]+=(vals<=CFG['comparison_low_quality_threshold']).sum()
            processed+=1
        event('common_support_progress',row=y,total_rows=HEIGHT,land_blocks=processed)
        atomic_json(OUT/'progress/current.json',{'status':'running','stage':'common_support_and_threshold',
            'completed_grid_rows':min(y+512,HEIGHT),'total_grid_rows':HEIGHT,'land_blocks':processed})
    common=None
    (OUT/'inputs/common_valid_land.partial.tif').replace(OUT/'inputs/common_valid_land.tif')
    if not counts[0]:raise ValueError('No globally common valid land')
    np.savez_compressed(OUT/'tables/background_histograms.npz',hist=hist,sums=sums,counts=counts,scopes=np.array(scopes),years=np.array(years))
    breaks=jenks_histogram(hist[0].sum(axis=0),CFG['jenks_classes']);threshold=breaks[0]
    rows=[]
    for group,scope in enumerate(scopes):
        n=int(counts[group])
        if not n:continue
        for yi,year in enumerate(years):
            cumulative=np.cumsum(hist[group,yi]);q05=(np.searchsorted(cumulative,.05*n)+1)/bins
            low_count=hist[group,yi,:round(threshold*bins)].sum()
            rows.append({'scope':scope,'year':year,'area_km2':n,'grid_cells':n,'mean_hq':sums[group,yi]/n,
                'q05_hq_histogram_upper_edge':q05,'low_quality_pct':100*low_count/n,
                'low_quality_pct_threshold_0239':100*original_low[group,yi]/n,'low_quality_threshold':threshold})
    write_csv(OUT/'tables/background_annual.csv',rows)
    atomic_json(OUT/'tables/global_low_quality_threshold.json',{'threshold':threshold,'jenks_five_class_upper_bounds':breaks,
        'method':'weighted Fisher-Jenks on pooled global cell-year histogram, not a random subsample',
        'histogram_bin_width':1/bins,'years':years,'common_land_cells':int(counts[0]),
        'common_land_area_km2':int(counts[0]),'initial_land_cells':initial_land,
        'land_cells_excluded_for_missing_any_year':initial_land-int(counts[0]),'comparison_threshold':.239,
        'threshold_rule':'HQ strictly below first class upper bin edge; all years fixed'})
    atomic_json(OUT/'progress/background_complete.json',{'status':'complete','threshold':threshold,'common_land_cells':int(counts[0])})
    event('background_complete',threshold=threshold,common_land_cells=int(counts[0]))

if __name__=='__main__':build_background()
