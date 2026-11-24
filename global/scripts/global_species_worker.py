"""Read-only per-species extraction workers; the parent owns all SQLite writes."""
import numpy as np
from global_species import (setup,ogr,gdal,OUT,from_wkb,make_valid,union_all,
    species_weights,sample_chunks,realm_names,range_quantile)

CONTEXT=None

def initialize(parts,threshold,years,support_path):
    global CONTEXT
    setup()
    CONTEXT={'sources':[ogr.Open(p,0) for p in parts],
        'threshold':threshold,'names':realm_names(),
        'common':gdal.Open(str(support_path)),
        'realms':gdal.Open(str(OUT/'inputs/realm.tif')),
        'china':gdal.Open(str(OUT/'inputs/china.tif')),
        'hqs':{year:gdal.Open(str(OUT/'rasters'/f'HQ_{year}.tif')) for year in years},
        'hfps':{year:gdal.Open(str(OUT/'rasters'/f'HFP_{year}.tif')) for year in (2000,2020)}}

def compute(item):
    context=CONTEXT
    if context is None:raise RuntimeError('Species worker has not been initialized')
    key=int(item['species_id']);pieces=[]
    for source in context['sources']:
        layer=source.GetLayerByName('parts');layer.SetAttributeFilter(f'species_id={key}')
        for feature in layer:pieces.append(from_wkb(bytes(feature.GetGeometryRef().ExportToWkb())))
    if not pieces:
        return {'key':key,'metrics':[],'hfp':[],'completed':(key,'no_retained_polygon',0,0)}
    geom=make_valid(union_all(pieces))
    chunks=species_weights(geom,context['common'],context['realms'],context['china'])
    if not chunks:
        return {'key':key,'metrics':[],'hfp':[],'completed':(key,'no_common_valid_land',geom.area/1e6,0)}
    weights=np.concatenate([c['weights'] for c in chunks]).astype(np.float64)
    r=np.concatenate([c['realm'] for c in chunks]);c=np.concatenate([c['china'] for c in chunks])
    groups={'Global':np.ones(len(weights),bool)}
    if c.any():groups['China']=c
    for code in np.unique(r):groups['Realm:'+context['names'][int(code)]]=r==code
    records=[]
    for year,ds in context['hqs'].items():
        values=sample_chunks(ds,chunks)
        if not np.all(np.isfinite(values)&(values>=0)&(values<=1)):raise ValueError('Invalid value on common support')
        for scope,select in groups.items():
            w=weights[select];v=values[select];total=w.sum()
            records.append((key,item['scientific_name'],item['taxon'],item['category'],scope,year,
                float(np.average(v,weights=w)),range_quantile(v,w,.05),float(100*w[v<context['threshold']].sum()/total),
                float(100*w[v<=.239].sum()/total),float(total),len(w)))
    hfp_rows=[]
    for year,ds in context['hfps'].items():
        values=sample_chunks(ds,chunks);valid=np.isfinite(values)&(values!=-9999)
        hfp_rows.append((key,year,float(np.average(values[valid],weights=weights[valid])) if valid.any() else None,float(weights[valid].sum())))
    return {'key':key,'metrics':records,'hfp':hfp_rows,'completed':(key,'analyzed',geom.area/1e6,float(weights.sum()))}
