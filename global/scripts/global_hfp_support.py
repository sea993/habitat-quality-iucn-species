"""Put both HFP endpoints on identical species-range support.

Only pixels valid in one endpoint are re-extracted. Their weighted sums and
areas are subtracted from the original endpoint totals; HQ is never recomputed.
"""
import concurrent.futures
import csv
import portalocker
import multiprocessing
import sqlite3
import sys
import numpy as np
from global_common import *
import global_species_worker as worker
from global_species import from_wkb,make_valid,union_all,species_weights,sample_chunks

def build_support_difference():
    setup()
    target=OUT/'inputs/hfp_support_difference.tif'
    marker=OUT/'inputs/hfp_support_difference.json'
    with (OUT/'inputs/hfp_support_difference.lock').open('a+') as lock:
        portalocker.lock(lock,portalocker.LOCK_EX)
        if target.exists() and marker.exists():return target
        datasets=[gdal.Open(str(OUT/'rasters'/f'HFP_{year}.tif')) for year in (2000,2020)]
        land=gdal.Open(str(OUT/'inputs/common_valid_land.tif'))
        partial=target.with_suffix('.partial.tif')
        if partial.exists():partial.unlink()
        output=new_raster(partial)
        counts={'both_valid_cells':0,'only_2000_cells':0,'only_2020_cells':0,'neither_valid_cells':0}
        for y in range(0,HEIGHT,512):
            for x in range(0,WIDTH,512):
                nx,ny=min(512,WIDTH-x),min(512,HEIGHT-y)
                use=land.ReadAsArray(x,y,nx,ny)>0
                if not use.any():continue
                values=[ds.ReadAsArray(x,y,nx,ny) for ds in datasets]
                a,b=[np.isfinite(v)&(v!=-9999)&use for v in values]
                output.GetRasterBand(1).WriteArray((a^b).astype(np.uint8),x,y)
                counts['both_valid_cells']+=int((a&b).sum())
                counts['only_2000_cells']+=int((a&~b).sum())
                counts['only_2020_cells']+=int((b&~a).sum())
                counts['neither_valid_cells']+=int((use&~a&~b).sum())
        output=None;partial.replace(target)
        atomic_json(marker,{**counts,'scope':'common HQ land mask','endpoint_years':[2000,2020]})
        event('hfp_support_difference_ready',**counts)
        return target

def excluded_contributions(item):
    context=worker.CONTEXT;key=int(item['species_id']);pieces=[]
    for source in context['sources']:
        layer=source.GetLayerByName('parts');layer.SetAttributeFilter(f'species_id={key}')
        for feature in layer:pieces.append(from_wkb(bytes(feature.GetGeometryRef().ExportToWkb())))
    geom=make_valid(union_all(pieces))
    chunks=species_weights(geom,context['common'],context['realms'],context['china'])
    result={year:(0.,0.) for year in (2000,2020)}
    if chunks:
        weights=np.concatenate([c['weights'] for c in chunks]).astype(float)
        validity=[]
        for year,ds in context['hfps'].items():
            values=sample_chunks(ds,chunks);valid=np.isfinite(values)&(values!=-9999)
            validity.append(valid)
            result[year]=(float(weights[valid].sum()),float(np.sum(weights[valid]*values[valid],dtype=np.float64)))
        if not np.all(validity[0]^validity[1]):raise ValueError('HFP difference mask does not match endpoint validity')
    return key,result

def adjust_paired_means(original,excluded):
    """original: year -> (mean, area); excluded: year -> (area, weighted sum)."""
    if all(excluded[y][0]==0 for y in (2000,2020)):
        if abs(original[2000][1]-original[2020][1])>1e-6:raise ValueError('Unequal areas without exclusive cells')
        return original.copy(),0.
    remaining={y:original[y][1]-excluded[y][0] for y in (2000,2020)}
    discrepancy=abs(remaining[2000]-remaining[2020])
    if discrepancy>1e-6 or min(remaining.values())< -1e-6:
        raise ValueError(f'Common HFP areas disagree: {remaining}')
    area=remaining[2000] if excluded[2000][0]==0 else remaining[2020] if excluded[2020][0]==0 else (remaining[2000]+remaining[2020])/2
    area=max(0.,area)
    result={}
    for year in (2000,2020):
        mean,original_area=original[year]
        if original_area>0 and (mean is None or not np.isfinite(mean)):raise ValueError('Invalid original HFP mean')
        total=(mean*original_area if original_area>0 else 0.)-excluded[year][1]
        result[year]=(mean if excluded[year][0]==0 and area>0 else float(total/area) if area>0 else None,area)
    return result,discrepancy

def harmonize():
    setup()
    marker=OUT/'progress/hfp_common_support_complete.json'
    if marker.exists():return
    if not (OUT/'progress/species_complete.json').exists():raise RuntimeError('Wait for the complete HQ/species extraction')
    difference=build_support_difference()
    with (OUT/'progress/hfp_harmonize.lock').open('a+') as lock:
        portalocker.lock(lock,portalocker.LOCK_EX)
        if marker.exists():return
        db=sqlite3.connect(OUT/'species/metrics.sqlite')
        db.executescript('''CREATE TABLE IF NOT EXISTS hfp_original_support AS SELECT * FROM hfp;
            CREATE UNIQUE INDEX IF NOT EXISTS idx_original_hfp ON hfp_original_support(species_id,year);
            CREATE TABLE IF NOT EXISTS hfp_harmonized (species_id INTEGER PRIMARY KEY, excluded_area_2000 REAL,
                excluded_area_2020 REAL, area_discrepancy REAL);''')
        support=json.loads((OUT/'inputs/hfp_support_difference.json').read_text())
        if support['only_2000_cells']==0:
            # The raster scan proves 2000 support is nested in 2020 support.
            # Equal endpoint areas therefore need no polygon re-extraction.
            with db:
                db.execute('''INSERT OR IGNORE INTO hfp_harmonized
                    SELECT a.species_id,0,0,0 FROM hfp_original_support a
                    JOIN hfp_original_support b ON a.species_id=b.species_id
                    WHERE a.year=2000 AND b.year=2020 AND a.valid_area_km2=b.valid_area_km2''')
                db.execute('''UPDATE hfp SET mean_hfp=NULL,valid_area_km2=0
                    WHERE species_id IN (SELECT species_id FROM hfp_original_support WHERE year=2000 AND valid_area_km2=0)''')
                db.execute('''INSERT OR IGNORE INTO hfp_harmonized
                    SELECT a.species_id,0,b.valid_area_km2,0 FROM hfp_original_support a
                    JOIN hfp_original_support b ON a.species_id=b.species_id
                    WHERE a.year=2000 AND b.year=2020 AND a.valid_area_km2=0''')
        ids={r[0] for r in db.execute('SELECT DISTINCT species_id FROM hfp_original_support')}
        done={r[0] for r in db.execute('SELECT species_id FROM hfp_harmonized')}
        pending_ids=ids-done
        items=[r for r in csv.DictReader((ROOT/'tables/global_species_screened.csv').open(encoding='utf-8-sig'))
               if int(r['species_id']) in pending_ids]
        parts=json.loads((OUT/'species/parts_complete.json').read_text())['files']
        workers=CFG.get('species_workers',8)
        atomic_json(OUT/'progress/current.json',{'status':'running','stage':'hfp_common_support',
            'completed_species':len(done),'total_species':len(ids)})
        iterator=iter(items)
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
                initializer=worker.initialize,initargs=(parts,0.,[],str(difference))) as pool:
            jobs={pool.submit(excluded_contributions,next(iterator)) for _ in range(min(len(items),workers*2))}
            while jobs:
                finished,_=concurrent.futures.wait(jobs,return_when=concurrent.futures.FIRST_COMPLETED)
                for job in finished:
                    jobs.remove(job)
                    try:key,excluded=job.result()
                    except BaseException:
                        for queued in jobs:queued.cancel()
                        raise
                    original={y:(mean,area) for y,mean,area in db.execute(
                        'SELECT year,mean_hfp,valid_area_km2 FROM hfp_original_support WHERE species_id=?',(key,))}
                    corrected,discrepancy=adjust_paired_means(original,excluded)
                    with db:
                        db.executemany('UPDATE hfp SET mean_hfp=?,valid_area_km2=? WHERE species_id=? AND year=?',
                            [(mean,area,key,y) for y,(mean,area) in corrected.items()])
                        db.execute('INSERT INTO hfp_harmonized VALUES (?,?,?,?)',
                            (key,excluded[2000][0],excluded[2020][0],discrepancy))
                    done.add(key)
                    if len(done)%100==0:
                        event('hfp_common_support_progress',completed=len(done),total=len(ids))
                        atomic_json(OUT/'progress/current.json',{'status':'running','stage':'hfp_common_support',
                            'completed_species':len(done),'total_species':len(ids)})
                    item=next(iterator,None)
                    if item is not None:jobs.add(pool.submit(excluded_contributions,item))
        if len(done)!=len(ids):raise ValueError('Incomplete HFP support harmonization')
        cursor=db.execute('SELECT * FROM hfp ORDER BY species_id,year')
        with (OUT/'tables/species_hfp_endpoints.csv').open('w',newline='',encoding='utf-8-sig') as stream:
            writer=csv.writer(stream);writer.writerow([d[0] for d in cursor.description]);writer.writerows(cursor)
        affected,max_error=db.execute('SELECT SUM(excluded_area_2000>0 OR excluded_area_2020>0),MAX(area_discrepancy) FROM hfp_harmonized').fetchone()
        db.close()
        atomic_json(marker,{'status':'complete','species':len(ids),'species_with_endpoint_exclusions':affected,
            'max_common_area_discrepancy_km2':max_error,'method':'Same spatial support at both HFP endpoints; subtract exclusive-cell weighted totals',
            'original_endpoint_results':'hfp_original_support table in metrics.sqlite','hq_metrics_changed':False})
        event('hfp_common_support_complete',species=len(ids),affected=affected)

if __name__=='__main__':
    if '--mask-only' in sys.argv:build_support_difference()
    else:harmonize()
