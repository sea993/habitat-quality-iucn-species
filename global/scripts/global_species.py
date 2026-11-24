"""Complete species ranges, exact boundary fractions, annual/realm metrics.

Each species is dissolved before rasterization. SQLite commits a whole species
at a time, so a restart never duplicates observations or splits a trajectory.
"""
import csv
import json
import math
import sqlite3
import portalocker
import concurrent.futures
import multiprocessing
import numpy as np
import shapely
from shapely import box, make_valid, union_all, from_wkb, to_wkb, intersection, area
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.affinity import translate
from global_common import *
from global_background import realm_names,clean_realms
from spatial_rules import retain_species_record,weighted_quantile
from weighted_range_quantile import range_quantile

def polygon_parts(geom):
    if geom.geom_type=='Polygon':return [geom]
    if hasattr(geom,'geoms'):
        return [part for g in geom.geoms for part in polygon_parts(g)]
    return []

def split_dateline(geom):
    """Unwrap rings that cross +/-180, split, and return canonical longitudes."""
    pieces=[]
    for poly in polygon_parts(geom):
        exterior=np.asarray(poly.exterior.coords)
        if np.any(abs(np.diff(exterior[:,0]))>180) and poly.bounds[1]>-85 and poly.bounds[3]<85:
            ext=exterior.copy();ext[:,0]=np.rad2deg(np.unwrap(np.deg2rad(ext[:,0])))
            anchor=np.mean(ext[:,0]);holes=[]
            for ring in poly.interiors:
                coords=np.asarray(ring.coords).copy();coords[:,0]=np.rad2deg(np.unwrap(np.deg2rad(coords[:,0])))
                coords[:,0]+=360*round((anchor-np.mean(coords[:,0]))/360);holes.append(coords)
            poly=make_valid(Polygon(ext,holes))
            lo,_,hi,_=poly.bounds
            for shift in range(math.floor((lo+180)/360),math.floor((hi+180)/360)+1):
                clipped=poly.intersection(box(-180+shift*360,-90,180+shift*360,90))
                pieces.extend(polygon_parts(translate(clipped,xoff=-shift*360)))
        else:pieces.extend(polygon_parts(make_valid(poly)))
    return MultiPolygon(pieces) if pieces else MultiPolygon()

def prepare_species_parts():
    setup()
    with (OUT/'species/parts.lock').open('a+') as lock:
        portalocker.lock(lock,portalocker.LOCK_EX)
        return _prepare_species_parts()

def _prepare_species_parts():
    setup();folder=OUT/'species/parts';folder.mkdir(exist_ok=True)
    expected={int(r['species_id']) for r in csv.DictReader((ROOT/'tables/global_species_screened.csv').open(encoding='utf-8-sig'))}
    outputs=[]
    for srcpath in sorted(path(CFG['species_root']).glob('*/*.shp')):
        target=folder/(srcpath.stem+'.gpkg');marker=target.with_suffix('.json')
        if target.exists() and marker.exists():outputs.append(str(target));continue
        partial=target.with_suffix('.partial.gpkg')
        if partial.exists():partial.unlink()
        output=ogr.GetDriverByName('GPKG').CreateDataSource(str(partial))
        layerout=output.CreateLayer('parts',srs=GLOBAL_SRS,geom_type=ogr.wkbMultiPolygon)
        layerout.CreateField(ogr.FieldDefn('species_id',ogr.OFTInteger64))
        src=ogr.Open(str(srcpath),0);layer=src.GetLayer(0);tx=ogr.osr.CoordinateTransformation(WGS,GLOBAL_SRS)
        count=0;ids=set();repairs=0
        output.StartTransaction()
        for feature in layer:
            values=dict(feature.items())
            if not retain_species_record(values):continue
            key=int(values['id_no'])
            if key not in expected:raise ValueError('Species screen differs from archived species roster')
            geom=from_wkb(bytes(feature.GetGeometryRef().ExportToWkb()))
            if not geom.is_valid:repairs+=1
            geom=split_dateline(geom)
            if geom.is_empty:continue
            g=ogr.CreateGeometryFromWkb(to_wkb(geom));g.Segmentize(.05);g.Transform(tx)
            if not g.IsValid():g=g.MakeValid();repairs+=1
            clean=make_valid(from_wkb(bytes(g.ExportToWkb())))
            polys=polygon_parts(clean)
            if not polys:continue
            g=ogr.CreateGeometryFromWkb(to_wkb(MultiPolygon(polys)))
            feat=ogr.Feature(layerout.GetLayerDefn());feat.SetField('species_id',key);feat.SetGeometry(g);layerout.CreateFeature(feat)
            count+=1;ids.add(key)
            if count%1000==0:
                output.CommitTransaction();output.StartTransaction()
                event('species_part_progress',file=srcpath.name,retained_polygons=count)
        output.CommitTransaction();output.ExecuteSQL('CREATE INDEX idx_species ON parts (species_id)')
        output=None;src=None;partial.replace(target)
        atomic_json(marker,{'source':str(srcpath),'retained_polygon_records':count,'species':len(ids),
            'geometry_repairs':repairs,'dateline_policy':'unwrap crossing rings, split at +/-180 before projection','crs':CFG['global_crs']})
        outputs.append(str(target));event('species_part_complete',source=srcpath.name,polygons=count,species=len(ids))
    atomic_json(OUT/'species/parts_complete.json',{'files':outputs,'expected_species':len(expected)})
    return outputs

def burn_geometry(geom,gt,shape,boundary=False):
    memory=ogr.GetDriverByName('Memory').CreateDataSource('')
    layer=memory.CreateLayer('g',srs=GLOBAL_SRS,geom_type=ogr.wkbUnknown)
    feature=ogr.Feature(layer.GetLayerDefn());feature.SetGeometry(ogr.CreateGeometryFromWkb(to_wkb(geom)));layer.CreateFeature(feature)
    ds=gdal.GetDriverByName('MEM').Create('',shape[1],shape[0],1,gdal.GDT_Byte);ds.SetGeoTransform(gt);ds.SetProjection(WKT)
    gdal.RasterizeLayer(ds,[1],layer,burn_values=[1],options=['ALL_TOUCHED=TRUE'] if boundary else [])
    return ds.ReadAsArray()

def exact_fractions(geom,gt,shape):
    """Only boundary cells need polygon intersection; full interior cells weigh 1."""
    shapely.prepare(geom)
    filled=burn_geometry(geom,gt,shape)>0
    edge=burn_geometry(geom.boundary,gt,shape,boundary=True)>0
    fractions=filled.astype(np.float32)
    rr,cc=np.where(edge)
    for start in range(0,len(rr),20000):
        r,c=rr[start:start+20000],cc[start:start+20000]
        cells=box(gt[0]+c*PIXEL,gt[3]-(r+1)*PIXEL,gt[0]+(c+1)*PIXEL,gt[3]-r*PIXEL)
        covered=shapely.covers(geom,cells)
        fractions[r[covered],c[covered]]=1
        if not covered.all():
            use=~covered
            fractions[r[use],c[use]]=np.clip(area(intersection(cells[use],geom))/PIXEL**2,0,1)
    return fractions

def species_weights(geom,common,realms,china):
    shapely.prepare(geom)
    xmin,ymin,xmax,ymax=geom.bounds
    x0=max(0,int(math.floor((xmin-XMIN)/PIXEL)));x1=min(WIDTH,int(math.ceil((xmax-XMIN)/PIXEL)))
    y0=max(0,int(math.floor((YMAX-ymax)/PIXEL)));y1=min(HEIGHT,int(math.ceil((YMAX-ymin)/PIXEL)))
    chunks=[]
    for y in range(y0,y1,512):
        for x in range(x0,x1,512):
            nx,ny=min(512,x1-x),min(512,y1-y);gt=window_gt(x,y)
            footprint=box(gt[0],gt[3]-ny*PIXEL,gt[0]+nx*PIXEL,gt[3])
            if not geom.intersects(footprint):continue
            valid=common.ReadAsArray(x,y,nx,ny)>0
            if not valid.any():continue
            if geom.covers(footprint):
                fractions=np.ones((ny,nx),dtype=np.float32)
            else:
                piece=geom.intersection(footprint)
                if piece.is_empty or piece.area<=0:continue
                fractions=exact_fractions(piece,gt,(ny,nx))
            rr,cc=np.where(valid&(fractions>1e-10))
            if not len(rr):continue
            r=clean_realms(realms.ReadAsArray(x,y,nx,ny));c=china.ReadAsArray(x,y,nx,ny)
            chunks.append({'x':x,'y':y,'nx':nx,'ny':ny,'rr':rr.astype(np.uint16),'cc':cc.astype(np.uint16),
                'weights':fractions[rr,cc],'realm':r[rr,cc],'china':c[rr,cc]>0})
    return chunks

def prepare_hfp():
    setup()
    with (OUT/'inputs/hfp.lock').open('a+') as lock:
        portalocker.lock(lock,portalocker.LOCK_EX)
        _prepare_hfp()

def _prepare_hfp():
    for year in (2000,2020):
        target=OUT/'rasters'/f'HFP_{year}.tif'
        metadata_path=OUT/'inputs'/f'hfp_{year}_metadata.json'
        if target.exists() and metadata_path.exists():continue
        source=path(CFG['hfp_root'])/f'hfp{year}.tif'
        if not source.exists():raise FileNotFoundError(source)
        ds=gdal.Open(str(source))
        metadata={'source':str(source),'geotransform':ds.GetGeoTransform(),'crs':ds.GetProjection(),
            'nodata':ds.GetRasterBand(1).GetNoDataValue(),'scale':ds.GetRasterBand(1).GetScale(),'offset':ds.GetRasterBand(1).GetOffset()}
        if metadata['scale'] not in (None,1) or metadata['offset'] not in (None,0):raise ValueError('HFP scale/offset requires explicit application')
        if metadata['nodata'] is not None and not np.isfinite(metadata['nodata']):
            metadata['nodata']=str(metadata['nodata'])
        ds=None
        if not target.exists():
            partial=target.with_suffix('.partial.tif')
            if partial.exists():partial.unlink()
            gdal.Warp(str(target.with_suffix('.partial.tif')),str(source),dstSRS=WKT,outputBounds=CFG['global_bounds'],
                width=WIDTH,height=HEIGHT,resampleAlg='average',dstNodata=-9999,outputType=gdal.GDT_Float32,
                creationOptions=['COMPRESS=DEFLATE','TILED=YES','BIGTIFF=IF_SAFER','SPARSE_OK=YES'],warpMemoryLimit=256*1024**2)
            target.with_suffix('.partial.tif').replace(target)
        atomic_json(metadata_path,metadata)
        event('hfp_aligned',year=year)

def sample_chunks(ds,chunks):
    return np.concatenate([ds.ReadAsArray(c['x'],c['y'],c['nx'],c['ny'])[c['rr'],c['cc']] for c in chunks])

def extract_species(limit=None):
    setup()
    if (OUT/'progress/species_complete.json').exists():return
    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'species_preparation'})
    parts=prepare_species_parts();prepare_hfp()
    metadata=list(csv.DictReader((ROOT/'tables/global_species_screened.csv').open(encoding='utf-8-sig')))
    threshold=json.loads((OUT/'tables/global_low_quality_threshold.json').read_text())['threshold']
    db=sqlite3.connect(OUT/'species/metrics.sqlite')
    db.executescript('''
      CREATE TABLE IF NOT EXISTS metrics (species_id INTEGER, scientific_name TEXT, taxon TEXT, category TEXT, scope TEXT, year INTEGER,
        mean_hq REAL, q05_hq REAL, low_quality_pct REAL, low_quality_pct_0239 REAL, valid_area_km2 REAL, intersecting_cells INTEGER,
        PRIMARY KEY(species_id,scope,year));
      CREATE TABLE IF NOT EXISTS hfp (species_id INTEGER, year INTEGER, mean_hfp REAL, valid_area_km2 REAL, PRIMARY KEY(species_id,year));
      CREATE TABLE IF NOT EXISTS completed (species_id INTEGER PRIMARY KEY, status TEXT, source_area_km2 REAL, analyzed_area_km2 REAL);
    ''')
    done={row[0] for row in db.execute('SELECT species_id FROM completed')};processed=0
    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'species_extraction',
        'completed_species':len(done),'screened_species':len(metadata)})
    from global_species_worker import initialize,compute
    pending=[item for item in metadata if int(item['species_id']) not in done]
    if limit:pending=pending[:limit]
    workers=CFG.get('species_workers',2)
    iterator=iter(pending)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn'),
            initializer=initialize,initargs=(parts,threshold,CFG['years'],str(OUT/'inputs/common_valid_land.tif'))) as pool:
        jobs=set()
        for _ in range(min(len(pending),workers*2)):jobs.add(pool.submit(compute,next(iterator)))
        while jobs:
            finished,_=concurrent.futures.wait(jobs,return_when=concurrent.futures.FIRST_COMPLETED)
            for job in finished:
                jobs.remove(job)
                try:result=job.result()
                except BaseException:
                    for queued in jobs:queued.cancel()
                    raise
                key=result['key']
                with db:
                    db.executemany('INSERT INTO metrics VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',result['metrics'])
                    db.executemany('INSERT INTO hfp VALUES (?,?,?,?)',result['hfp'])
                    db.execute('INSERT INTO completed VALUES (?,?,?,?)',result['completed'])
                processed+=1;done.add(key)
                if processed%25==0:
                    event('species_metrics_progress',completed=len(done),total=len(metadata),last_species=key)
                    atomic_json(OUT/'progress/current.json',{'status':'running','stage':'species_extraction',
                        'completed_species':len(done),'screened_species':len(metadata),'workers':workers})
                item=next(iterator,None)
                if item is not None:jobs.add(pool.submit(compute,item))
    counts=dict(db.execute('SELECT status,COUNT(*) FROM completed GROUP BY status'))
    if sum(counts.values())==len(metadata):
        for table,name in [('metrics','species_annual_metrics.csv'),('hfp','species_hfp_endpoints.csv'),('completed','species_analysis_status.csv')]:
            cursor=db.execute('SELECT * FROM '+table+' ORDER BY species_id')
            target=OUT/'tables'/name
            with target.open('w',newline='',encoding='utf-8-sig') as stream:
                writer=csv.writer(stream);writer.writerow([d[0] for d in cursor.description]);writer.writerows(cursor)
        atomic_json(OUT/'progress/species_complete.json',{'status':'complete','screened_species':len(metadata),'analysis_counts':counts})
        event('species_extraction_complete',**counts)
    db.close()

if __name__=='__main__':extract_species()
