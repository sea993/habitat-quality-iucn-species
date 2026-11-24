"""Annual tiled InVEST runs. Actual InVEST 3.14.2, not a replacement model."""
import concurrent.futures
import csv
import json
import logging
import time
import zipfile
import traceback
import numpy as np
import natcap.invest
from natcap.invest import habitat_quality as hq
from global_common import *
from spatial_rules import distance_for_invest

def extract_year(year):
    existing=CFG.get('netcdf', {}).get(str(year))
    if existing and path(existing).exists():return path(existing),False
    matches=[]
    for archive in path(CFG['esa_archives']).glob('*.zip'):
        with zipfile.ZipFile(archive) as source:
            for item in source.infolist():
                if item.filename.endswith('.nc') and f'-{year}-' in item.filename:
                    matches.append((archive,item.filename,item.file_size))
    if len(matches)!=1:raise ValueError(f'Year {year}: {len(matches)} NetCDF matches')
    archive,member,size=matches[0]
    folder=OUT/'temp/netcdf';folder.mkdir(exist_ok=True)
    target=folder/Path(member).name
    if not target.exists() or target.stat().st_size!=size:
        space_check();event('extract_year',year=year,source=str(archive),bytes=size)
        partial=target.with_suffix('.partial')
        with zipfile.ZipFile(archive) as source,source.open(member) as src,partial.open('wb') as dst:
            shutil.copyfileobj(src,dst,16*1024**2)
        if partial.stat().st_size!=size:raise ValueError('Incomplete extraction')
        partial.replace(target)
    return target,True

def roads_for_tile(tile,local,bounds,buffer=None):
    folder=OUT/'inputs/roads';folder.mkdir(exist_ok=True)
    distance=CFG['buffer_m'] if buffer is None else buffer
    suffix='' if distance==CFG['buffer_m'] else f'_halo{distance}'
    target=folder/f'{tile["id"]:04d}{suffix}.tif'
    if target.exists():return target
    x0,y0,x1,y1=bounds
    nx,ny=round((x1-x0)/PIXEL),round((y1-y0)/PIXEL)
    ds=new_raster(target.with_suffix('.partial.tif'),nx,ny,gdal.GDT_Byte,255,(x0,PIXEL,0,y1,0,-PIXEL),local.ExportToWkt())
    ds.GetRasterBand(1).Fill(0)
    vector=ogr.Open(str(path(CFG['roads'])),0);layer=vector.GetLayer(0)
    tx=ogr.osr.CoordinateTransformation(WGS,local)
    memory=ogr.GetDriverByName('Memory').CreateDataSource('')
    lines=memory.CreateLayer('roads',srs=local,geom_type=ogr.wkbUnknown)
    count=0;seen=set()
    for bbox in buffered_wgs_boxes(tile,distance):
        layer.SetSpatialFilterRect(*bbox)
        for feature in layer:
            fid=feature.GetFID()
            if fid in seen:continue
            seen.add(fid)
            geom=feature.GetGeometryRef().Clone()
            if geom.IsEmpty():continue
            geom.Transform(tx)
            f=ogr.Feature(lines.GetLayerDefn());f.SetGeometry(geom);lines.CreateFeature(f);count+=1
    gdal.RasterizeLayer(ds,[1],lines,burn_values=[1],options=['ALL_TOUCHED=TRUE'])
    ds=None;vector=None;memory=None
    target.with_suffix('.partial.tif').replace(target)
    atomic_json(folder/f'{tile["id"]:04d}{suffix}.json',{'tile':tile['id'],'source':CFG['roads'],'features':count,
        'fixed_all_years':True,'binary_presence':True,'all_touched':True,'buffer_m':distance})
    return target

def run_local_model(tile,year,nc,folder,buffer=None):
    local,bounds=tile_geometry(tile,buffer)
    x0,y0,x1,y1=bounds
    folder.mkdir(parents=True,exist_ok=True)
    roads=roads_for_tile(tile,local,bounds,buffer)
    lulc=folder/'lulc.tif'
    src=f'NETCDF:"{nc}":lccs_class'
    gdal.Warp(str(lulc),src,format='GTiff',srcSRS='EPSG:4326',dstSRS=local.ExportToWkt(),
        outputBounds=bounds,xRes=PIXEL,yRes=PIXEL,resampleAlg='mode',errorThreshold=0.,srcNodata=0,dstNodata=0,
        outputType=gdal.GDT_Byte,creationOptions=['COMPRESS=DEFLATE','TILED=YES'],warpMemoryLimit=128*1024**2)
    ds=gdal.Open(str(lulc));a=ds.ReadAsArray();gt=ds.GetGeoTransform();wkt=ds.GetProjection();ds=None
    lookup=list(csv.DictReader((ROOT/'config/esa_landcover_crosswalk.csv').open(encoding='utf-8-sig')))
    allowed={int(r['esa_code']) for r in lookup}
    unknown=set(np.unique(a))-allowed-{0}
    if unknown:raise ValueError(f'Unmapped ESA classes: {unknown}')
    paths={'rr':str(roads)}
    crop=[int(r['esa_code']) for r in lookup if int(r['crop_threat'])]
    for threat,mask in [('urb',a==190),('crp',np.isin(a,crop))]:
        target=folder/(threat+'.tif');write_array(target,mask.astype(np.uint8),gt,wkt,255);paths[threat]=str(target)
    unit=str(hq.MODEL_SPEC['args']['threats_table_path']['columns']['max_dist']['units'])
    if natcap.invest.__version__!=CFG['invest_version']:raise ValueError('Unexpected InVEST version')
    rows=[{'threat':name,'max_dist':distance_for_invest(t['distance_m'],unit),'weight':t['weight'],
            'decay':t['decay'],'cur_path':paths[name]} for name,t in CFG['threats'].items()]
    write_csv(folder/'threats.csv',rows)
    args={'workspace_dir':str(folder/'model'),'results_suffix':str(year),'lulc_cur_path':str(lulc),
        'threats_table_path':str(folder/'threats.csv'),'sensitivity_table_path':str(ROOT/'config/sensitivity_invest3142.csv'),
        'half_saturation_constant':CFG['half_saturation_constant'],'n_workers':-1}
    issues=hq.validate(args)
    if issues:raise ValueError(issues)
    hq.execute(args)
    return folder/'model'/f'quality_c_{year}.tif',lulc,{'max_dist_unit':unit,'threats':rows,'args':args}

def worker_init():
    setup()
    logging.basicConfig(filename=OUT/'logs'/f'invest_worker_{os.getpid()}.log',level=logging.WARNING,
                        format='%(asctime)s %(name)s %(levelname)s %(message)s')

def sample_model_on_global_grid(quality,tile):
    """Nearest-neighbour reprojection of owned land-cell centers only.

    A rectangular polar output window also contains off-globe coordinates.
    Asking GDAL to estimate source bounds for those corners can skip valid
    land at the dateline. Explicit center transforms avoid those empty corners.
    """
    x,y,nx,ny=tile['window']
    owners=gdal.Open(str(OUT/'inputs/tile_id.tif')).ReadAsArray(x,y,nx,ny)==tile['id']
    rr,cc=np.where(owners)
    source=gdal.Open(str(quality));gt=source.GetGeoTransform();values=source.ReadAsArray()
    transform=Transformer.from_crs(WKT,source.GetProjection(),always_xy=True)
    xx,yy=transform.transform(XMIN+(x+cc+.5)*PIXEL,YMAX-(y+rr+.5)*PIXEL)
    finite=np.isfinite(xx)&np.isfinite(yy)
    col=np.zeros(len(rr),np.int64);row=np.zeros(len(rr),np.int64)
    col[finite]=np.floor((xx[finite]-gt[0])/gt[1]).astype(np.int64)
    row[finite]=np.floor((yy[finite]-gt[3])/gt[5]).astype(np.int64)
    inside=finite&(col>=0)&(row>=0)&(col<values.shape[1])&(row<values.shape[0])
    data=np.full((ny,nx),-1,np.float32)
    data[rr[inside],cc[inside]]=values[row[inside],col[inside]]
    source=None
    return data,owners

def process_tile(tile,year,nc):
    setup();space_check();started=time.monotonic()
    scratch=RAM_SCRATCH if CFG.get('scratch_backend')=='memory' else OUT/'temp'
    scratch.mkdir(parents=True,exist_ok=True)
    if scratch==RAM_SCRATCH and shutil.disk_usage(scratch).free<256*1024**2:
        raise RuntimeError('Insufficient free memory scratch space for the next tile')
    folder=scratch/f'{year}_{tile["id"]:04d}'
    # A prior failed attempt can leave InVEST taskgraph artifacts. Scratch only.
    remove_generated(folder)
    quality,lulc,pars=run_local_model(tile,year,nc,folder)
    x,y,nx,ny=tile['window'];gt=window_gt(x,y)
    data,owner=sample_model_on_global_grid(quality,tile)
    valid=owner&np.isfinite(data)&(data>=0)&(data<=1)
    if np.any(owner&np.isfinite(data)&(data>1+1e-6)):raise ValueError('HQ outside [0,1]')
    data[~valid]=-1
    target=OUT/'tiles'/f'{year}_{tile["id"]:04d}.tif'
    write_array(target.with_suffix('.partial.tif'),data,gt,WKT,-1)
    target.with_suffix('.partial.tif').replace(target)
    audit={'year':year,'tile':tile['id'],'land_cells':int(owner.sum()),'valid_cells':int(valid.sum()),
        'missing_land_cells':int((owner&~valid).sum()),'mean_hq':float(data[valid].mean()) if valid.any() else None,
        'elapsed_seconds':time.monotonic()-started,'invest_version':natcap.invest.__version__,
        'max_dist_unit':pars['max_dist_unit'],'buffer_m':CFG['buffer_m'],'reprojection':'nearest cell-center transform on owned land',
        'scratch_backend':CFG.get('scratch_backend','disk')}
    atomic_json(OUT/'progress'/f'{year}_{tile["id"]:04d}.json',audit)
    remove_generated(folder)
    return audit

def assemble(year,tiles):
    target=OUT/'rasters'/f'HQ_{year}.tif';partial=target.with_suffix('.partial.tif')
    if target.exists() and (OUT/'progress'/f'year_{year}.json').exists():return
    ds=new_raster(partial,dtype=gdal.GDT_Float32,nodata=-1);ds.GetRasterBand(1).Fill(-1)
    ownership=gdal.Open(str(OUT/'inputs/tile_id.tif'))
    for tile in tiles:
        x,y,nx,ny=tile['window'];src=OUT/'tiles'/f'{year}_{tile["id"]:04d}.tif'
        if not src.exists():raise FileNotFoundError(src)
        values=gdal.Open(str(src)).ReadAsArray();current=ds.ReadAsArray(x,y,nx,ny)
        use=ownership.ReadAsArray(x,y,nx,ny)==tile['id']
        current[use]=values[use];ds.GetRasterBand(1).WriteArray(current,x,y)
    ds=None;ownership=None;partial.replace(target)
    reports=[json.loads((OUT/'progress'/f'{year}_{tile["id"]:04d}.json').read_text()) for tile in tiles]
    atomic_json(OUT/'progress'/f'year_{year}.json',{'year':year,'status':'complete','tiles':len(tiles),
        'valid_cells':sum(r['valid_cells'] for r in reports),'missing_land_cells':sum(r['missing_land_cells'] for r in reports),
        'global_includes_china':True,'roads_fixed':True,'raster':str(target)})
    event('year_complete',year=year,tiles=len(tiles),raster=str(target))
    # Tile mosaics are reproducible intermediates. Remove only named run outputs
    # after the complete annual raster and its completion record both exist.
    for tile in tiles:
        src=(OUT/'tiles'/f'{year}_{tile["id"]:04d}.tif').resolve()
        if src.parent!=(OUT/'tiles').resolve():raise ValueError('Unexpected tile path')
        src.unlink(missing_ok=True)

def run_years(years=None):
    setup()
    tiles=json.loads((OUT/'inputs/tiles.json').read_text())
    priority={1249:0,801:1,877:2}
    ordered=sorted(tiles,key=lambda t:(priority.get(t['id'],3),-t['land_cells'],t['id']))
    for year in years or CFG['years']:
        if (OUT/'progress'/f'year_{year}.json').exists() and (OUT/'rasters'/f'HQ_{year}.tif').exists():continue
        nc,temporary=extract_year(year)
        pending=[t for t in ordered if not ((OUT/'progress'/f'{year}_{t["id"]:04d}.json').exists() and (OUT/'tiles'/f'{year}_{t["id"]:04d}.tif').exists())]
        event('year_start',year=year,pending_tiles=len(pending),all_tiles=len(tiles),workers=CFG['workers'])
        done=len(tiles)-len(pending)
        with concurrent.futures.ProcessPoolExecutor(max_workers=CFG['workers'],initializer=worker_init) as pool:
            jobs={pool.submit(process_tile,t,year,str(nc)):t['id'] for t in pending}
            for job in concurrent.futures.as_completed(jobs):
                try:
                    audit=job.result()
                except BaseException:
                    for queued in jobs:queued.cancel()
                    raise
                done+=1
                event('tile_complete',year=year,tile=audit['tile'],completed=done,total=len(tiles),seconds=round(audit['elapsed_seconds'],2),valid_cells=audit['valid_cells'])
                atomic_json(OUT/'progress/current.json',{'status':'running','stage':'annual_habitat_quality','year':year,
                    'completed_tiles_this_year':done,'total_tiles_per_year':len(tiles),'years':CFG['years'],'last_tile':audit['tile'],
                    'free_disk_gb':round(space_check(),2)})
        assemble(year,tiles)
        if temporary:
            nc=nc.resolve()
            if (OUT/'temp/netcdf').resolve() not in nc.parents:raise ValueError('Unexpected extraction path')
            nc.unlink()

if __name__=='__main__':
    run_years()
