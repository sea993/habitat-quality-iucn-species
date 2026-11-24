"""Shared global grid, paths and resumable-run utilities. Sources are read only."""
import csv
import json
import math
import os
import shutil
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from osgeo import gdal, ogr, osr

gdal.UseExceptions()
ogr.UseExceptions()
ROOT = Path(__file__).resolve().parents[1]
CONFIG = Path(os.environ.get('HQ_CONFIG', ROOT / 'config/run_config.json')).resolve()
CFG = json.loads(CONFIG.read_text(encoding='utf-8-sig'))
OUT = Path(os.environ.get('HQ_OUTPUT', ROOT.parent.parent / 'review-output/global-upstream')).resolve()
RAM_SCRATCH = Path(os.environ.get('HQ_SCRATCH', OUT / 'temp/scratch')).resolve()
PIXEL = CFG['pixel_size_m']
XMIN, YMIN, XMAX, YMAX = CFG['global_bounds']
WIDTH, HEIGHT = round((XMAX-XMIN)/PIXEL), round((YMAX-YMIN)/PIXEL)
GT = (XMIN, PIXEL, 0, YMAX, 0, -PIXEL)

def path(value):
    return Path(os.environ.get('HQ_INPUTS', ROOT.parent.parent / 'external-inputs/global')).resolve() / value

def srs(value):
    result = osr.SpatialReference()
    result.SetFromUserInput(value)
    result.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return result

class Transformer:
    """Array wrapper around the already-installed GDAL/PROJ transformer."""
    @classmethod
    def from_crs(cls, source, target, always_xy=True):
        obj=cls()
        obj.tx=osr.CoordinateTransformation(srs(source),srs(target))
        return obj

    def transform(self,x,y):
        x,y=np.broadcast_arrays(np.asarray(x,dtype=float),np.asarray(y,dtype=float))
        shape=x.shape
        xx=np.empty(x.size);yy=np.empty(y.size)
        for start in range(0,x.size,100000):
            stop=min(x.size,start+100000)
            points=np.column_stack((x.ravel()[start:stop],y.ravel()[start:stop]))
            transformed=np.asarray(self.tx.TransformPoints(points))
            xx[start:stop]=transformed[:,0];yy[start:stop]=transformed[:,1]
        return xx.reshape(shape),yy.reshape(shape)

GLOBAL_SRS = srs(CFG['global_crs'])
WGS = srs('EPSG:4326')
WKT = GLOBAL_SRS.ExportToWkt()
TO_WGS = Transformer.from_crs(WKT, 'EPSG:4326', always_xy=True)

def setup():
    for name in ('inputs', 'rasters', 'tiles', 'temp', 'logs', 'progress', 'tables', 'figures', 'species'):
        (OUT/name).mkdir(parents=True, exist_ok=True)
    gdal.SetCacheMax(CFG['gdal_cache_mb_per_worker']*1024**2)

def atomic_json(target, value):
    target = Path(target)
    temp = target.with_suffix(target.suffix+'.partial')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(target)

def write_csv(target, rows):
    rows = list(rows)
    if not rows:
        return
    target = Path(target)
    temp = target.with_suffix(target.suffix+'.partial')
    with temp.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temp.replace(target)

def event(stage, **kwargs):
    record = {'time_utc': datetime.now(timezone.utc).isoformat(), 'stage': stage, **kwargs}
    print(json.dumps(record, ensure_ascii=False), flush=True)
    with (OUT/'logs/events.jsonl').open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(record, ensure_ascii=False)+'\n')
    return record

def space_check():
    available = shutil.disk_usage(OUT).free/1e9
    if available < CFG['minimum_free_disk_gb']:
        raise RuntimeError(f'Only {available:.1f} GB free; computation paused before writing new data')
    return available

def remove_generated(folder):
    """Only named scratch descendants of this project may be removed."""
    folder = Path(folder).resolve()
    allowed = [(OUT/'temp').resolve(),RAM_SCRATCH.resolve()]
    if folder in allowed or not any(root in folder.parents for root in allowed):
        raise ValueError(f'Unsafe scratch cleanup target: {folder}')
    if folder.exists():
        shutil.rmtree(folder)

def new_raster(target, width=WIDTH, height=HEIGHT, dtype=gdal.GDT_Byte, nodata=0, gt=GT, wkt=WKT):
    ds = gdal.GetDriverByName('GTiff').Create(str(target), width, height, 1, dtype,
        options=['TILED=YES', 'COMPRESS=DEFLATE', 'ZLEVEL=4', 'BIGTIFF=IF_SAFER', 'SPARSE_OK=YES', 'NUM_THREADS=1'])
    ds.SetGeoTransform(gt)
    ds.SetProjection(wkt)
    ds.GetRasterBand(1).SetNoDataValue(nodata)
    return ds

def write_array(target, array, gt, wkt, nodata, dtype=None):
    if dtype is None:
        dtype = gdal.GDT_Byte if array.dtype == np.uint8 else gdal.GDT_Float32
    ds = new_raster(target, array.shape[1], array.shape[0], dtype, nodata, gt, wkt)
    ds.GetRasterBand(1).WriteArray(array)
    ds = None

def window_gt(x, y):
    return (XMIN+x*PIXEL, PIXEL, 0, YMAX-y*PIXEL, 0, -PIXEL)

def tile_id(lon, lat):
    # Longitude 180 belongs to the last interval, not to an extra tile.
    col = np.clip(np.floor((np.asarray(lon)+180)/6).astype(int), 0, 59)
    row = np.clip(np.floor((np.asarray(lat)+90)/6).astype(int), 0, 29)
    return (row*60+col+1).astype(np.uint16)

def tile_geometry(tile, buffer=None):
    """Local equal-area metric domain enclosing a 6-degree core plus halo."""
    buffer = CFG['buffer_m'] if buffer is None else buffer
    west, south, east, north = tile['wgs_bounds']
    lon, lat = (west+east)/2, (south+north)/2
    local = srs(f'+proj=laea +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m +no_defs')
    transform = Transformer.from_crs('EPSG:4326', local.ExportToWkt(), always_xy=True)
    # Dense perimeter, plus center, captures curved edges near the poles.
    lons = np.r_[np.linspace(west,east,121), np.linspace(west,east,121), np.full(121,west), np.full(121,east)]
    lats = np.r_[np.full(121,south), np.full(121,north), np.linspace(south,north,121), np.linspace(south,north,121)]
    xx, yy = transform.transform(lons,lats)
    bounds = [math.floor((min(xx)-buffer)/PIXEL)*PIXEL, math.floor((min(yy)-buffer)/PIXEL)*PIXEL,
              math.ceil((max(xx)+buffer)/PIXEL)*PIXEL, math.ceil((max(yy)+buffer)/PIXEL)*PIXEL]
    return local, bounds

def buffered_wgs_boxes(tile, buffer=None):
    """Conservative lon/lat search bounds, split explicitly at the dateline."""
    buffer = CFG['buffer_m'] if buffer is None else buffer
    west,south,east,north = tile['wgs_bounds']
    # 100 km/degree is deliberately below the true meridional scale.
    dy = buffer/100000 + .05
    south,north = max(-90,south-dy),min(90,north+dy)
    extreme = max(abs(south),abs(north))
    if extreme>89:
        return [(-180,south,180,north)]
    dx = dy/max(math.cos(math.radians(extreme)),.01)
    west,east = west-dx,east+dx
    if west < -180:
        return [(-180,south,east,north),(west+360,south,180,north)]
    if east > 180:
        return [(west,south,180,north),(-180,south,east-360,north)]
    return [(west,south,east,north)]
