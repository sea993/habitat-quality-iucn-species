"""Create fixed land/realm masks and tile ownership on the global 1 km grid."""
import collections
import json
import numpy as np
from osgeo import gdal, ogr
from global_common import *

def without_holes(geometry):
    """Keep inland water in the land footprint; ESA supplies its water class."""
    flat = ogr.GT_Flatten(geometry.GetGeometryType())
    if flat == ogr.wkbPolygon:
        result=ogr.Geometry(ogr.wkbPolygon)
        if geometry.GetGeometryCount(): result.AddGeometry(geometry.GetGeometryRef(0))
        return result
    result=ogr.Geometry(ogr.wkbMultiPolygon)
    for index in range(geometry.GetGeometryCount()):
        item=without_holes(geometry.GetGeometryRef(index))
        if ogr.GT_Flatten(item.GetGeometryType())==ogr.wkbPolygon:result.AddGeometry(item)
        else:
            for j in range(item.GetGeometryCount()):result.AddGeometry(item.GetGeometryRef(j))
    return result

def project_boundaries():
    target=OUT/'inputs/boundaries.gpkg'
    if target.exists() and (OUT/'inputs/boundaries_complete.json').exists():return target
    if target.exists():target.unlink()
    out=ogr.GetDriverByName('GPKG').CreateDataSource(str(target))
    land=out.CreateLayer('land',srs=GLOBAL_SRS,geom_type=ogr.wkbMultiPolygon)
    land.CreateField(ogr.FieldDefn('china',ogr.OFTInteger))
    src=ogr.Open(str(path(CFG['land_boundaries'])),0);layer=src.GetLayer(0)
    tx=ogr.osr.CoordinateTransformation(srs('EPSG:4326'),GLOBAL_SRS)
    country_names=[]
    out.StartTransaction()
    for feature in layer:
        name=feature.GetField('NAME_0');country_names.append(name)
        geom=without_holes(feature.GetGeometryRef())
        # Preserve the coastline while densifying very long straight segments.
        geom.Segmentize(.1)
        geom.Transform(tx)
        if not geom.IsValid():geom=geom.MakeValid()
        feat=ogr.Feature(land.GetLayerDefn())
        feat.SetField('china',int(name in {'China','Taiwan','Hong Kong','Macao','Macau'}))
        feat.SetGeometry(ogr.ForceToMultiPolygon(geom));land.CreateFeature(feat)
    out.CommitTransaction();src=None
    realms=out.CreateLayer('realms',srs=GLOBAL_SRS,geom_type=ogr.wkbMultiPolygon)
    realms.CreateField(ogr.FieldDefn('realm_id',ogr.OFTInteger))
    src=ogr.Open(str(path(CFG['ecoregions'])),0);layer=src.GetLayer(0)
    names=sorted({feature.GetField('REALM') for feature in layer if feature.GetField('REALM')})
    mapping={name:i+1 for i,name in enumerate(names)};layer.ResetReading()
    out.StartTransaction()
    for feature in layer:
        geom=feature.GetGeometryRef().Clone();geom.Segmentize(.1);geom.Transform(tx)
        if not geom.IsValid():geom=geom.MakeValid()
        feat=ogr.Feature(realms.GetLayerDefn());feat.SetField('realm_id',mapping.get(feature.GetField('REALM'),0))
        feat.SetGeometry(ogr.ForceToMultiPolygon(geom));realms.CreateFeature(feat)
    out.CommitTransaction();out=None;src=None
    atomic_json(OUT/'inputs/realm_names.json',{'0':'Unclassified',**{str(v):k for k,v in mapping.items()}})
    atomic_json(OUT/'inputs/boundaries_complete.json',{'countries':len(country_names),'china_names':[n for n in country_names if n in {'China','Taiwan','Hong Kong','Macao','Macau'}],
        'land_source':CFG['land_boundaries'],'land_policy':'country exterior rings, inland water retained; center-based raster mask',
        'realm_source':CFG['ecoregions'],'realm_count':len(mapping),'crs':CFG['global_crs']})
    return target

def prepare():
    setup()
    if (OUT/'inputs/grid_complete.json').exists():return
    event('prepare_boundaries')
    source=project_boundaries()
    for filename,layername,attribute,where in [('land.tif','land',None,None),('china.tif','land','china',None),('realm.tif','realms','realm_id',None)]:
        ds=new_raster(OUT/'inputs'/filename)
        vector=ogr.Open(str(source));layer=vector.GetLayerByName(layername)
        options=['ATTRIBUTE='+attribute] if attribute else []
        gdal.RasterizeLayer(ds,[1],layer,burn_values=[1] if not attribute else [],options=options)
        ds=None;vector=None
        event('mask_complete',mask=filename)
    land=gdal.Open(str(OUT/'inputs/land.tif'))
    owners=new_raster(OUT/'inputs/tile_id.tif',dtype=gdal.GDT_UInt16)
    info={};total=0
    for y in range(0,HEIGHT,128):
        a=land.ReadAsArray(0,y,WIDTH,min(128,HEIGHT-y))
        rr,cc=np.where(a>0)
        ids=np.zeros(a.shape,np.uint16)
        if len(rr):
            lon,lat=TO_WGS.transform(XMIN+(cc+.5)*PIXEL,YMAX-(rr+y+.5)*PIXEL)
            if not np.isfinite(lon).all() or not np.isfinite(lat).all():raise ValueError('Land outside valid global projection')
            tid=tile_id(lon,lat);ids[rr,cc]=tid;total+=len(rr)
            for key in np.unique(tid):
                select=tid==key;x0=int(cc[select].min());x1=int(cc[select].max())+1;y0=int(rr[select].min())+y;y1=int(rr[select].max())+y+1
                key=int(key)
                if key not in info:info[key]=[x0,y0,x1,y1,0]
                old=info[key];old[:4]=[min(old[0],x0),min(old[1],y0),max(old[2],x1),max(old[3],y1)];old[4]+=int(select.sum())
        owners.GetRasterBand(1).WriteArray(ids,0,y)
    owners=None;land=None
    tiles=[]
    for key,(x0,y0,x1,y1,n) in sorted(info.items()):
        row,col=divmod(key-1,60)
        tiles.append({'id':key,'wgs_bounds':[-180+6*col,-90+6*row,-174+6*col,-84+6*row],
                      'window':[x0,y0,x1-x0,y1-y0],'land_cells':n})
    atomic_json(OUT/'inputs/tiles.json',tiles)
    atomic_json(OUT/'inputs/grid_complete.json',{'grid':{'crs':CFG['global_crs'],'geotransform':GT,'width':WIDTH,'height':HEIGHT},
                'land_cells':total,'land_area_km2_at_1km_grid':total,'tiles':len(tiles),'global_includes_china':True})
    event('global_grid_complete',tiles=len(tiles),land_cells=total)

if __name__=='__main__':prepare()
