"""Actual model verification before permitting a long production run."""
import json
import logging
import time
import numpy as np
from global_common import *
from global_rasters import run_local_model,extract_year

def validate():
    setup()
    logging.basicConfig(filename=OUT/'logs/global_validation_model.log',level=logging.WARNING)
    output=OUT/'progress/model_validation.json'
    if output.exists() and json.loads(output.read_text()).get('status')=='passed':return
    original=json.loads((ROOT/'config/model_validation_reference.json').read_text())
    if original['status']!='passed':raise AssertionError('Original numerical model test failed')
    pilot={'threats': {'urb': {'distance_m': 10000, 'weight': 1.0, 'decay': 'exponential'}, 'crp': {'distance_m': 2000, 'weight': 0.25, 'decay': 'linear'}, 'rr': {'distance_m': 3000, 'weight': 1.0, 'decay': 'linear'}}, 'half_saturation_constant': 0.5}
    if pilot['threats']!=CFG['threats'] or pilot['half_saturation_constant']!=CFG['half_saturation_constant']:
        raise AssertionError('Scientific parameters diverge from pilot')
    tile={'id':1249,'wgs_bounds':[108,30,114,36]}
    nc,_=extract_year(2000)
    arrays=[]
    _,core=tile_geometry(tile,0)
    for halo in (15000,25000):
        folder=OUT/'temp'/f'validation_halo{halo}'
        remove_generated(folder)
        event('validate_buffer_start',halo_m=halo)
        q,l,pars=run_local_model(tile,2000,str(nc),folder,buffer=halo)
        ds=gdal.Open(str(q));gt=ds.GetGeoTransform()
        x=round((core[0]-gt[0])/PIXEL);y=round((gt[3]-core[3])/PIXEL)
        nx=round((core[2]-core[0])/PIXEL);ny=round((core[3]-core[1])/PIXEL)
        arrays.append(ds.ReadAsArray(x,y,nx,ny));ds=None
    common=(arrays[0]>=0)&(arrays[1]>=0)
    difference=abs(arrays[0][common]-arrays[1][common])
    maximum=float(difference.max())
    if maximum>1e-6:raise AssertionError(f'Buffer difference {maximum}')
    if not np.array_equal(arrays[0]>=0,arrays[1]>=0):raise AssertionError('Buffer valid masks differ')
    report={'status':'passed','original_scientific_tests':original,'scientific_parameters_match_pilot':True,
            'full_tile_buffer_comparison_m':[15000,25000],'tested_valid_cells':int(common.sum()),
            'maximum_absolute_hq_difference':maximum,'model':'natcap.invest.habitat_quality 3.14.2',
            'global_includes_china':True,'static_grip4':True}
    atomic_json(output,report)
    for halo in (15000,25000):remove_generated(OUT/'temp'/f'validation_halo{halo}')
    event('global_model_validation_passed',**report)

if __name__=='__main__':validate()
