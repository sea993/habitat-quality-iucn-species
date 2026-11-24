"""Small, tested scientific rules for land-cover parameters and species screening."""
import numpy as np

FIELDS = ('habitat','urb','crp','rr')
ALLOWED_CATEGORIES = {'LC','NT','VU','EN','CR'}

def distance_for_invest(meters, unit):
    if unit in ('kilometer','kilometers','km'): return meters / 1000.
    if unit in ('meter','meters','m'): return float(meters)
    raise ValueError(f'Unsupported distance unit: {unit}')

def exclusion_reason(r):
    if str(r.get('terrestria','')).lower() not in ('true','1'):return 'not_terrestrial'
    if r.get('category') not in ALLOWED_CATEGORIES:return 'category_not_LC_NT_VU_EN_CR'
    if r.get('presence') != 1:return 'presence_not_extant'
    if r.get('origin') not in (1,2):return 'origin_not_native_or_reintroduced'
    if r.get('seasonal') not in (1,2,3):return 'season_not_resident_breeding_nonbreeding'
    if r.get('id_no') is None:return 'missing_species_id'
    return None

def retain_species_record(r):return exclusion_reason(r) is None

def mixture_parameters(weights, base):
    if not np.isclose(sum(weights.values()),1.) or min(weights.values())<0:
        raise ValueError('Mixture weights must be nonnegative and sum to one')
    return {f:sum(w*base[c][f] for c,w in weights.items()) for f in FIELDS}

def weighted_quantile(values, weights, q):
    v=np.asarray(values,dtype=float);w=np.asarray(weights,dtype=float)
    keep=np.isfinite(v)&np.isfinite(w)&(w>0)
    v=v[keep];w=w[keep]
    if not len(v):return float('nan')
    idx=np.argsort(v);v=v[idx];w=w[idx]
    return float(v[min(np.searchsorted(np.cumsum(w),q*w.sum(),side='left'),len(v)-1)])
