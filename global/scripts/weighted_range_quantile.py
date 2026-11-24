"""Exact weighted quantiles for pixel fractions, avoiding a full large sort."""
import math
import numpy as np
from spatial_rules import weighted_quantile

def range_quantile(values,weights,q):
    v=np.asarray(values)
    w=np.asarray(weights,dtype=float)
    valid=np.isfinite(v)&np.isfinite(w)&(w>0)
    if not valid.all():
        v=v[valid];w=w[valid]
    n=len(v)
    if n<50000 or not 0<q<1 or np.max(w)>1:
        return weighted_quantile(v,w,q)
    total=float(w.sum());target=q*total
    # Since each pixel weight is <=1, the crossing rank is bounded by
    # target and target plus the total fractional area absent from the pixels.
    # Two extra ranks keep integer and floating-point boundaries conservative.
    lo=max(0,math.floor(target)-2)
    hi=min(n-1,math.ceil(target+n-total)+2)
    order=np.argpartition(v,(lo,hi))
    preceding=float(w[order[:lo]].sum())
    middle=order[lo:hi+1]
    middle=middle[np.argsort(v[middle])]
    cumulative=np.cumsum(w[middle])
    index=int(np.searchsorted(cumulative,target-preceding,side='left'))
    if index>=len(middle):
        return weighted_quantile(v,w,q)
    return float(v[middle[index]])
