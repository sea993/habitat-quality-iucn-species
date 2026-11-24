import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from spatial_rules import weighted_quantile
from weighted_range_quantile import range_quantile

class RangeQuantileTests(unittest.TestCase):
    def test_large_fractional_pixels_equal_full_sort(self):
        rng=np.random.default_rng(20261117)
        for fractional in [.001,.1,1.0]:
            values=rng.random(60001).astype(np.float32)
            weights=np.ones(len(values))
            weights[:int(len(values)*fractional)]=rng.random(int(len(values)*fractional))
            for q in [.01,.05,.5,.95]:
                self.assertEqual(range_quantile(values,weights,q),weighted_quantile(values,weights,q))

    def test_tied_values_and_integer_rank(self):
        values=np.repeat(np.arange(10,dtype=float),6000)/9
        weights=np.ones(len(values))
        for q in [.05,.1,.5,.9]:
            self.assertEqual(range_quantile(values,weights,q),weighted_quantile(values,weights,q))
        values=values[::-1];weights[::100]=.25
        self.assertEqual(range_quantile(values,weights,.05),weighted_quantile(values,weights,.05))

    def test_nonfinite_and_nonfractional_weights(self):
        values=np.arange(60003,dtype=float);weights=np.ones(len(values))
        values[0]=np.nan;weights[1]=np.nan;weights[2]=0;weights[3]=3
        self.assertEqual(range_quantile(values,weights,.05),weighted_quantile(values,weights,.05))

if __name__=='__main__':unittest.main()
