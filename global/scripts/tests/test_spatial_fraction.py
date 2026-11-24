"""Executed in the existing WSL GIS/InVEST runtime."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from global_species import exact_fractions
from shapely.geometry import box

class FractionTests(unittest.TestCase):
    def test_fractional_edges_conserve_polygon_area(self):
        geom=box(250,250,1750,1750)
        weights=exact_fractions(geom,(0,1000,0,2000,0,-1000),(2,2));w=weights[weights>0]
        self.assertEqual(len(w),4)
        self.assertAlmostEqual(float(w.sum()),2.25,places=10)
        np.testing.assert_allclose(w,[.5625]*4,atol=1e-12)

    def test_subpixel_range_is_not_lost(self):
        geom=box(10,10,110,110)
        weights=exact_fractions(geom,(0,1000,0,2000,0,-1000),(2,2));w=weights[weights>0]
        self.assertEqual(w.sum(), np.float32(.01))  # production weights use float32

if __name__=='__main__':unittest.main()

