import sys
import unittest
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon,box
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from global_species import split_dateline,exact_fractions
from global_background import jenks_histogram
from global_statistics import ols_slopes,bootstrap_means

class GlobalAnalysisTests(unittest.TestCase):
    def test_constant_trajectories_have_exactly_zero_slope(self):
        years=np.arange(2000,2021,dtype=float)
        before=years.copy()
        constants=np.array([[.1,.78316572],[.999999,.333333333333]])
        values=np.repeat(constants[:,None,:],21,axis=1)
        self.assertTrue(np.array_equal(ols_slopes(values,years),np.zeros_like(constants)))
        self.assertTrue(np.array_equal(years,before))

    def test_dateline_range_not_transworld(self):
        polygon=Polygon([(179,-1),(-179,-1),(-179,1),(179,1),(179,-1)])
        parts=split_dateline(polygon)
        self.assertAlmostEqual(parts.area,4)
        self.assertEqual(len(parts.geoms),2)
        self.assertFalse(parts.intersects(box(0,-1,1,1)))

    def test_fractional_area_and_tiny_species(self):
        gt=(0,1000,0,2000,0,-1000)
        for geom in [box(100,100,1500,1800),box(100,100,200,200),box(0,0,2000,2000).difference(box(300,300,500,500))]:
            fractions=exact_fractions(geom,gt,(2,2))
            self.assertAlmostEqual(float(fractions.sum(dtype=float)),geom.area/1e6,places=6)

    def test_jenks_five_known_clusters(self):
        hist=np.zeros(100)
        hist[[1,20,40,60,95]]=[5,10,20,5,10]
        self.assertEqual(jenks_histogram(hist,5),[.02,.21,.41,.61,.96])

    def test_paired_bootstrap_keeps_trajectories(self):
        years=np.arange(2000,2021)
        values=np.empty((4,21,2))
        for i in range(4):
            values[i,:,0]=i+.01*(years-2000)
            values[i,:,1]=2*values[i,:,0]-.2
        slopes=ols_slopes(values,years)
        self.assertTrue(np.allclose(slopes,np.array([[.01,.02]]*4)))
        boot=bootstrap_means(values,100,20261117)
        self.assertTrue(np.allclose(boot[:,:,1],2*boot[:,:,0]-.2))
        self.assertTrue(np.allclose(np.diff(boot[:,:,0],axis=1),.01))

if __name__=='__main__':unittest.main()
