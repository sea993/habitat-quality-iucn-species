import sys
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from spatial_rules import distance_for_invest, retain_species_record, weighted_quantile, mixture_parameters

class PilotCoreTests(unittest.TestCase):
    def test_distance_units(self):
        self.assertEqual(distance_for_invest(10000, 'kilometer'), 10)
        self.assertEqual(distance_for_invest(10000, 'meter'), 10000)
        with self.assertRaises(ValueError): distance_for_invest(10000, 'degree')

    def test_species_rules_keep_china(self):
        r=dict(id_no=1, category='VU', presence=1, origin=1, seasonal=1, terrestria='true', country='China')
        self.assertTrue(retain_species_record(r))
        for key,value in [('presence',4),('origin',3),('seasonal',4),('category','DD'),('terrestria','false')]:
            self.assertFalse(retain_species_record(dict(r, **{key:value})))
        self.assertTrue(retain_species_record(dict(r, origin=2, seasonal=3)))

    def test_weighted_quantile_area(self):
        self.assertAlmostEqual(weighted_quantile([.1,.9],[1,99],.05),.9)
        self.assertAlmostEqual(weighted_quantile([.1,.9],[10,90],.05),.1)

    def test_mixture(self):
        base={1:dict(habitat=.5,urb=.7,crp=.6,rr=.4),2:dict(habitat=1.,urb=.8,crp=.6,rr=.5)}
        self.assertAlmostEqual(mixture_parameters({1:.75,2:.25},base)['habitat'],.625)
        with self.assertRaises(ValueError):mixture_parameters({1:.5},base)

if __name__=='__main__':unittest.main()
