import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from global_hfp_support import adjust_paired_means

class HfpSupportTests(unittest.TestCase):
    def test_fractional_weights_match_direct_common_cells(self):
        original={2000:(4/1.5,1.5),2020:(7.6/1.2,1.2)}
        excluded={2000:(.5,2.),2020:(.2,1.6)}
        result,error=adjust_paired_means(original,excluded)
        self.assertAlmostEqual(result[2000][0],2.)
        self.assertAlmostEqual(result[2020][0],6.)
        self.assertEqual(result[2000][1],result[2020][1])
        self.assertAlmostEqual(error,0.)

    def test_identical_support_keeps_original_means(self):
        original={2000:(1.234567,12.5),2020:(3.456789,12.5)}
        result,error=adjust_paired_means(original,{2000:(0,0),2020:(0,0)})
        self.assertEqual(result,original);self.assertEqual(error,0)

    def test_disjoint_support_has_no_mean(self):
        result,error=adjust_paired_means({2000:(2.,1.),2020:(3.,2.)},{2000:(1.,2.),2020:(2.,6.)})
        self.assertEqual(result,{2000:(None,0.),2020:(None,0.)})

    def test_inconsistent_area_is_rejected(self):
        with self.assertRaises(ValueError):
            adjust_paired_means({2000:(2.,1.),2020:(3.,2.)},{2000:(0,0),2020:(0,0)})

if __name__=='__main__':unittest.main()
