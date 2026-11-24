import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from global_common import *

class GlobalRules(unittest.TestCase):
    def test_grid_includes_china_and_dateline(self):
        self.assertEqual(int(tile_id(112.5,30.5)),1249)
        self.assertNotEqual(int(tile_id(-179.9,0)),int(tile_id(179.9,0)))
        self.assertEqual(int(tile_id(180,90)),1800)

    def test_halo_dateline_splits(self):
        boxes=buffered_wgs_boxes({'wgs_bounds':[174,0,180,6]})
        self.assertEqual(len(boxes),2)
        self.assertEqual(boxes[1][0],-180)
        self.assertGreater(boxes[1][2],-180)

    def test_local_buffer_and_grid(self):
        tile={'wgs_bounds':[108,30,114,36]}
        local,core=tile_geometry(tile,0)
        _,halo=tile_geometry(tile,15000)
        self.assertEqual(core[0]-halo[0],15000)
        self.assertEqual(halo[2]-core[2],15000)
        self.assertEqual(core[1]-halo[1],15000)
        self.assertEqual(halo[3]-core[3],15000)
        self.assertTrue(all(x%1000==0 for x in halo))

    def test_safe_cleanup_rejects_sources(self):
        with self.assertRaises(ValueError):remove_generated(path(CFG['esa_archives']))
        with self.assertRaises(ValueError):remove_generated(OUT/'temp')

if __name__=='__main__':unittest.main()
