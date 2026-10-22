import pathlib,sys,unittest
import numpy as np
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from millstone_archive import recorded_power

class MeterCadence(unittest.TestCase):
 def test_missing_per_second_files_preserve_all_five_seconds(self):
  got=[recorded_power(t,[0,5],[400e3,600e3])for t in range(6)]
  np.testing.assert_allclose(got,[400e3,440e3,480e3,520e3,560e3,600e3])
 def test_gaps_and_invalid_series_remain_unavailable(self):
  self.assertTrue(np.isnan(recorded_power(100,[],[])))
  self.assertTrue(np.isnan(recorded_power(100,[0,200],[1e6,1e6])))
  with self.assertRaises(ValueError):recorded_power(0,[5,0],[1.,1.])
if __name__=='__main__':unittest.main()
