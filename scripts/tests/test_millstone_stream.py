"""Integration boundaries and atomic checkpoint readiness for live native fits."""
import datetime as dt, importlib.util, json, pathlib, sys, tempfile, unittest
import h5py
import numpy as np
SCRIPTS=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SCRIPTS))
spec=importlib.util.spec_from_file_location('millstone_stream',SCRIPTS/'stream-millstone-direct.py')
stream=importlib.util.module_from_spec(spec);spec.loader.exec_module(stream)
class IntegrationWindows(unittest.TestCase):
    def test_legacy_anchor_strict_upper_bound_and_unknown_window(self):
        start=int(dt.datetime(2024,4,8,tzinfo=dt.timezone.utc).timestamp())
        with tempfile.TemporaryDirectory() as tmp:
            day=pathlib.Path(tmp);(day/'audit').mkdir()
            available=list(range(10,210,10))+list(range(400,610,10))+list(range(800,1010,10))
            for t in range(0,1100,10):
                with h5py.File(day/'audit'/f'{start+t}.h5','w') as h:
                    h['outputs_json']=np.frombuffer(json.dumps([{'path':'reduced'}] if t in available else []).encode(),np.uint8)
            groups=stream.partitions('2024-04-08',day)
            # Independent legacy loop: the first available source anchors each
            # interval; starts at or beyond +300 s begin a new integration.
            expected=[];current=[];anchor=available[0]
            for t in available:
                if t>=anchor+300:expected.append(current);current=[];anchor=t
                current.append(t)
            if anchor+300<=1100:expected.append(current)
            actual=[[t-start for t in group if t-start in available] for group in groups]
            self.assertEqual(actual,expected)
            self.assertEqual(groups[0][0],start+10)
            self.assertEqual(groups[0][-1],start+300)
            (day/'audit'/f'{start+500}.h5').unlink()
            self.assertEqual(len(stream.partitions('2024-04-08',day)),1)
    def test_incomplete_integration_is_not_ready(self):
        start=int(dt.datetime(2024,4,8,tzinfo=dt.timezone.utc).timestamp())
        with tempfile.TemporaryDirectory() as tmp:
            day=pathlib.Path(tmp);(day/'audit').mkdir()
            for t in range(0,290,10):
                with h5py.File(day/'audit'/f'{start+t}.h5','w') as h:
                    h['outputs_json']=np.frombuffer(b'[{"path":"reduced"}]',np.uint8)
            self.assertEqual(stream.partitions('2024-04-08',day),[])
if __name__=='__main__':unittest.main()
