"""Numerical I/O checks for the native-fit browser bridge."""
import importlib.util, json, pathlib, tempfile, unittest
import h5py
import numpy as np

script=pathlib.Path(__file__).resolve().parents[1]/'export-millstone-native-view.py'
spec=importlib.util.spec_from_file_location('native_export',script)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def store(group,name,value):
    group[name]=np.frombuffer(json.dumps(value).encode(),dtype=np.uint8)

class NativeView(unittest.TestCase):
    def test_absolute_density_error_units_and_axis_packing(self):
        with tempfile.TemporaryDirectory() as tmp:
            public=pathlib.Path(tmp);historical=public/'old';historical.mkdir()
            (historical/'manifest.json').write_text(json.dumps(dict(antenna='ZENITH',date='2024-04-08',raw_url='raw.f32',plasma_url='/guisdap/data/plasma.f32',calibration=dict(gain=3e8))))
            fits=public/'fits';fits.mkdir()
            source=public/'source.h5'
            with h5py.File(source,'w') as h:
                h.create_group('records/a/ion_line')['system_temperature_k']=[170.]
            expected=[]
            for i in range(2):
                parameters=np.array([[1e11*(i+1),1000+i,1800+i,-30-i],[2e11*(i+1),1100+i,1900+i,40+i]])
                expected.append(parameters)
                with h5py.File(fits/f'fit-{i:06}.h5','w') as h:
                    r=h.create_group('records/a');r['range_m']=[300e3,400e3]
                    store(r,'configuration_json',dict(antenna='ZENITH',start_unix_s=1712599200+i*180,end_unix_s=1712599380+i*180))
                    fit=r.create_group('fit');fit['parameters']=parameters;fit['standard_error']=parameters*.01
                    ion=r.create_group('ion_line');ion['lag_s']=[1e-5,2e-5];ion['acf']=np.arange(8.).reshape(2,2,2)*1e-6;ion['variance']=np.full((2,2),2e-14)
                    rs=[]
                    for rg in [300e3,400e3]:
                        rs.append(dict(range_m=rg,result=dict(status='converged',reduced_chi_square=1.2,model_curve=[dict(lag_s=1e-5,value=[2e-6,3e-6]),dict(lag_s=2e-5,value=[4e-6,5e-6])]),error=None))
                    native=r.create_group('ion_fit');store(native,'results_json',rs);store(native,'assumptions_json',dict(calibration=6e11,calibration_units='W m^-3 K^-1 km^-1'))
                    sources=h.create_group('sources/combined');store(sources,'records_json',[dict(path=str(source),start_unix_s=0,end_unix_s=1,transmitter_w=1e6,azimuth_deg=0,elevation_deg=90)])
            output=public/'new'
            module.export(historical,fits,output,public)
            packed=np.fromfile(output/'parameters.f32','<f4').reshape(4,2,2)
            m=json.loads((output/'manifest.json').read_text())
            for i in range(2):
                np.testing.assert_array_equal(packed[:,:,i].T,expected[i].astype('f4'))
                np.testing.assert_array_equal(m['fits'][i]['parameters'],expected[i])
                np.testing.assert_array_equal(m['fits'][i]['standard_error'],expected[i]*.01)
            d=json.loads((output/'fit-000000-acf.json').read_text())
            np.testing.assert_allclose(d['gates'][0]['sigma'],[.1,.1])
            np.testing.assert_allclose(d['gates'][0]['observed'],[[0,1],[2,3]])
            np.testing.assert_allclose(d['gates'][0]['model'],[[2,3],[4,5]])
            self.assertEqual(d['units'],'µK W⁻¹ km⁻¹')
            self.assertEqual(m['raw_url'],'/guisdap/old/raw.f32')
            self.assertEqual(m['native_calibration']['gain'],6e11)
            self.assertEqual((output/'fit-000000.h5').read_bytes(),(fits/'fit-000000.h5').read_bytes())
            with self.assertRaises(AssertionError):
                module.export(historical,fits,output,public)

if __name__=='__main__':unittest.main()
