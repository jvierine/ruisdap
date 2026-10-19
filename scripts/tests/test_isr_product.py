import importlib.util
from pathlib import Path
from types import SimpleNamespace
import tempfile,unittest
import h5py,numpy as np
spec=importlib.util.spec_from_file_location('product',Path(__file__).parents[1]/'isr_product.py');product=importlib.util.module_from_spec(spec);spec.loader.exec_module(product)
class ProductTest(unittest.TestCase):
 def test_millstone_frequency_units_order_and_signed_values(self):
  with tempfile.TemporaryDirectory()as directory:
   src=Path(directory)/'source.h5';dest=Path(directory)/'product.h5'
   with h5py.File(src,'w')as h:
    g=h.create_group('1712599325')
    for k,v in dict(spec=np.array([[-1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]),freqs=np.array([0.,12.5,-12.5]),ranges=np.array([100.,200.,300.]),antenna=np.bytes_('ZENITH'),channel=np.bytes_('zenith-h'),sr=25000000,t1=1712599385,pl=480).items():g[k]=v
   args=SimpleNamespace(radar_frequency_hz=None)
   product.atomic_hdf5(dest,lambda h:product.millstone(h,src,args))
   with h5py.File(dest)as h:
    self.assertEqual(product.validate(h),1);r=next(iter(h['records'].values()))
    np.testing.assert_array_equal(r['plasma_line/frequency_hz'][()],[-12.5e6,0.,12.5e6])
    np.testing.assert_array_equal(r['plasma_line/spectrum'][()],np.array([[7.,-1.,4.],[8.,2.,5.],[9.,3.,6.]]))
    self.assertEqual(r['plasma_line/spectrum'].dtype,np.dtype('float64'))
    original=next(iter(h['sources'].values()));np.testing.assert_array_equal(original['1712599325/spec'][()],np.array([[-1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]))
   with h5py.File(dest,'r+')as h:
    r=next(iter(h['records'].values()));r['plasma_line/frequency_hz'][1]=-12.5e6
    with self.assertRaises(ValueError):product.validate(h)
 def test_fit_ratio_uncertainty_is_not_temperature_error(self):
  with tempfile.TemporaryDirectory()as directory:
   src=Path(directory)/'fit.h5';dest=Path(directory)/'out.h5'
   with h5py.File(src,'w')as h:
    for key,value in dict(ne=[1.,2.],Ti=[1000.,1000.],Te=[1500.,1500.],vi=[0.,0.],rgs=[200.,300.],t0=1.,t1=201.,dne=[.1,.2],dTi=[10.,10.],dvi=[5.,5.]).items():h[key]=value
    h['dTe/Ti']=[.1,.2]
   args=SimpleNamespace(antenna='ZENITH',mode='lpi-30',channel=None,radar_frequency_hz=None)
   product.atomic_hdf5(dest,lambda h:product.millstone(h,src,args))
   with h5py.File(dest)as h:
    r=next(iter(h['records'].values()));self.assertTrue(np.isnan(r['fit/standard_error'][:,2]).all());np.testing.assert_array_equal(r['fit/temperature_ratio_standard_error'][()],[.1,.2])
    self.assertNotIn('altitude_m',r)
 def test_atomic_import_does_not_replace_existing_or_publish_partial(self):
  with tempfile.TemporaryDirectory()as directory:
   path=Path(directory)/'v1.h5';path.write_bytes(b'original')
   with self.assertRaises(ValueError):product.atomic_hdf5(path,lambda h:None)
   self.assertEqual(path.read_bytes(),b'original')
   path2=Path(directory)/'bad.h5'
   with self.assertRaises(ValueError):product.atomic_hdf5(path2,lambda h:None)
   self.assertFalse(path2.exists())
if __name__=='__main__':unittest.main()
