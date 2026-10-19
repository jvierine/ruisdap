#!/usr/bin/env python3
"""Compare the native read_vector API with upstream DigitalRF on synthetic
compressed multi-directory fixtures and original Millstone RF, without staging
any campaign RF. Numerical results are retained only in an HDF5 test report.
"""
import argparse,json,pathlib,subprocess,tempfile
import digital_rf,h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--binary',type=pathlib.Path,required=True);p.add_argument('--report',type=pathlib.Path,required=True);p.add_argument('--millstone-root',type=pathlib.Path);a=p.parse_args();assert not a.report.exists()
results=[]
with tempfile.TemporaryDirectory(prefix='digital-rf-synthetic-test-')as tmp:
 t=pathlib.Path(tmp);roots=[]
 dtype=np.dtype([('r','<i2'),('i','<i2')])
 for i,start in enumerate([35990,36000]):
  root=t/f'synthetic-{i}';channel=root/'ch';channel.mkdir(parents=True);roots.append(root)
  with digital_rf.DigitalRFWriter(str(channel),dtype,3600,1000,start,10,1,compression_level=1,checksum=True,is_complex=True,num_subchannels=2,is_continuous=False)as w:
   z=np.zeros((10,2),dtype);z['r'][:,0]=np.arange(i*10,(i+1)*10);z['r'][:,1]=z['r'][:,0]+1000;z['i']=-z['r'];w.rf_write(z)
 def check(name,paths,channel,start,length,should_fail=False):
  reference=digital_rf.DigitalRFReader(list(map(str,paths)),rdcc_nbytes=64*1024*1024)
  try:y=reference.read_vector(start,length,channel);ref_error=False
  except IOError:ref_error=True
  out=t/(name+'.h5');q=subprocess.run([str(a.binary),str(out),channel,str(start),str(length),*map(str,paths)],capture_output=True,text=True);native_error=q.returncode!=0
  assert ref_error==native_error==(should_fail),dict(name=name,reference=ref_error,native=native_error,error=q.stderr)
  if not should_fail:
   with h5py.File(out)as h:
    v=h['vector'][:];z=(v[:,0]+1j*v[:,1]).reshape(tuple(h['shape'][:]))
    np.testing.assert_array_equal(z,y);cache=int(h.attrs['cache_bytes'])
   results.append(dict(name=name,samples=length,shape=list(y.shape),reference_dtype=str(y.dtype),cache_bytes=cache,exact=True))
  else:results.append(dict(name=name,matched_gap_error=True))
 check('two_roots_hour_boundary',roots,'ch',35995,10)
 check('one_file',roots,'ch',35991,4)
 check('missing_file',roots,'ch',36005,10,True)
 if a.millstone_root:
  r=a.millstone_root/'rf_data';check('real_pulse',[r],'zenith-l',1712599210000000,10000)
  check('real_file_boundary',[r],'zenith-l',1712599210999950,100)
  check('real_missing_interval',[r],'zenith-l',1712564270000000,10000,True)
  check('real_after_gap',[r],'zenith-l',1712564298000000,10000)
with h5py.File(a.report,'w')as h:h['checks_json']=np.frombuffer(json.dumps(results).encode(),np.uint8);h.attrs['reference']='digital_rf '+digital_rf.__version__
print(json.dumps(results,indent=2))
