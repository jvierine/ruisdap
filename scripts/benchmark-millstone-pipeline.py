#!/usr/bin/env python3
"""Bounded benchmark on the frozen 1122-pulse real-data pilot. Numerical outputs HDF5.
No claim of 128-rank scaling: this pilot contains only two independent packs.
"""
import hashlib,json,os,pathlib,subprocess,time,h5py,numpy as np
base=pathlib.Path('/home/j/data/guisdap/20261019-native-campaign-pilot-v2');out=pathlib.Path('/home/j/data/guisdap/20261019-lpi-benchmark-v1')
binary=pathlib.Path('/home/j/data/guisdap/20261019-native-eclipse-all-pulses-v1/guisdap-millstone');report=out/'pipeline.h5';assert not report.exists()
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
with h5py.File(report,'w')as h:
 h.attrs['kernel_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest();h.attrs['workload']='same two actual ten-second packs, no replicated observations'
 times=[]
 for ranks in [1,2]:
  for rep in range(2):
   work=out/f'mpi-{ranks}-repeat-{rep}';m=json.loads((base/'manifest.json').read_text());m.update(output=str(work),live_directory=None);path=out/f'mpi-{ranks}-{rep}.json';path.write_text(json.dumps(m))
   command=['/usr/bin/mpirun','-np',str(ranks),'--cpu-set','122,123','--bind-to','core',str(binary),'process',str(path)]
   start=time.perf_counter()
   with(out/f'mpi-{ranks}-{rep}.log').open('w')as log:subprocess.run(command,env=env,stdout=log,stderr=log,check=True)
   elapsed=time.perf_counter()-start;times.append([ranks,rep,elapsed])
   if 'process_seconds'in h:del h['process_seconds']
   h['process_seconds']=times;h.flush();print(json.dumps(dict(ranks=ranks,repeat=rep,process_seconds=elapsed)),flush=True)
   # Numerical products must agree between scheduling choices and stage harness.
   for i,p in enumerate(sorted(work.glob('*.h5'))):
    with h5py.File(p)as v,h5py.File(out/f'bench-pack-{i}-repeat-0.h5')as ref:
     rv=next(iter(v['records'].values()))['ion_line'];rr=next(iter(ref['records'].values()))['ion_line']
     for name in ['acf','ground_clutter_acf','variance','retained_fraction','receiver_gain']:
      np.testing.assert_array_equal(rv[name][:],rr[name][:])
 # One conditional fitting workload, separate from decoding.
 f=json.loads((base/'fit.json').read_text());f['output']=str(out/'fit');path=out/'fit.json';path.write_text(json.dumps(f));start=time.perf_counter()
 with(out/'fit.log').open('w')as log:subprocess.run(['taskset','-c','123',str(binary),'fit',str(path)],env=env,stdout=log,stderr=log,check=True)
 elapsed=time.perf_counter()-start;h['fit_seconds']=elapsed;h.flush();print(json.dumps(dict(fit_seconds=elapsed)),flush=True)
