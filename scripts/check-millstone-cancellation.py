#!/usr/bin/env python3
"""Matched original-DigitalRF cancellation controls; HDF5 only, no RF staging.
Recompute only selected source integrations without subtraction, verify exact
pulse identity, then fit both branches using the identical frozen Rust kernel.
"""
import argparse,concurrent.futures,hashlib,json,os,pathlib,shutil,subprocess
import h5py,numpy as np

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def data_json(h,key):return json.loads(bytes(h[key][:]))
def main():
 p=argparse.ArgumentParser();p.add_argument('--campaign',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--binary',type=pathlib.Path,required=True);p.add_argument('--anchor',type=int,action='append',required=True);p.add_argument('--workers',type=int,default=16);a=p.parse_args()
 assert 1<=a.workers<=128
 # Counterfactual controls share the campaign's 128 physical cores. Child
 # workers inherit this affinity; never acquire an additional core pool.
 rank_cpus=set()
 for proc in pathlib.Path('/proc').iterdir():
  if not proc.name.isdigit():continue
  try:
   argv=(proc/'cmdline').read_bytes().split(b'\0')
   if len(argv)>1 and argv[1]==b'process-archive' and pathlib.Path(os.fsdecode(argv[0])).parent==a.campaign:
    rank_cpus.update(os.sched_getaffinity(int(proc.name)))
  except (OSError,UnicodeError):continue
 if not rank_cpus:
  physical={}
  for row in subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True).splitlines():
   if row.startswith('#'):continue
   cpu,core,socket=map(int,row.split(','))
   if cpu in os.sched_getaffinity(0):physical.setdefault((socket,core),cpu)
  rank_cpus=set(list(physical.values())[:128])
 assert rank_cpus,'No available campaign cores'
 os.sched_setaffinity(0,rank_cpus)
 a.output.mkdir(exist_ok=True,parents=True);binary=a.output/'guisdap-millstone-controls'
 if not binary.exists():shutil.copy2(a.binary,binary);binary.chmod(0o700)
 assert digest(binary)==digest(a.binary),'Control kernel changed'
 day=a.campaign/'2024-04-08';base=json.loads((day/'archive.json').read_text());plans={}
 for anchor in a.anchor:
  m=json.loads((day/'stream-jobs-legacy300-v3'/str(anchor)/'fit.json').read_text());plans[anchor]=m
 inputs=sorted(set(pathlib.Path(p) for m in plans.values() for p in m['inputs']))
 windows=sorted(set(int(p.name.split('-')[0]) for p in inputs));env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'};env.pop('DISPLAY',None)
 def reduce(begin):
  folder=a.output/str(begin);folder.mkdir(exist_ok=True)
  m=dict(base,start_unix_s=begin,end_unix_s=begin+10,output=str(folder/'no-cancel'),live_directory=None,cancel_debris=False)
  path=folder/'no-cancel.json'
  if path.exists():assert json.loads(path.read_text())==m
  else:path.write_text(json.dumps(m))
  with (folder/'no-cancel.log').open('a') as log:subprocess.run(['nice','-n','10',str(binary),'process-archive',str(path)],check=True,env=env,stdout=log,stderr=log)
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:list(pool.map(reduce,windows))
 controls={};ledger=[]
 for source in inputs:
  begin=int(source.name.split('-')[0]);control=a.output/str(begin)/'no-cancel/products'/source.name
  with h5py.File(source) as x,h5py.File(control) as y:
   r=x['records'][next(iter(x['records']))];s=y['records'][next(iter(y['records']))]
   assert not data_json(s,'configuration_json')['cancel_debris'],'Control did not disable subtraction'
   assert np.array_equal(r['pulse_sample'][:],s['pulse_sample'][:]),'Control pulse selection changed'
   assert np.array_equal(r['range_m'][:],s['range_m'][:])
   assert np.array_equal(r['ion_line/lag_s'][:],s['ion_line/lag_s'][:])
  controls[str(source)]=str(control);ledger.append(dict(before=str(control),before_sha256=digest(control),after=str(source),after_sha256=digest(source)))
 def fit(job):
  anchor,kind=job;folder=a.output/'integrations'/str(anchor);folder.mkdir(parents=True,exist_ok=True)
  m=dict(plans[anchor],output=str(folder/(kind+'-fits')),range_km=[200.,550.])
  if kind=='no-cancel':m['inputs']=[controls[p] for p in m['inputs']]
  path=folder/(kind+'-fit.json')
  if path.exists():assert json.loads(path.read_text())==m
  else:path.write_text(json.dumps(m))
  with (folder/(kind+'-fit.log')).open('a') as log:subprocess.run(['nice','-n','10',str(binary),'fit',str(path)],check=True,env=env,stdout=log,stderr=log)
  print('MATCHED INTEGRATION READY',anchor,kind,flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=min(a.workers,4)) as pool:list(pool.map(fit,[(t,k) for t in a.anchor for k in ['with-cancel','no-cancel']]))
 with h5py.File(a.output/'matched-controls.h5','w') as h:
  h.attrs['kernel_sha256']=digest(binary);h.attrs['control']='identical pulses, calibration, lag/range response, integration and fitter; subtraction only differs'
  h['sources_json']=np.frombuffer(json.dumps(ledger).encode(),np.uint8)
  for anchor in a.anchor:
   group=h.create_group(str(anchor))
   for kind in ['with-cancel','no-cancel']:
    with h5py.File(a.output/'integrations'/str(anchor)/(kind+'-fits/fit-000000.h5')) as f:
     r=f['records'][next(iter(f['records']))];g=group.create_group(kind)
     for src,dst in [('range_m','range_m'),('fit/parameters','parameters'),('fit/standard_error','standard_error'),('ion_fit/results_json','fit_results_json')]:g[dst]=r[src][:]
 print('MATCHED CONTROL AUDIT',a.output/'matched-controls.h5',flush=True)
if __name__=='__main__':main()
