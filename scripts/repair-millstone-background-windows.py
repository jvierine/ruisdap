#!/usr/bin/env python3
"""Repair only explicitly logged background-gate failures using the verified guard.
Never restart or replace active reduction workers; valid products stay immutable.
"""
import argparse,hashlib,json,pathlib,re,subprocess,os
import h5py
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);a=p.parse_args()
with h5py.File(a.root/'kernel.h5') as h:
 binary=a.root/h.attrs['reduction_binary_name'];sha=h.attrs['reduction_kernel_sha256']
assert hashlib.sha256(binary.read_bytes()).hexdigest()==sha
for date in ['2024-04-08','2024-04-07','2024-04-09']:
 day=a.root/date;log=day/'reduce.log'
 if not log.exists():continue
 windows=sorted(set(map(int,re.findall(r'archive window (\d+): cancelling pulse sample \d+: detection overlaps background window',log.read_text()))))
 for begin in windows:
  if (day/'audit'/f'{begin}.h5').exists():continue
  plan=json.loads((day/'archive.json').read_text());plan.update(start_unix_s=begin,end_unix_s=begin+10,live_directory=None)
  folder=a.root/'repairs'/date;folder.mkdir(parents=True,exist_ok=True);path=folder/f'{begin}.json'
  if path.exists():assert json.loads(path.read_text())==plan
  else:path.write_text(json.dumps(plan))
  with (folder/f'{begin}.log').open('a') as output:
   subprocess.run(['nice','-n','5',str(binary),'process-archive',str(path)],stdout=output,stderr=output,check=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
  assert (day/'audit'/f'{begin}.h5').exists()
  print(date,begin,'REPAIRED',flush=True)
