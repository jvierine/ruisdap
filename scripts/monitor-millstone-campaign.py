#!/usr/bin/env python3
"""Progress heartbeat for one explicitly launched campaign, not a scheduler.
Exits when the owned run stops. Never modifies a numerical product. ACF
previews remain subject to the GUI lease. Controller phase changes win races.
"""
import argparse
import json
import pathlib
import subprocess
import tempfile
import time
p=argparse.ArgumentParser()
p.add_argument('--root',type=pathlib.Path,required=True)
p.add_argument('--unit',required=True)
a=p.parse_args()
offsets={}
counts={}
while True:
 state=subprocess.run(['systemctl','--user','is-active',a.unit],capture_output=True,text=True).stdout.strip()
 active=state in ['active','activating']
 path=a.root/'live/rank-0.json'
 if path.exists():
  try:q=json.loads(path.read_text())
  except (OSError,ValueError):time.sleep(5);continue
  old_phase=q.get('phase','')
  if not active:
   if 'complete' not in old_phase:q['phase']='Campaign stopped; inspect processing log before continuing'
  elif 'preparing' in old_phase:
   date=old_phase[:10]
   for log in (a.root/date).glob('prepare-*.log'):
    with log.open()as f:
     f.seek(offsets.get(log,0))
     while True:
      begin=f.tell();line=f.readline()
      if not line or not line.endswith('\n'):f.seek(begin);break
      try:row=json.loads(line)
      except ValueError:continue
      # Final per-hour dictionaries summarize the scalar pulse rows already counted.
      if isinstance(row.get('science_pulses'),int):counts[date]=counts.get(date,0)+row['science_pulses']
     offsets[log]=f.tell()
   q['phase']=date+' preparing inputs · '+format(counts.get(date,0),',')+' science pulses packed'
   q['packed_science_pulses']=counts.get(date,0)
  q['updated_unix_s']=time.time()
  with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False)as f:
   json.dump(q,f);tmp=pathlib.Path(f.name)
  if json.loads(path.read_text()).get('phase','')==old_phase:tmp.replace(path)
  else:tmp.unlink()
 if not active:break
 time.sleep(5)
