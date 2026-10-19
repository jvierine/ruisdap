#!/usr/bin/env python3
"""Progress only for the explicitly authorized direct campaign. Never restarts
jobs or modifies RF or reduced scientific products; exits when the run stops.
"""
import argparse,json,pathlib,subprocess,tempfile,time
import h5py
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--unit',required=True);a=p.parse_args();seen=set();counts={}
while True:
 state=subprocess.run(['systemctl','--user','is-active',a.unit],capture_output=True,text=True).stdout.strip();active=state in ['active','activating'];path=a.root/'live/rank-0.json'
 if path.exists():
  try:q=json.loads(path.read_text())
  except (OSError,ValueError):time.sleep(10);continue
  phase=q.get('phase','');date=phase[:10]
  if active and 'reduction' in phase:
   audit=list((a.root/date/'audit').glob('*.h5'));counter=counts.setdefault(date,[0,0])
   for file in audit:
    if file in seen:continue
    with h5py.File(file)as h:counter[0]+=len(h['science_pulse_sample']);counter[1]+=len(h['candidate_science_pulse_sample'])
    seen.add(file)
   q.update(completed=len(audit),total=8640,science_pulse_count=counter[0],candidate_science_pulse_count=counter[1])
  elif not active and not phase.startswith('all available ZENITH eclipse days complete'):
   q['phase']='Direct campaign stopped; inspect processing log before continuing'
  q['updated_unix_s']=time.time()
  with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False)as h:json.dump(q,h);tmp=pathlib.Path(h.name)
  if json.loads(path.read_text()).get('phase','')==phase:tmp.replace(path)
  else:tmp.unlink()
 if not active:break
 time.sleep(10)
