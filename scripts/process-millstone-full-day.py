#!/usr/bin/env python3
"""Reduce every eligible 30-us ZENITH integration in a complete UTC day.
One input reader per UTC hour precedes 128 MPI ranks. HDF5 stays authoritative. Publication
follows successful reduction and fitting. No automatic restart is installed.
"""
import argparse, concurrent.futures, datetime as dt, json, os, pathlib, subprocess, sys, time
p=argparse.ArgumentParser()
for name in ['root','output','historical','publish']:
 p.add_argument('--'+name,type=pathlib.Path,required=True)
p.add_argument('--date',required=True)
a=p.parse_args()
assert not (a.output/'products').exists(),'Resume is limited to preparation; do not overwrite reduced products'
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'live').mkdir(exist_ok=True)
start=int(dt.datetime.fromisoformat(a.date).replace(tzinfo=dt.timezone.utc).timestamp())
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'};env.pop('DISPLAY',None)
source=json.loads((a.historical/'manifest.json').read_text())
assert source['date']==a.date and source['antenna']=='ZENITH'
gain=source['calibration']['gain']*480*(299792.458*30/(2*1e6))

def status(phase,done,total):
 value=dict(host='avaruus',rank=0,ranks=1,phase=phase,completed=done,total=total,elapsed_s=0,updated_unix_s=time.time())
 path=a.output/'live/rank-0.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value));tmp.replace(path)

def prepare(hour):
 out=a.output/f'hour-{hour:02d}'
 if (out/'manifest.json').exists():
  rows=json.loads((out/'manifest.json').read_text())['inputs']
  assert all(pathlib.Path(p).is_file() for p in rows)
  return hour,rows
 cmd=[sys.executable,'scripts/prepare-millstone-native-run.py','--root',str(a.root),'--start',str(start+hour*3600),'--scan-seconds','3600','--jobs','0','--antenna','ZENITH','--output',str(out),'--resume']
 with (a.output/f'prepare-{hour:02d}.log').open('a') as log:
  subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 return hour,json.loads((out/'manifest.json').read_text())['inputs']

status('preparing full-day inputs',0,24)
inputs=[];completed=0
with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:
 for f in concurrent.futures.as_completed([pool.submit(prepare,h)for h in range(24)]):
  hour,rows=f.result();inputs.extend(rows);completed+=1
  print(f'Prepared hour {hour:02d}: {len(rows)} integrations; {completed}/24 hours',flush=True)
  status('preparing full-day inputs',completed,24)
inputs.sort()
assert inputs,'No supported coded integrations'
manifest=dict(inputs=inputs,output=str(a.output/'products'),live_directory=str(a.output/'live'))
(a.output/'manifest.json').write_text(json.dumps(manifest))
print(f'Processing all {len(inputs)} eligible ZENITH integrations on 128 cores',flush=True)
mpi=['/usr/bin/mpirun','-np','128','--map-by','core','--bind-to','core']
binary=str(pathlib.Path('target/release/guisdap-millstone').resolve())
subprocess.run([*mpi,binary,'process',str(a.output/'manifest.json')],env=env,check=True)
products=sorted((a.output/'products').glob('*.h5'))
assert len(products)==len(inputs),'Native product count differs from eligible input count'
fit=dict(inputs=list(map(str,products)),output=str(a.output/'fits'),integration_s=180.,magic_constant_by_antenna={'ZENITH':gain},ions=[dict(mass_amu=16.,fraction=1.)],range_km=[200.,540.],initial=dict(electron_density_m3=1e11,ion_temperature_k=1000.,electron_temperature_k=2000.,velocity_m_s=0.,ion_collision_s=0.),settings=dict(tie_electron_ion_temperature=False,fit_collision_rate=False,fixed_collision_s=0.,max_iterations=80))
(a.output/'fit.json').write_text(json.dumps(fit))
status('fitting full-day profiles',0,1)
subprocess.run([*mpi,binary,'fit',str(a.output/'fit.json')],env=env,check=True)
subprocess.run([sys.executable,'scripts/export-millstone-native-view.py','--historical',str(a.historical),'--fits',str(a.output/'fits'),'--output',str(a.publish)],env=env,check=True)
config=pathlib.Path('/srv/guisdap/catalog/millstone-config.json')
c=json.loads(config.read_text());d=next(d for d in c['datasets']if d['antenna']=='ZENITH'and d['date']==a.date)
choices=d.setdefault('analyses',[]);choices[:]=[x for x in choices if x['id']!='native-full-20261019']
choices.append(dict(id='native-full-20261019',label='Rust experimental · full-day 30 µs code',data_root='data/'+a.publish.name+'/'))
tmp=config.with_suffix('.native.tmp');tmp.write_text(json.dumps(c));tmp.chmod(0o644);tmp.replace(config)
status('complete',1,1)
print('Full-day native analysis published',flush=True)
