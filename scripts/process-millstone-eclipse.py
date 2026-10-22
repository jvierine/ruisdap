#!/usr/bin/env python3
"""Run all available eclipse days with 128 physical MPI cores and explicit audits.
No automatic restart. The immutable numerical executable is frozen before input
preparation. Historical products remain available throughout the campaign run.
"""
import argparse,concurrent.futures,datetime as dt,fcntl,hashlib,json,os,pathlib,shutil,subprocess,sys,time
import h5py,numpy as np
from isr_product import validate,text
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--context-apr7',type=pathlib.Path,required=True);p.add_argument('--context-apr9',type=pathlib.Path,required=True);a=p.parse_args()
assert not a.output.exists(),'Choose a fresh, immutable campaign version'
a.output.mkdir(parents=True);(a.output/'live').mkdir();binary=a.output/'guisdap-millstone';shutil.copy2('target/release/guisdap-millstone',binary);binary.chmod(0o700);kernel_sha=hashlib.sha256(binary.read_bytes()).hexdigest()
public=pathlib.Path('/srv/guisdap');reference=public/'data/millstone-2024-04-08-zenith-cal-d9ac4d6a1088';ref=json.loads((reference/'manifest.json').read_text());gain=ref['calibration']['gain']*480*(299792.458*30/(2*1e6))
with h5py.File(a.output/'kernel.h5','w')as h:h.attrs['binary_sha256']=kernel_sha;h['calibration_reference_utf8']=np.frombuffer(json.dumps(ref['calibration']).encode(),np.uint8)
available=sorted({p.name[:10]for p in (a.root/'metadata/id_metadata').iterdir()if p.is_dir()});assert set(available)=={'2024-04-07','2024-04-08','2024-04-09'},'Unexpected campaign dates; inspect archive before selecting scope'
order=['2024-04-08','2024-04-07','2024-04-09'];env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'};env.pop('DISPLAY',None)
def status(phase,completed,total):
 q=dict(host='avaruus',rank=0,ranks=1,phase=phase,completed=completed,total=total,updated_unix_s=time.time(),elapsed_s=0);path=a.output/'live/rank-0.json';tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(q));tmp.replace(path)
mpi=['/usr/bin/mpirun','-np','128','--map-by','core','--bind-to','core']
for date in order:
 day=a.output/date;day.mkdir();(day/'live').mkdir();start=int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp())
 def prepare(hour):
  out=day/f'hour-{hour:02d}'
  with (day/f'prepare-{hour:02d}.log').open('w')as log:subprocess.run([sys.executable,'scripts/prepare-millstone-campaign.py','--root',str(a.root),'--start',str(start+3600*hour),'--scan-seconds','3600','--output',str(out)],env=env,stdout=log,stderr=log,check=True)
  return hour,json.loads((out/'manifest.json').read_text())
 inputs=[];status(date+' preparing all pulses',0,24)
 with concurrent.futures.ThreadPoolExecutor(max_workers=12)as pool:
  for done,f in enumerate(concurrent.futures.as_completed([pool.submit(prepare,h)for h in range(24)]),1):
   hour,m=f.result();inputs.extend(m['inputs']);print(date,'prepared',hour,len(m['inputs']),'integrations',m.get('pulse_counts'),flush=True);status(date+' preparing all pulses',done,24)
 assert inputs,'No supported ZENITH science pulses on campaign day'
 weights={};keys=[];hashes=[];counts={'coded':0,'long':0}
 for name in inputs:
  with h5py.File(name)as h:
   science=h['pulse_sample'][:][h['science_pulse'][:].astype(bool)];keys.append(science);kind='long'if 'long' in pathlib.Path(name).stem else'coded';counts[kind]+=len(science);weights[name]=len(h['pulse_sample'])+120*len(h['detections'])
  hashes.append(hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest().encode())
 keys=np.concatenate(keys);assert len(np.unique(keys))==len(keys),'Duplicate science pulse across integrations'
 with h5py.File(day/'input-audit.h5','w')as h:
  h['science_pulse_sample']=keys;h['source_sha256']=np.array(hashes,dtype='S64');h.attrs['binary_sha256']=kernel_sha;h.attrs['counts_json']=json.dumps(counts);h.attrs['frequency_hz']=440.2e6;h.attrs['lag_count']=45;h.attrs['range_gate_count']=240
 inputs.sort(key=lambda n:(-weights[n],n));plan=dict(inputs=inputs,output=str(day/'products'),live_directory=str(day/'live'));path=day/'process.json';path.write_text(json.dumps(plan))
 print(date,'all prepared',len(inputs),'integrations',counts,'first',float(np.min(keys))/1e6,'last',float(np.max(keys))/1e6,flush=True)
 lock=pathlib.Path('/home/j/data/guisdap/2022-11-02/calibration/native-compute.lock')
 lock.parent.mkdir(parents=True,exist_ok=True)
 with lock.open('a')as lease:
  fcntl.flock(lease,fcntl.LOCK_EX);assert hashlib.sha256(binary.read_bytes()).hexdigest()==kernel_sha
  status(date+' reducing all pulses on 128 cores',0,len(inputs));subprocess.run([*mpi,str(binary),'process',str(path)],env=env,check=True)
  products=sorted((day/'products').glob('*.h5'));assert len(products)==len(inputs),'Missing native reduction products'
  for product in products:
   with h5py.File(product)as h:validate(h)
  if date=='2024-04-08':context=reference
  else:
   source_context=a.context_apr7 if date=='2024-04-07'else a.context_apr9
   assert (source_context/'zenith/manifest.json').exists(),'Per-day source context not ready'
   subprocess.run([sys.executable,'scripts/publish-millstone-campaign-context.py','--source',str(a.root),'--context',str(source_context),'--date',date,'--reference',str(reference)],env=env,check=True)
   context=public/'data'/('millstone-'+date+'-zenith-campaign-context-v1')
  old=json.loads((context/'manifest.json').read_text());assert old['date']==date
  windows=[[f['start_unix_s'],f['end_unix_s']+10]for f in old['fits']]
  fit=dict(inputs=list(map(str,products)),output=str(day/'fits'),integration_s=200.,windows=windows,magic_constant_by_antenna={'ZENITH':gain},ions=[dict(mass_amu=16.,fraction=1.)],range_km=[90.,1075.],historical_settings=True,initial=dict(electron_density_m3=1e11,ion_temperature_k=1000.,electron_temperature_k=2000.,velocity_m_s=0.,ion_collision_s=0.),settings=dict(tie_electron_ion_temperature=False,fit_collision_rate=False,fixed_collision_s=0.,max_iterations=160,millstone_historical_bounds=True))
  path=day/'fit.json';path.write_text(json.dumps(fit));status(date+' fitting historical windows and range averages',0,1);subprocess.run([*mpi,str(binary),'fit',str(path)],env=env,check=True)
  destination=public/'data'/('millstone-'+date+'-zenith-native-complete-20261019-v1')
  subprocess.run([sys.executable,'scripts/export-millstone-native-view.py','--historical',str(context),'--fits',str(day/'fits'),'--output',str(destination)],env=env,check=True)
  native=json.loads((destination/'manifest.json').read_text());native['analysis_default_window']=[start,start+86400];native['science_pulse_count']=len(keys);native['pulse_counts']=counts;native['kernel_sha256']=kernel_sha;native['analysis_note']='All supported powered monostatic ZENITH coded and 480 us long pulses available in the archive. Ten-second 45-lag inversions on the original 0–1075 km grid; historical integration windows, below-300 km clutter subtraction, range averaging and Holt molecular composition. Above 700 km O+/H+ is profiled on a grid. The physical finite-response, finite-Debye Rust fitter and propagated diagonal errors differ from the old interpolated free-amplitude fitter. Calibration initially transfers the reviewed April 8 constant; per-analysis plasma-line recalibration is available.'
  (destination/'manifest.json').write_text(json.dumps(native));(destination/'manifest.json').chmod(0o644)
  config=public/'catalog/millstone-config.json';c=json.loads(config.read_text());d=next(d for d in c['datasets']if d['date']==date and d['antenna']=='ZENITH');d.setdefault('analyses',[]).append(dict(id='native-complete-20261019',label='Rust · all coded + long pulses · 45 lags',data_root='data/'+destination.name+'/'))
  import tempfile
  with tempfile.NamedTemporaryFile('w',dir=config.parent,delete=False)as h:json.dump(c,h);tmp=pathlib.Path(h.name)
  tmp.chmod(0o644);tmp.replace(config)
 print(date,'PUBLISHED',len(native['fits']),'profiles',counts,flush=True)
status('all available ZENITH eclipse days complete',3,3)
