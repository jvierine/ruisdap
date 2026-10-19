#!/usr/bin/env python3
"""Source-verified native calibration; invoked only by the authenticated server.
Estimation evaluates the native complex ACF response. Applying C creates new
canonical fits and a new immutable viewer analysis, preserving the source.
"""
import fcntl, hashlib, importlib.util, json, os, pathlib, subprocess, sys, tempfile, time
import h5py
import numpy as np

def main(q):
 def phase(message):
  if q.get('job_status_file'):
   path=pathlib.Path(q['job_status_file'])
   with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False)as f:json.dump(dict(status='running',phase=message),f);tmp=pathlib.Path(f.name)
   tmp.chmod(0o600);tmp.replace(path)
 source=pathlib.Path(q['source']);m=json.loads((source/'manifest.json').read_text());review=q['review']
 assert m['analysis_kind']=='native' and m['date']==q['date'] and m['antenna']==q['antenna']
 checked=[];requests=[];rejected=[]
 for point in review['points']:
  if not point.get('included'):continue
  p=dict(point);p.pop('native_unavailable_reason',None)
  assert p['date']==m['date'] and p['antenna']==m['antenna'],'Calibration mixes datasets'
  obs=next((r for r in m['plasma_records'] if r['id']==p['observation_id']),None)
  assert obs is not None,'Unknown plasma observation'
  t=(obs['start']+obs['end'])/2
  assert abs(p['time_unix_s']-t)<.01,'Plasma observation and time disagree'
  assert np.isfinite(p['frequency_standard_error_hz']) and p['frequency_standard_error_hz']>0
  f=min(m['fits'],key=lambda f:abs((f['start_unix_s']+f['end_unix_s'])/2-t))
  i=int(np.argmin(abs(np.array(m['range_m'])-p['range_m'])))
  if abs((f['start_unix_s']+f['end_unix_s'])/2-t)>300 or abs(m['range_m'][i]-p['range_m'])>2500 or f['fit_status'][i]!='converged':
   rejected.append(dict(id=p['id'],reason='No converged native fit within 300 s and 2.5 km'));continue
  path=(source/f['hdf5_url']).resolve()
  assert path.is_relative_to(source.resolve()) and hashlib.sha256(path.read_bytes()).hexdigest()==f['sha256'],'Native fit identity changed'
  p.update(fit_id=f['id'],gate=i,range_m=m['range_m'][i]);checked.append(p)
  requests.append(dict(fit=str(path),range_m=p['range_m'],frequency_hz=p['frequency_hz']))
 assert requests,'Include at least one point with a converged native fit'
 binary=pathlib.Path('target/release/guisdap-millstone').resolve()
 with tempfile.TemporaryDirectory() as temp:
  req=pathlib.Path(temp)/'points.json';req.write_text(json.dumps(requests))
  result=subprocess.run([str(binary),'estimate-calibration',str(req)],capture_output=True,text=True,check=True)
  estimates=json.loads(result.stdout)
 groups={}
 for p,value in zip(checked,estimates):
  p.update(value);groups.setdefault(p['fit_id'],[]).append(p['multiplier'])
 independent=[float(np.median(v))for v in groups.values()];gain=float(np.median(independent));scatter=float(np.median(abs(np.array(independent)-gain))*1.4826)
 answer=dict(magic_constant=gain,magic_constant_units='W m^-3 K^-1 km^-1',revision=review['revision'],independent_fit_count=len(independent),point_count=len(checked),points=checked,rejected_points=rejected,between_fit_scatter=scatter if len(independent)>1 else None,method='Weighted complex ACF projection at Bohm-Gross Ne through the native receiver response; median within fits, then across fits. Ti, Te, velocity and composition are conditional on the source fit.')
 if not q['apply']:return answer
 lock=pathlib.Path(q['private']).parent/'native-compute.lock';lock.parent.mkdir(parents=True,exist_ok=True)
 with lock.open('a') as lease:
  phase('Waiting for the Millstone 128-core compute slot')
  fcntl.flock(lease,fcntl.LOCK_EX)
  phase('Refitting native ACFs on 128 cores')
  digest=hashlib.sha256((json.dumps(review,sort_keys=True)+json.dumps(m['native_calibration'])+hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest()).encode()).hexdigest()
  analysis='native-cal-'+digest[:12];dest=pathlib.Path(q['public'])/'data'/('millstone-'+m['date']+'-'+m['antenna'].lower()+'-'+analysis)
  final=dest
  if dest.exists():
   assert (dest/'calibration.h5').exists() and json.loads((dest/'manifest.json').read_text())['native_calibration'].get('review_sha256')==digest,'Existing calibration publication is incomplete; inspect it before retrying'
  if not dest.exists():
   work=pathlib.Path(q['private'])/digest/str(time.time_ns());work.mkdir(parents=True)
   dest=dest.with_name('.'+dest.name+'-'+str(time.time_ns()))
   inputs=set();options=None
   for f in m['fits']:
    with h5py.File(source/f['hdf5_url']) as h:
     r=h['records'][next(iter(h['records']))]
     sources=json.loads(bytes(h['sources/combined/records_json'][:]))
     for row in sources:
      path=pathlib.Path(row['path']);assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],'Reduced input changed';inputs.add(str(path))
     if options is None:
      model=json.loads(bytes(r['ion_fit/model_json'][:]));settings=json.loads(bytes(r['ion_fit/settings_json'][:]))
      gates=json.loads(bytes(r['ion_fit/results_json'][:]));rg=[g['range_m']/1000 for g in gates]
      options=dict(ions=model['ions'],settings=settings,range_km=[min(rg),max(rg)])
   fit=dict(inputs=sorted(inputs),output=str(work/'fits'),integration_s=m.get('native_integration_s',180),historical_settings=m.get('native_historical_settings',False),windows=[[f['start_unix_s'],f['end_unix_s']]for f in m['fits']],magic_constant_by_antenna={m['antenna']:gain},initial=dict(electron_density_m3=1e11,ion_temperature_k=1000.,electron_temperature_k=2000.,velocity_m_s=0.,ion_collision_s=0.),**options)
   plan=work/'fit.json';plan.write_text(json.dumps(fit))
   env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
   with (work/'fit.log').open('w') as log:
    subprocess.run(['/usr/bin/mpirun','-np','128','--map-by','core','--bind-to','core',str(binary),'fit',str(plan)],env=env,stdout=log,stderr=log,check=True)
   phase('Publishing recalibrated native HDF5 products')
   spec=importlib.util.spec_from_file_location('native_export',pathlib.Path(__file__).with_name('export-millstone-native-view.py'));module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
   # Keep stdout reserved for the structured API response.
   import contextlib
   with contextlib.redirect_stdout(sys.stderr):module.export(source,work/'fits',dest,pathlib.Path(q['public']))
   n=json.loads((dest/'manifest.json').read_text());n['analysis_id']=final.name;n['native_calibration'].update(answer);n['native_calibration']['gain']=gain;n['native_calibration']['review_points']=review['points'];n['native_calibration']['source_analysis']=m['analysis_id'];n['native_calibration']['review_sha256']=digest
   (dest/'manifest.json').write_text(json.dumps(n));(dest/'manifest.json').chmod(0o644)
   with h5py.File(dest/'calibration.h5','w') as h:
    h['review_json_utf8']=np.frombuffer(json.dumps(review).encode(),np.uint8);h['estimate_json_utf8']=np.frombuffer(json.dumps(answer).encode(),np.uint8)
    h['constant']=gain;h['constant'].attrs['units']=answer['magic_constant_units'];h.attrs['source_manifest_sha256']=hashlib.sha256((source/'manifest.json').read_bytes()).hexdigest()
   (dest/'calibration.h5').chmod(0o644)
   dest.rename(final)
  dest=final
  config=pathlib.Path(q['config']);c=json.loads(config.read_text());ds=next(d for d in c['datasets']if d['date']==m['date']and d['antenna']==m['antenna']);choices=ds.setdefault('analyses',[])
  if not any(a['id']==analysis for a in choices):choices.append(dict(id=analysis,label='Rust · plasma-line recalibrated',data_root='data/'+dest.name+'/'))
  with tempfile.NamedTemporaryFile('w',dir=config.parent,delete=False)as f:json.dump(c,f);tmp=pathlib.Path(f.name)
  tmp.chmod(0o644);tmp.replace(config)
  return dict(answer,gain=gain,analysis=analysis,hdf5_url='/guisdap/data/'+dest.name+'/calibration.h5',message='Native ACFs refitted with the reviewed system constant. The source analysis is preserved.')

if __name__=='__main__':
 try:print(json.dumps(main(json.load(sys.stdin)),allow_nan=False))
 except subprocess.CalledProcessError as e:sys.exit(e.stderr or str(e))
 except Exception as e:sys.exit(str(e))
