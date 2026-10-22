#!/usr/bin/env python3
"""Publish actual per-day context using a frozen, explicitly dated calibration.
This transfers the reviewed same-instrument campaign constant; it does not claim
an independently reviewed calibration on days without plasma-line points.
"""
import argparse,importlib.util,json,pathlib,shutil,subprocess,sys,tempfile
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--source',type=pathlib.Path,required=True);p.add_argument('--context',type=pathlib.Path,required=True);p.add_argument('--date',required=True);p.add_argument('--reference',type=pathlib.Path,required=True);a=p.parse_args()
public=pathlib.Path('/srv/guisdap');state=pathlib.Path('/home/j/data/guisdap/catalog-state');out=public/'data'/('millstone-'+a.date+'-zenith-campaign-context-v1')
reference=json.loads((a.reference/'manifest.json').read_text());gain=reference['calibration']['gain'];source=a.context/'zenith';m=json.loads((source/'manifest.json').read_text());assert m['date']==a.date and m['antenna']=='ZENITH'
index=a.context/'source-index.json'
subprocess.run([sys.executable,'scripts/index-millstone-viewer.py',str(a.source),str(a.context),str(index),'--date',a.date],check=True)
subprocess.run([sys.executable,'scripts/export-millstone-raw-calibration.py',str(index)],check=True)
if not out.exists():
 destination=out;out=pathlib.Path(tempfile.mkdtemp(prefix='.millstone-context-',dir=out.parent));out.rmdir()
 shutil.copytree(source,out);out.chmod(0o755)
 with h5py.File(source/'raw-calibration-v2.h5') as h:alpha=h['receiver_gain'][:];az=h['azimuth_deg'][:];el=h['elevation_deg'][:]
 raw=np.fromfile(source/m['raw_url'],'<f4').reshape(len(m['raw_range_m']),len(m['raw_times']));raw=2*gain*raw/alpha[None,:]
 m['raw_units']='m^-3';m['raw_definition']='Raw Ne from measured 10 us ACF, recorded receiver gain and Ptx; Te/Ti=1; short-lag and small-Debye approximation'
 k=4*np.pi*m['radar_frequency_hz']/299792458
 parameters=[]
 for i,f in enumerate(m['fits']):
  values=np.array(f['parameters'],dtype=float);proxy=values[:,0].copy();ne=proxy*gain;te=values[:,2];ti=values[:,1];b=k*k*8.8541878188e-12*1.380649e-23*te/(1.602176634e-19**2);ratio=1+te/ti
  lo=np.zeros_like(ne);hi=np.maximum(ne,1)*100
  for _ in range(80):
   mid=(lo+hi)/2;model=mid*ratio/((1+b/mid)*(ratio+b/mid));below=model<ne;lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
  values[:,0]=np.where((ne>0)&(te>0)&(ti>0),(lo+hi)/2,np.nan);f['uncalibrated_density_proxy']=proxy.tolist();f['parameters']=values.tolist();f['azimuth_deg']=float(az[i]);f['elevation_deg']=float(el[i]);parameters.append(values)
  if 'standard_error'in f:
   errors=np.array(f['standard_error'],dtype=float);f['uncalibrated_standard_error']=errors.tolist();density=values[:,0];power=density*ratio/((1+b/density)*(ratio+b/density));derivative=power*(3/density-1/(density+b)-ratio/(ratio*density+b));errors[:,0]=gain*errors[:,0]/derivative;f['standard_error']=errors.tolist()
 m['parameter_names'][0]='Ne / m⁻³';m['density_semantics']='plasma_line_calibrated_density';m['calibration']=dict(gain=gain,revision=reference['calibration']['revision'],magic_constant_units='W m^-3 K^-1 km^-2',points=[],reference_date=reference['date'],method='Frozen same-antenna April 8 plasma-line calibration transferred to this campaign day; no independent per-day point review claimed')
 np.stack(parameters).transpose(2,1,0).astype('<f4').tofile(out/m['parameters_url']);raw.astype('<f4').tofile(out/m['raw_url'])
 with h5py.File(out/'campaign-calibration.h5','w')as h:
  h['gain']=gain;h['gain'].attrs['units']='W m^-3 K^-1 km^-2';h.attrs['reference_date']=reference['date'];h['raw_ne_m3']=raw;h['fit_parameters']=np.array(parameters);h['receiver_gain']=alpha
 spec=importlib.util.spec_from_file_location('cleaner','scripts/export-millstone-native-view.py');cleaner=importlib.util.module_from_spec(spec);spec.loader.exec_module(cleaner);cleaner.write_json(out/'manifest.json',m)
 for f in out.iterdir():f.chmod(0o644)
 out.rename(destination);out=destination
# Source verification uses the uncalibrated immutable imports, never a second gain.
merged=json.loads((state/'millstone-sources.json').read_text());new=json.loads(index.read_text());merged['plasma'].update(new['plasma']);merged['datasets'].update(new['datasets'])
def atomic(path,value,mode):
 with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False)as f:json.dump(value,f);tmp=pathlib.Path(f.name)
 tmp.chmod(mode);tmp.replace(path)
atomic(state/'millstone-sources.json',merged,0o600)
config=public/'catalog/millstone-config.json';c=json.loads(config.read_text());existing=next((d for d in c['datasets']if d['date']==a.date and d['antenna']=='ZENITH'),None)
if existing is None:c['datasets'].append(dict(date=a.date,antenna='ZENITH',data_root='data/'+out.name+'/',analyses=[]))
atomic(config,c,0o644);print('Published actual day context',a.date,flush=True)
