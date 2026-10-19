#!/usr/bin/env python
"""Read existing products; export a bounded day of linked scientific quicklooks.
The HDF5 output records units, source identities and all displayed signed values.
"""
import argparse,datetime as dt,hashlib,json,os
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);p.add_argument('--date',default='2024-04-08');a=p.parse_args()
start=dt.datetime.fromisoformat(a.date).replace(tzinfo=dt.timezone.utc).timestamp();end=start+86400
if a.output.exists():raise ValueError('choose a new export version')
a.output.mkdir(parents=True)
def finite(x):
 if isinstance(x,np.ndarray):return finite(x.tolist())
 if isinstance(x,list):return [finite(v)for v in x]
 if isinstance(x,dict):return {k:finite(v)for k,v in x.items()}
 if isinstance(x,(float,np.floating))and not np.isfinite(x):return None
 if isinstance(x,np.generic):return x.item()
 return x
def utc(x):return dt.datetime.fromtimestamp(float(x),dt.timezone.utc).strftime('%H:%M:%S')
def f32(path,data):np.asarray(data,dtype='<f4').tofile(path)
def source(path):return dict(path=str(path),size=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns)
for antenna,channel in [('ZENITH','zenith-l'),('MISA','misa-l')]:
 out=a.output/antenna.lower();out.mkdir();fits=[];raw=[];sources=[]
 folder=a.root/'lpi_30'/channel
 for path in sorted(folder.glob('pp-*.h5')):
  t=float(path.stem.split('-')[1])
  if not start<=t<end:continue
  with h5py.File(path)as h:
   if not all(k in h for k in ['ne','Ti','Te','vi','rgs','t0','t1']):continue
   r=h['rgs'][()]*1000.;parameters=np.stack([h[k][()]for k in ['ne','Ti','Te','vi']],axis=1)
   errors=np.full(parameters.shape,np.nan)
   for i,k in enumerate(['dne','dTi','dTe','dvi']):
    if k in h and isinstance(h[k],h5py.Dataset):errors[:,i]=h[k][()]
   digest=hashlib.sha256(path.read_bytes()).hexdigest();rid=hashlib.sha256((digest+':/').encode()).hexdigest()[:24]
   fit=dict(id=rid,start_unix_s=float(h['t0'][()]),end_unix_s=float(h['t1'][()]),range_m=r,parameters=parameters,standard_error=errors,transmitter_w=float(h['P_tx'][()]),system_temperature_k=float(h['T_sys'][()]),elevation_deg=float(h['el'][()]),space_object_count=h['space_object_count'][()]if 'space_object_count'in h else [])
   fits.append(fit);sources.append(source(path)|dict(sha256=digest))
 print(antenna,len(fits),'fit integrations',flush=True)
 if not fits:continue
 ranges=fits[0]['range_m']
 if any(not np.array_equal(f['range_m'],ranges)for f in fits):raise ValueError('fit range grids differ; export modes separately')
 for path in sorted(folder.glob('lpi-*.h5')):
  t=float(path.stem.split('-')[1])
  if not start<=t<end:continue
  with h5py.File(path)as h:
   if not all(k in h for k in ['acfs_e','rgs_km','lags','i0','P_tx']):continue
   r=h['rgs_km'][()];lags=h['lags'][()];lag=int(np.argmin(abs(lags-1e-5)))
   if abs(lags[lag]-1e-5)>1e-9:continue
   power=float(h['P_tx'][()]);values=h['acfs_e'][:,lag].real*r*r/power if power>0 else np.full(len(r),np.nan)
   raw.append(dict(time=float(h['i0'][()]),values=values,range_m=r*1000.));sources.append(source(path))
 print(antenna,len(raw),'raw 10 us profiles',flush=True)
 pl=[];plranges=None;freqs=None
 root=a.root/'metadata'/'integrated_plasma_line_metadata_hires'
 for directory in sorted(root.glob(a.date+'T*')):
  for path in sorted(directory.glob('*.h5')):
   try:
    with h5py.File(path)as h:
     for key,g in h.items():
      if not isinstance(g,h5py.Group)or not all(k in g for k in ['spec','ranges','freqs','antenna','t1']):continue
      name=g['antenna'][()];name=name.decode()if isinstance(name,bytes)else str(name)
      if name!=antenna or not start<=float(key)<end:continue
      rg=g['ranges'][()];use=np.flatnonzero((rg>=185)&(rg<560))
      if len(use)<50:continue
      fq=g['freqs'][()];sr=float(g['sr'][()])
      if np.isclose(np.ptp(fq)*1e6,sr,rtol=.002):fq=fq*1e6
      elif not np.isclose(np.ptp(fq),sr,rtol=.002):raise ValueError('frequency units disagree with receiver sample rate')
      order=np.argsort(fq);fq=fq[order]
      spec=g['spec'][:,int(use[0]):int(use[-1])+1].astype(float)
      spec-=np.nanmedian(spec,axis=0,keepdims=True)
      collapsed=np.nanmax(spec,axis=1)[order];collapsed-=np.nanmedian(collapsed)
      # Average 5 adjacent frequency bins: 10.417 kHz, not interpolation.
      step=5;n=len(fq)//step*step;cf=fq[:n].reshape(-1,step).mean(axis=1);cv=collapsed[:n].reshape(-1,step).mean(axis=1)
      cv[abs(cf)<1e6]=np.nan
      if freqs is None:freqs=cf
      if not np.allclose(freqs,cf):continue
      rid=hashlib.sha256((hashlib.sha256(path.read_bytes()).hexdigest()+':'+key).encode()).hexdigest()[:24]
      pl.append(dict(start=float(key),end=float(g['t1'][()]),values=cv,source=str(path),source_key=key,id=rid))
      sources.append(source(path))
   except (OSError,ValueError)as e:print(path.name,str(e),flush=True)
  print(antenna,directory.name,len(pl),'plasma integrations',flush=True)
 times=np.array([(f['start_unix_s']+f['end_unix_s'])/2 for f in fits]);parameters=np.array([f['parameters']for f in fits]);fitpixels=parameters.transpose(2,1,0)
 rawranges=raw[0]['range_m']if raw else ranges
 if any(not np.array_equal(r['range_m'],rawranges)for r in raw):raise ValueError('raw range grids differ; export modes separately')
 rawtime=np.array([r['time']for r in raw]);rawpixels=np.array([r['values']for r in raw]).T
 pltimes=np.array([(r['start']+r['end'])/2 for r in pl]);plpixels=np.array([r['values']for r in pl]).T
 f32(out/'parameters.f32',fitpixels);f32(out/'raw-ne.f32',rawpixels);f32(out/'plasma.f32',plpixels)
 manifest=dict(schema='ruisdap-millstone-quicklook-v1',station='Millstone Hill',antenna=antenna,date=a.date,radar_frequency_hz=440.2e6,fit_times=times,range_m=ranges,parameters_url='parameters.f32',parameter_names=['Uncalibrated Ne proxy','Ti / K','Te / K','Velocity / m s⁻¹'],fits=fits,raw_times=rawtime,raw_range_m=rawranges,raw_url='raw-ne.f32',raw_definition='Re ACF(10 microseconds) * range_km^2 / transmitter_W; signed uncalibrated proxy; one source integration per time',plasma_times=pltimes,plasma_frequency_hz=freqs if freqs is not None else [],plasma_url='plasma.f32',plasma_records=[{k:v for k,v in r.items()if k not in ['values','source']}for r in pl],collapse_range_km=[185,560],collapse_definition='Per-range frequency median removed; maximum over 185–560 km; then frequency median removed; 5 frequency bins averaged; |offset| < 1 MHz masked in display export',density_semantics='uncalibrated_ne_proxy',frequency_convention='signed receiver offset; no sign reversal',frequency_sign_verified=False,velocity_convention='source convention; sign not independently verified')
 (out/'manifest.json').write_text(json.dumps(finite(manifest),allow_nan=False))
 with h5py.File(out/'quicklook.h5','w')as h:
  h.attrs['schema']=manifest['schema'];h.attrs['metadata_json']=json.dumps(finite({k:v for k,v in manifest.items()if k not in ['fits','plasma_records']}),allow_nan=False);h.attrs['sources_json']=json.dumps(sources)
  for k,v in [('fit_time_unix_s',times),('range_m',ranges),('parameters',parameters),('raw_time_unix_s',rawtime),('raw_range_m',rawranges),('raw_ne_proxy',rawpixels),('plasma_time_unix_s',pltimes),('plasma_frequency_hz',freqs if freqs is not None else []),('collapsed_plasma_spectrum',plpixels)]:h.create_dataset(k,data=v,compression='gzip'if np.asarray(v).size else None)
  for k,unit in [('fit_time_unix_s','UTC Unix s'),('raw_time_unix_s','UTC Unix s'),('plasma_time_unix_s','UTC Unix s'),('range_m','m'),('raw_range_m','m'),('parameters','source_proxy,K,K,m s^-1'),('raw_ne_proxy','source_proxy'),('plasma_frequency_hz','Hz'),('collapsed_plasma_spectrum','source_spectral_power')]:h[k].attrs['units']=unit
 print('Export complete:',antenna,len(fits),len(raw),len(pl),flush=True)
