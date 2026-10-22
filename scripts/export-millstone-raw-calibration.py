#!/usr/bin/env python
"""Export the recorded receiver gains for an existing immutable day export."""
import argparse,json,os
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('index',type=Path);a=p.parse_args()
for name,d in json.loads(a.index.read_text())['datasets'].items():
 m=json.loads(Path(d['manifest']).read_text());dest=Path(d['manifest']).with_name('raw-calibration-v2.h5')
 if dest.exists():print(name,'already exported',flush=True);continue
 lookup={round(v['time'],3):v['path'] for v in d['raw']};gains=[];sources=[]
 for i,t in enumerate(m['raw_times']):
  path=Path(lookup[round(t,3)])
  with h5py.File(path) as h:
   assert abs(float(h['i0'][()])-t)<.01
   gain=float(h['alpha'][()]);gain=gain if np.isfinite(gain) and gain>0 else np.nan
   gains.append(gain);sources.append(str(path))
  if i%500==0:print(name,i,'/',len(m['raw_times']),flush=True)
 az=[];el=[]
 for f in m['fits']:
  path=Path(d['raw'][0]['path']).parent/('pp-%d.h5'%f['start_unix_s'])
  with h5py.File(path) as h:
   assert abs(float(h['t0'][()])-f['start_unix_s'])<.01
   az.append(float(h['az'][()]));el.append(float(h['el'][()]))
  sources.append(str(path))
 temp=dest.with_suffix('.tmp.h5')
 with h5py.File(temp,'x') as h:
  h.attrs['schema']='ruisdap-millstone-raw-calibration-v1'
  h.attrs['source_paths_json']=json.dumps(sources)
  h.create_dataset('time_unix_s',data=m['raw_times']).attrs['units']='UTC Unix s'
  h.create_dataset('receiver_gain',data=gains).attrs['units']='source power K^-1'
  h.create_dataset('fit_start_unix_s',data=[f['start_unix_s'] for f in m['fits']])
  h.create_dataset('azimuth_deg',data=az).attrs['units']='degrees clockwise from north'
  h.create_dataset('elevation_deg',data=el).attrs['units']='degrees above horizon'
  h.attrs['missing_receiver_gain_count']=int(np.count_nonzero(~np.isfinite(gains)))
  h.attrs['raw_definition']='Re ACF(10 us) * range_km^2 / recorded transmitter_W; receiver gain NOT yet removed'
 os.rename(temp,dest);print(name,len(gains),'receiver gains exported',flush=True)
