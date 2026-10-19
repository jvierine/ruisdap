#!/usr/bin/env python
"""Prepare distinct, observed coded-pulse integrations for the native MPI worker.
RF samples and CFAR records remain unchanged. Only monostatic, powered pulses
with real adjacent same-code clutter partners enter an integration.
"""
import argparse,bisect,datetime as dt,hashlib,importlib.util,json
from pathlib import Path
import h5py,numpy as np
spec=importlib.util.spec_from_file_location('pulse_io',Path(__file__).with_name('pack-millstone-pulses.py'));io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--start',type=int,required=True);p.add_argument('--scan-seconds',type=int,default=600);p.add_argument('--jobs',type=int,default=128);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
if a.output.exists():raise ValueError('choose a new run directory')
a.output.mkdir(parents=True);(a.output/'inputs').mkdir()
def hour(t):return dt.datetime.fromtimestamp(t,dt.timezone.utc).strftime('%Y-%m-%dT%H-00-00')
ants=[]
for hdir in sorted({hour(t)for t in range(a.start-3600,a.start+a.scan_seconds+3600,3600)}):
 for f in sorted((a.root/'metadata/antenna_control_metadata'/hdir).glob('*.h5')):
  with h5py.File(f)as h:
   for k,g in h.items():ants.append((int(k),{n:g[n][()]for n in g}))
ants.sort(key=lambda r:r[0]);anttimes=[r[0]for r in ants]
def antenna(k):
 i=bisect.bisect_right(anttimes,k)-1
 if i<0:return None
 r=ants[i][1]
 if r['rx_antenna']!=r['tx_antenna']:return None
 return r['rx_antenna'].decode(),r
inputs=[];catalogues={};rejected=[];ndet=0
for second in range(a.start,a.start+a.scan_seconds):
 if len(inputs)>=a.jobs:break
 name=hour(second);idsfile=a.root/'metadata/id_metadata'/name/f'md@{second}.h5'
 if not idsfile.exists():continue
 with h5py.File(idsfile)as h:ids={int(k):int(g['sweepid'][()])for k,g in h.items()if 'sweepid'in g}
 if not any(1<=v<=32 for v in ids.values()):continue
 powerfile=a.root/'metadata/powermeter'/name/f'power@{second}.h5'
 if not powerfile.exists():continue
 with h5py.File(powerfile)as h:
  powers=[(int(k),{n:float(g[n][()])for n in g})for k,g in h.items()]
 if not powers:continue
 power=min(powers,key=lambda r:abs(r[0]-second))[1]
 # A job contains one antenna and one-second RF source; skip switching pulses.
 for ant,channel in [('ZENITH','zenith-l'),('MISA','misa-l')]:
  pt=power['zenith_power'if ant=='ZENITH'else'misa_power']
  if not np.isfinite(pt)or pt<400e3:continue
  keys=sorted(k for k,v in ids.items()if 1<=v<=32 and second*1000000<=k and k+10000<=(second+1)*1000000 and (antenna(k)or('',))[0]==ant)
  paired={}
  for i,k in enumerate(keys):
   for j in [i-1,i+1]:
    if 0<=j<len(keys)and ids[keys[j]]==ids[k]and abs(keys[j]-k)<20000:paired[k]=keys[j];break
  keys=sorted(paired)
  if len(keys)<16:continue
  rf=a.root/'rf_data'/channel/name/f'rf@{second}.000.h5';tf=a.root/'rf_data/tx-h'/name/f'rf@{second}.000.h5'
  if not rf.exists()or not tf.exists():continue
  base,rx=io.rf(rf);tbase,txref=io.rf(tf)
  if base!=tbase:raise ValueError('RF clocks disagree')
  y=np.array([rx[k-base:k-base+10000]for k in keys]);transmit=np.zeros_like(y);transmit[:,76:624]=y[:,76:624]
  delays=[io.channel_delay(y[i],txref[k-base:k-base+10000],[76,624])for i,k in enumerate(keys)]
  catalogue=a.root/'metadata/satellite_detections'/name/f'satellite@{second//3600*3600}.h5'
  if catalogue not in catalogues:catalogues[catalogue]=h5py.File(catalogue,'r')
  cat=catalogues[catalogue];detections=[];cfar=[]
  for i,k in enumerate(keys):
   if str(k)not in cat:continue
   g=cat[str(k)]
   if g['channel'][()].decode()!=channel:continue
   for j,(raw,dop)in enumerate(zip(g['raw_delay_sample'][()],g['doppler_hz'][()])):
    detections.append([i,float(raw)-76.-delays[i],float(dop)])
    cfar.append(dict(pulse_sample=k,raw_delay_sample=int(raw),doppler_hz=float(dop),snr_db=float(g['snr_db'][()][j]),cfar_ratio=float(g['cfar_ratio'][()][j]),cfar_threshold=float(g['cfar_threshold'][()][j])))
  cfg=dict(station='Millstone Hill',antenna=ant,channel=channel,mode='alternating-code-16',radar_frequency_hz=440.2e6,start_unix_s=keys[0]/1e6,end_unix_s=(keys[-1]+10000)/1e6,
   lpi=dict(code_energy_samples=[480.]*len(keys),sample_rate_hz=1e6,decimation=30,lag_samples=list(range(10,161,10)),range_start=28,range_end=120,measurement_start=34,measurement_end=240,tx_gate=[76,624],echo_gate=[1000,8200],noise_gate=[7700,8200],injection_gate=[8400,8850],filter_half=100,passband_hz=18000.,outlier_threshold=10.,background_column=True),cancel_debris=True,cancellation=dict(refine_delay_samples=4.,refine_doppler_hz=100.),injection_temperature_k=1172.,range_doppler_fft=1024,range_doppler_samples=[],
   instrument_metadata=dict(transmitter_w=pt,bit_length_us=30.,code_bits=16,pulse_energy_definition='16 equal-amplitude 30-us bauds at 1 MHz sample rate; measured template retains the received pulse shape',tx_to_receiver_delay_samples=delays,cfar=cfar,azimuth_deg=0. if ant=='ZENITH'else float(antenna(keys[0])[1]['misa_azimuth']),elevation_deg=90. if ant=='ZENITH'else float(antenna(keys[0])[1]['misa_elevation'])))
  out=a.output/'inputs'/f'{second}-{ant}.h5'
  with h5py.File(out,'w')as h:
   h['schema_utf8']=np.frombuffer(b'ruisdap-pulses-v1',np.uint8);h['configuration_json']=np.frombuffer(json.dumps(cfg).encode(),np.uint8)
   h['pulse_sample']=np.array(keys,np.uint64);h['sweep_id']=np.array([ids[k]for k in keys],np.uint16);h['pair_index']=np.array([keys.index(paired[k])for k in keys],np.uint64)
   for label,z in [('echo',y),('transmit',transmit)]:h.create_dataset(label,data=np.stack([z.real,z.imag],axis=-1),compression='gzip',compression_opts=1)
   h['detections']=np.asarray(detections,float).reshape(-1,3)
   for label,f in [('receive',rf),('transmit',tf),('id_metadata',idsfile),('power',powerfile)]:h.attrs[label+'_sha256']=hashlib.sha256(f.read_bytes()).hexdigest();h.attrs[label+'_path']=str(f)
   h.attrs['catalogue_path']=str(catalogue);h.attrs['scope']='experimental native reduction; 16 measured ion lags, 126-540 km; no fitted plasma parameters in this worker stage'
  inputs.append(str(out));ndet+=len(detections)
  print(json.dumps(dict(job=len(inputs),second=second,antenna=ant,pulses=len(keys),detections=len(detections))),flush=True)
for h in catalogues.values():h.close()
if len(inputs)<a.jobs:raise ValueError(f'only {len(inputs)} valid distinct integrations in requested scan')
(a.output/'manifest.json').write_text(json.dumps(dict(inputs=inputs,output=str(a.output/'products'),live_directory=str(a.output/'live'))))
(a.output/'single.json').write_text(json.dumps(dict(inputs=inputs[:1],output=str(a.output/'control'),live_directory=str(a.output/'control-live'))))
print(json.dumps(dict(jobs=len(inputs),detections=ndet,manifest=str(a.output/'manifest.json'))))
