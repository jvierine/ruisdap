#!/usr/bin/env python3
"""Pack all supported powered ZENITH coded and long pulses in ten-second windows.
Python performs archive I/O only; Rust does DC subtraction, CFAR-seeded complex
voltage cancellation, receiver filtering and the 45-lag inversion. Archive RF is
immutable. Unsupported sweep 800 is the separate plasma-line waveform.
"""
import argparse,bisect,collections,datetime as dt,hashlib,importlib.util,json
from pathlib import Path
import h5py,numpy as np
from millstone_archive import recorded_power
spec=importlib.util.spec_from_file_location('pulse_io',Path(__file__).with_name('pack-millstone-pulses.py'));io=importlib.util.module_from_spec(spec);spec.loader.exec_module(io)
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--start',type=int,required=True);p.add_argument('--scan-seconds',type=int,default=3600);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True);(a.output/'inputs').mkdir(exist_ok=True)
def hour(t):return dt.datetime.fromtimestamp(t,dt.timezone.utc).strftime('%Y-%m-%dT%H-00-00')
hours=sorted({hour(t)for t in range(a.start-3600,a.start+a.scan_seconds+3600,3600)})
ants=[];power=[]
for name in hours:
 for f in sorted((a.root/'metadata/antenna_control_metadata'/name).glob('*.h5')):
  with h5py.File(f)as h:
   for k,g in h.items():ants.append((int(k),{n:g[n][()]for n in g}))
 for f in sorted((a.root/'metadata/powermeter'/name).glob('*.h5')):
  with h5py.File(f)as h:
   for k,g in h.items():
    if 'zenith_power'in g:power.append((float(k),float(g['zenith_power'][()]),str(f)))
ants.sort();power.sort();anttimes=[x[0]for x in ants];ptimes=np.array([x[0]for x in power]);pvalues=np.array([x[1]for x in power])
cache=collections.OrderedDict();catalogues={};inputs=[];counts=collections.Counter()
def antenna(k):
 i=bisect.bisect_right(anttimes,k)-1
 return i>=0 and ants[i][1]['rx_antenna']==b'ZENITH' and ants[i][1]['tx_antenna']==b'ZENITH'
def watts(k):
 return recorded_power(k/1e6,ptimes,pvalues)
def rf(channel,second):
 key=(channel,second)
 if key not in cache:
  path=a.root/'rf_data'/channel/hour(second)/f'rf@{second}.000.h5'
  base,z=io.rf(path);cache[key]=(base,z,str(path),hashlib.sha256(path.read_bytes()).hexdigest())
  while len(cache)>24:cache.popitem(last=False)
 cache.move_to_end(key);return cache[key]
def vector(channel,k):
 second=k//1000000;base,z,path,sha=rf(channel,second);offset=k-base;result=z[offset:offset+10000];sources={path:sha}
 if len(result)<10000:
  b,nextz,path2,sha2=rf(channel,second+1);assert b==base+len(z),'RF discontinuity at second boundary';result=np.concatenate([result,nextz[:10000-len(result)]]);sources[path2]=sha2
 assert len(result)==10000;return result,sources
for begin in range(a.start,a.start+a.scan_seconds,10):
 ids={}
 for second in range(begin-1,begin+11):
  path=a.root/'metadata/id_metadata'/hour(second)/f'md@{second}.h5'
  if path.exists():
   with h5py.File(path)as h:ids.update({int(k):int(g['sweepid'][()])for k,g in h.items()if 'sweepid'in g})
 for family in ['coded','long']:
  valid=lambda k: (1<=ids[k]<=32 if family=='coded'else ids[k]==300) and antenna(k) and watts(k)>=400e3
  keys=sorted(k for k in ids if valid(k));paired={}
  for i,k in enumerate(keys):
   for j in [i-1,i+1]:
    if 0<=j<len(keys)and ids[keys[j]]==ids[k]and abs(keys[j]-k)<40000:paired[k]=keys[j];break
  selected=[k for k in keys if begin*1000000<=k<(begin+10)*1000000 and k in paired]
  if len(selected)<16:continue
  # Real partners just outside the integration are stored for clutter subtraction,
  # but are not themselves additional science measurements.
  keys=sorted(set(selected)|{paired[k]for k in selected});index={k:i for i,k in enumerate(keys)}
  for k in keys:
   if k not in selected:
    partner=next(x for x in selected if paired[x]==k);paired[k]=partner
  out=a.output/'inputs'/f'{begin}-ZENITH-{family}.h5'
  assert not out.exists(),'Choose a fresh campaign input version'
  tx_gate=[76,624]if family=='coded'else[76,645];echo_gate=[1000,8200]if family=='coded'else[1000,7700];noise_gate=[7700,8200]if family=='coded'else[7200,7700];injection_gate=[8400,8850]if family=='coded'else[7800,8371]
  echo=[];templates=[];delays=[];sources={};detections=[];cfar=[]
  for i,k in enumerate(keys):
   y,s=vector('zenith-l',k);sources.update(s);tx,s=vector('tx-h',k);sources.update(s)
   delay=io.channel_delay(y,tx,tx_gate);delays.append(delay);echo.append(y)
   template=np.zeros_like(y);template[tx_gate[0]:tx_gate[1]]=y[tx_gate[0]:tx_gate[1]];templates.append(template)
   second=k//1000000;path=a.root/'metadata/satellite_detections'/hour(second)/f'satellite@{second//3600*3600}.h5'
   if path not in catalogues:catalogues[path]=h5py.File(path,'r')
   cat=catalogues[path]
   if str(k)in cat:
    g=cat[str(k)]
    if g['channel'][()].decode()=='zenith-l':
     for j,(raw,dop)in enumerate(zip(g['raw_delay_sample'][()],g['doppler_hz'][()])):
      detections.append([i,float(raw)-tx_gate[0]-delay,float(dop)])
      cfar.append(dict(pulse_sample=k,raw_delay_sample=int(raw),doppler_hz=float(dop),snr_db=float(g['snr_db'][()][j]),cfar_ratio=float(g['cfar_ratio'][()][j]),cfar_threshold=float(g['cfar_threshold'][()][j])))
  lags=list(range(10,451,10))
  cfg=dict(station='Millstone Hill',antenna='ZENITH',channel='zenith-l',mode='alternating-code-16'if family=='coded'else'long-pulse-480us',radar_frequency_hz=440.2e6,start_unix_s=min(selected)/1e6,end_unix_s=(max(selected)+10000)/1e6,lpi=dict(code_energy_samples=[480.]*len(keys),science_pulse=[k in selected for k in keys],sample_rate_hz=1e6,decimation=30,lag_samples=lags,minimum_range_index=[int(np.rint((800-111-480*.2+lag)/30))for lag in lags],range_start=0,range_end=240,measurement_start=27,measurement_end=273,tx_gate=tx_gate,echo_gate=echo_gate,noise_gate=noise_gate,injection_gate=injection_gate,filter_half=100,passband_hz=18000.,outlier_threshold=10.,background_column=True),cancel_debris=True,cancellation=dict(refine_delay_samples=4.,refine_doppler_hz=100.),injection_temperature_k=1172.,range_doppler_fft=1024,range_doppler_samples=[],instrument_metadata=dict(transmitter_w=float(np.mean([watts(k)for k in selected])),pulse_transmitter_w=[watts(k)for k in keys],bit_length_us=30. if family=='coded'else 480.,code_bits=16 if family=='coded'else 1,dc_offset=[-.212,-.221],tx_to_receiver_delay_samples=delays,cfar=cfar,azimuth_deg=0.,elevation_deg=90.,science_pulse_count=len(selected),partner_pulse_count=len(keys)-len(selected),source_rf_sha256=sources,settings_reference='historical 10 s decode; 30 us range basis; 45 measured lags; full 7200 us delay grid; actual waveform-specific gates'))
  with h5py.File(out,'w')as h:
   h['schema_utf8']=np.frombuffer(b'ruisdap-pulses-v1',np.uint8);h['configuration_json']=np.frombuffer(json.dumps(cfg).encode(),np.uint8);h['pulse_sample']=np.array(keys,np.uint64);h['sweep_id']=np.array([ids[k]for k in keys],np.uint16);h['pair_index']=np.array([index[paired[k]]for k in keys],np.uint64);h['science_pulse']=np.array([k in selected for k in keys],np.uint8)
   for label,z in [('echo',echo),('transmit',templates)]:
    z=np.array(z);h.create_dataset(label,data=np.stack([z.real,z.imag],axis=-1).astype(np.int16),compression='gzip',compression_opts=1)
   h['detections']=np.array(detections,float).reshape(-1,3);h.attrs['rf_sha256_json']=json.dumps(sources);h.attrs['catalogue_paths_json']=json.dumps([str(x)for x in catalogues]);h.attrs['scope']='all supported observed powered ZENITH coded and 480 us long pulses; 45 measured lags; full historical range grid'
  inputs.append(str(out));counts[family]+=len(selected);print(json.dumps(dict(start=begin,family=family,science_pulses=len(selected),partners=len(keys)-len(selected),detections=len(detections))),flush=True)
for h in catalogues.values():h.close()
(a.output/'manifest.json').write_text(json.dumps(dict(inputs=inputs,output=str(a.output/'products'),live_directory=str(a.output/'live'),pulse_counts=counts)))
print(json.dumps(dict(inputs=len(inputs),science_pulses=counts)),flush=True)
