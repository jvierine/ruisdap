#!/usr/bin/env python
"""I/O bridge for a bounded, paired pulse train; DSP is delegated to Rust.
Requires explicit receiver, mode timing, and continuous DigitalRF input files.
"""
import argparse,importlib.util,json,hashlib
from pathlib import Path
import h5py,numpy as np
module=importlib.util.spec_from_file_location('pack',Path(__file__).with_name('pack-millstone-pulses.py'));pack=importlib.util.module_from_spec(module);module.loader.exec_module(pack)
p=argparse.ArgumentParser();p.add_argument('--receive',required=True);p.add_argument('--id-metadata',required=True);p.add_argument('--catalogue');p.add_argument('--output',required=True);p.add_argument('--channel',required=True);p.add_argument('--cancel',action='store_true');a=p.parse_args()
start,rx=pack.rf(a.receive)
with h5py.File(a.id_metadata)as h:
 ids={int(k):int(g['sweepid'][()])for k,g in h.items()if 'sweepid'in g}
keys=sorted(k for k,code in ids.items()if 1<=code<=32 and start<=k and k+10000<=start+len(rx))
# Adjacent same-code pairs are the Turunen clutter subtraction. Exclude pulses
# without a partner; do not pair unrelated codes or invent synthetic records.
paired={}
for i,k in enumerate(keys):
 for j in [i-1,i+1]:
  if 0<=j<len(keys)and ids[keys[j]]==ids[k]and abs(keys[j]-k)<20000:paired[k]=keys[j];break
keys=sorted(paired)
if len(keys)<16:raise ValueError('at least 16 real paired coded pulses are required')
y=np.array([rx[k-start:k-start+10000]for k in keys]);tx=np.zeros_like(y);tx[:,76:624]=y[:,76:624]
# Phase-aligned same-code averaging reduces template noise without fitting one
# independent amplitude to each baud edge. Only one complex gain is fit per echo.
for code in sorted(set(ids[k]for k in keys)):
 idx=[i for i,k in enumerate(keys)if ids[k]==code];ref=tx[idx[0]];aligned=[]
 for i in idx:
  norm=np.sqrt(np.vdot(tx[i],tx[i]).real);z=tx[i]/norm;rotation=np.vdot(ref,z);aligned.append(z*np.exp(-1j*np.angle(rotation)))
 mean=np.mean(aligned,axis=0)
 for i in idx:tx[i]=mean
config=dict(station='Millstone Hill',antenna='ZENITH'if a.channel.startswith('zenith')else'MISA',channel=a.channel,mode='alternating-code-16',radar_frequency_hz=440.2e6,start_unix_s=keys[0]/1e6,end_unix_s=(keys[-1]+10000)/1e6,
 lpi=dict(sample_rate_hz=1e6,decimation=30,lag_samples=[10,20,30,40],range_start=28,range_end=120,measurement_start=34,measurement_end=240,tx_gate=[76,624],echo_gate=[1000,8200],noise_gate=[7700,8200],injection_gate=[8400,8850],filter_half=100,passband_hz=18000.,outlier_threshold=10.,background_column=True),cancel_debris=a.cancel,cancellation=dict(refine_delay_samples=4.,refine_doppler_hz=100.),injection_temperature_k=1172.,range_doppler_fft=1024,range_doppler_samples=[])
detections=[]
if a.catalogue:
 with h5py.File(a.catalogue)as h:
  for i,k in enumerate(keys):
   if str(k)not in h:continue
   g=h[str(k)]
   if g['channel'][()].decode()!=a.channel:continue
   # The tx-to-receiver group delay is 11 samples in this catalogue's original
   # frame. Treat it as an initial value, allow bounded native refinement, and
   # retain it in metadata. A production recipe should use its measured delay.
   for raw,dop in zip(g['raw_delay_sample'][()],g['doppler_hz'][()]):detections.append([i,float(raw)-76.-11.,float(dop)])
if Path(a.output).exists():raise ValueError('output exists')
with h5py.File(a.output,'w')as h:
 h['schema_utf8']=np.frombuffer(b'ruisdap-pulses-v1',np.uint8);h['configuration_json']=np.frombuffer(json.dumps(config).encode(),np.uint8)
 h['pulse_sample']=np.array(keys,np.uint64);h['sweep_id']=np.array([ids[k]for k in keys],np.uint16);h['pair_index']=np.array([keys.index(paired[k])for k in keys],np.uint64)
 h['echo']=np.stack([y.real,y.imag],axis=-1);h['transmit']=np.stack([tx.real,tx.imag],axis=-1);h['detections']=np.array(detections,float).reshape(-1,3)
 h.attrs['scope']='bounded numerical migration pilot; antenna switching/power metadata must be validated before campaign production'
 h.attrs['receive_sha256']=hashlib.sha256(Path(a.receive).read_bytes()).hexdigest();h.attrs['id_metadata_sha256']=hashlib.sha256(Path(a.id_metadata).read_bytes()).hexdigest()
print(json.dumps(dict(pulses=len(keys),detections=len(detections),output=a.output)))
