#!/usr/bin/env python
"""I/O-only bridge from explicit DigitalRF HDF5 files and a detection catalogue.
All voltage cancellation, spectrum estimation and lag inversion run in Rust.
No imports from the legacy processing chain; no source files are modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.signal import correlate,correlation_lags


def rf(path):
    with h5py.File(path) as h:
        idx=h['rf_data_index'][()]
        if idx.shape!=(1,2):raise ValueError('only continuous one-block DigitalRF files supported; provide a contiguous extraction for gapped data')
        a=h['rf_data'][()]
        if a.dtype.names!=('r','i'):raise ValueError('expected complex int16 DigitalRF data')
        return int(idx[0,0]),a['r'].ravel().astype(float)+1j*a['i'].ravel().astype(float)


def channel_delay(echo,tx,gate):
    # Estimate tx-h -> same-channel leakage delay from the full measured pulse,
    # including edge guards. Known catalogue delay must be converted to this frame.
    a=echo[gate[0]-40:gate[1]+80]; b=tx[gate[0]-40:gate[1]+80]
    c=correlate(a,b,mode='full',method='fft');lags=correlation_lags(len(a),len(b));use=np.abs(lags)<=40
    k=np.flatnonzero(use)[np.argmax(np.abs(c[use]))]
    delta=0.
    if 0<k<len(c)-1:
        z=np.abs(c[k-1:k+2]);den=z[0]-2*z[1]+z[2]
        if den!=0:delta=0.5*(z[0]-z[2])/den
    return float(lags[k]+delta)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--receive',required=True);p.add_argument('--transmit',required=True);p.add_argument('--catalogue',required=True);p.add_argument('--channel',required=True);p.add_argument('--output',required=True)
    p.add_argument('--pulse-key',type=int,action='append',required=True)
    a=p.parse_args();rxstart,rx=rf(a.receive);txstart,tx=rf(a.transmit)
    if rxstart!=txstart:raise ValueError('transmit/receive files have different absolute time origins')
    records=[]
    with h5py.File(a.catalogue) as h:
        for key in a.pulse_key:
            g=h[str(key)];channel=g['channel'][()].decode()
            if channel!=a.channel:raise ValueError(f'catalogue receiver {channel} does not match {a.channel}')
            codes=g['sweep_id'][()];code=int(codes[0])
            if not np.all(codes==code):raise ValueError('mixed modes in one pulse')
            if code==800:n=40000;gate=[69,2171];echo=[3721,30000];noise=[29500,30000];injection=[30176,32033]
            elif code==300:n=10000;gate=[76,645];echo=[1000,7700];noise=[7200,7700];injection=[7850,8371]
            elif 1<=code<=32:n=10000;gate=[76,624];echo=[1000,8200];noise=[7700,8200];injection=[8400,8850]
            else:raise ValueError(f'unsupported sweep code {code}')
            start=key-rxstart
            if start<0 or start+n>len(rx):raise ValueError('pulse crosses file boundary; provide a longer contiguous source extraction')
            y=rx[start:start+n];t=tx[start:start+n]
            delay=channel_delay(y,t,gate)
            # Same-channel measured pulse already contains receiver phase and
            # bit-edge response. Keep its shape, zero only outside guarded gate.
            template=np.zeros(n,complex);template[gate[0]:gate[1]]=y[gate[0]:gate[1]]
            detections=[[float(raw)-gate[0]-delay,float(dop)] for raw,dop in zip(g['raw_delay_sample'][()],g['doppler_hz'][()])]
            records.append(dict(key=key,code=code,n=n,tx_gate=gate,echo_gate=echo,noise_gate=noise,injection_gate=injection,y=y,template=template,delay=delay,detections=detections))
    if len({r['n'] for r in records})!=1:raise ValueError('pack only modes with identical pulse geometry together')
    # Cancellation inspection packs contain exact detected pulses. A worker LPI
    # pack must additionally supply real paired pulses and at least 16 records;
    # never fabricate a clutter partner to make an inversion run.
    r=records[0];nr=r['n'];npulse=len(records)
    config=dict(station='Millstone Hill',antenna='MISA' if a.channel.startswith('misa')else'ZENITH',channel=a.channel,mode=str(r['code']),radar_frequency_hz=440.2e6,start_unix_s=min(x['key']for x in records)/1e6,end_unix_s=(max(x['key']for x in records)+nr)/1e6,
      lpi=dict(sample_rate_hz=1e6,decimation=30,lag_samples=[0,30,60,90],range_start=34,range_end=80,measurement_start=34,measurement_end=min(nr//30-5,200),tx_gate=r['tx_gate'],echo_gate=r['echo_gate'],noise_gate=r['noise_gate'],injection_gate=r['injection_gate'],filter_half=100,passband_hz=18000.,outlier_threshold=10.,background_column=True),cancel_debris=True,cancellation=dict(refine_delay_samples=4.,refine_doppler_hz=100.),injection_temperature_k=1172.,range_doppler_fft=4096,range_doppler_samples=[])
    destination=Path(a.output)
    if destination.exists():raise ValueError('output exists; choose a new version')
    with h5py.File(destination,'w') as h:
        h['schema_utf8']=np.frombuffer(b'ruisdap-pulses-v1',np.uint8);h['configuration_json']=np.frombuffer(json.dumps(config).encode(),np.uint8)
        h.attrs['purpose']='cancellation inspection only; no ground-clutter partner supplied'
        h['pulse_sample']=np.array([r['key']for r in records],np.uint64)
        h['pair_index']=np.roll(np.arange(npulse),1).astype(np.uint64)
        for name,key in [('echo','y'),('transmit','template')]:
            z=np.array([r[key]for r in records]);h.create_dataset(name,data=np.stack([z.real,z.imag],axis=-1),compression='gzip')
        h['detections']=np.array([[i,*d]for i,r in enumerate(records)for d in r['detections']],float).reshape(-1,3)
        h['tx_to_receiver_delay_samples']=np.array([r['delay']for r in records])
        sources=h.create_group('sources')
        for label,path in [('receive',a.receive),('transmit',a.transmit),('catalogue',a.catalogue)]:
            sources.attrs[label+'_path']=str(Path(path).resolve());sources.attrs[label+'_sha256']=hashlib.sha256(Path(path).read_bytes()).hexdigest()
    print(json.dumps(dict(pulses=npulse,channel=a.channel,delay_samples=[r['delay']for r in records],output=str(destination))))

if __name__=='__main__':main()
