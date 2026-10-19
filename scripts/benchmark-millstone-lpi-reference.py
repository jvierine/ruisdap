#!/usr/bin/env python3
"""Timed independent sparse NumPy/SciPy replay of the current native LPI operator.
This is NOT a benchmark of the historical DigitalRF campaign Python controller.
Input read and fixed ADC offset subtraction excluded; FFT FIR included in kernel.
"""
import argparse,hashlib,json,time,platform
from pathlib import Path
import h5py,numpy as np,scipy
from scipy.signal import fftconvolve
from scipy import sparse,linalg

def estimate(echo,tx,pair,c):
 s=c['lpi'];N=echo.shape[1];P=len(echo);half=s['filter_half'];j=np.arange(2*half);m=j-half+1e-6
 fir=np.hanning(2*half)*np.sin(2*np.pi*1.2*s['passband_hz']/s['sample_rate_hz']*m)/(np.pi*m);fir/=fir.sum()
 y=fftconvolve(echo,fir[None,:],mode='full',axes=1)[:,half:half+N]
 no=s['noise_gate'];inj=s['injection_gate'];dc=np.median(y[:,no[0]:no[1]].mean(axis=1).real)+1j*np.median(y[:,no[0]:no[1]].mean(axis=1).imag)
 noise=np.median(np.mean(abs(y[:,no[0]:no[1]]-dc)**2,axis=1));injection=np.median(np.mean(abs(y[:,inj[0]:inj[1]]-dc)**2,axis=1));gain=(injection-noise)/c['injection_temperature_k']
 tx=tx.copy();tx[:,:s['tx_gate'][0]]=0;tx[:,s['tx_gate'][1]:]=0;tx/=np.sqrt(np.sum(abs(tx)**2,axis=1))[:,None]
 dec=s['decimation'];nr=s['range_end']-s['range_start'];nl=len(s['lag_samples']);nc=nr+int(s['background_column']);result=np.full((nr,nl),np.nan+1j*np.nan);resultg=result.copy();variance=np.full((nr,nl),np.nan);retained=[]
 ms=np.arange(s['measurement_start'],s['measurement_end']);M=len(ms);ri=np.arange(s['range_start'],s['range_end']);indices=ms[:,None]-ri
 energy=np.asarray(s.get('code_energy_samples')or np.ones(P));energy=energy*(dec/s['sample_rate_hz']*299792.458/2) if s.get('code_energy_samples')else energy
 science=np.asarray(s.get('science_pulse')or [True]*P)
 for li,lag in enumerate(s['lag_samples']):
  def products(z):
   z=z[:,:N-lag]*np.conj(z[:,lag:]);n=z.shape[1]//dec;return z[:,:n*dec].reshape(P,n,dec).sum(axis=2)
  e=products(y);g=products(y-y[pair]);amb=products(tx)*energy[:,None]
  robust=np.sqrt(2*np.percentile(abs(e[:,ms])**2,34,axis=0));robustg=np.sqrt(2*np.percentile(abs(g[:,ms])**2,34,axis=0));local=np.sqrt(sum(np.roll(abs(e[:,ms])**2,j-5,axis=0)for j in range(10))/10);median=np.median(local);sigma=np.maximum(local,median)
  keep=(abs(e[:,ms])<=s['outlier_threshold']*robust)&(abs(g[:,ms])<=s['outlier_threshold']*robustg)&(sigma<=100*median)
  usable=(ms*dec>=s['echo_gate'][0]+half)&((ms+1)*dec-1+lag<s['echo_gate'][1]-half);keep&=usable;keep&=science[:,None]
  rows=[];columns=[];values=[];kept=np.zeros_like(keep);minimum=s.get('minimum_range_index',[s['range_start']]*nl)[li]
  for pi in np.flatnonzero(science):
   d=amb[pi,np.maximum(indices,0)].copy();d[(indices<0)|(ri[None,:]<minimum)]=0;d[abs(d)**2<=1e-30]=0
   if s['background_column']:d=np.column_stack([d,np.ones(M,complex)])
   selected=keep[pi]&np.any(d!=0,axis=1);kept[pi]=selected
   mr,cr=np.nonzero(d[selected]);rows.append(mr+int(kept[:pi].sum()));columns.append(cr);values.append(d[selected][mr,cr])
  retained.append(float(kept.sum()/(P*M)))
  if kept.sum()<nc:continue
  A=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(columns))),shape=(int(kept.sum()),nc)).tocsr();A=A.multiply(1/sigma[kept,None]).tocsr()
  scale=np.sqrt(np.asarray(abs(A).power(2).sum(axis=0)).ravel());active=scale>0;B=A[:,active].multiply(1/scale[active]).tocsr()
  normal=(B.conj().T@B).toarray();rhs=B.conj().T@np.column_stack([e[:,ms][kept]/sigma[kept],g[:,ms][kept]/sigma[kept]])
  try:factor=linalg.cho_factor(normal,lower=True);x=linalg.cho_solve(factor,rhs)/scale[active,None];vv=np.diag(linalg.cho_solve(factor,np.eye(len(normal)))).real/scale[active]**2
  except linalg.LinAlgError:continue
  allx=np.full((nc,2),np.nan+1j*np.nan);allvar=np.full(nc,np.nan);allx[active]=x;allvar[active]=vv
  result[:,li]=allx[:nr,0]/gain;resultg[:,li]=allx[:nr,1]/(2*gain);variance[:,li]=allvar[:nr]/gain**2
 return result,resultg,variance,np.array(retained),gain

def main():
 p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('native',type=Path);p.add_argument('report',type=Path);p.add_argument('--repeats',type=int,default=3);a=p.parse_args();assert not a.report.exists()
 with h5py.File(a.input)as h:
  c=json.loads(bytes(h['configuration_json'][:]));assert not c['cancel_debris']or h['detections'].shape[0]==0,'reference does not cancel echoes'
  v=h['echo'][:];echo=v[...,0]+1j*v[...,1];v=h['transmit'][:];tx=v[...,0]+1j*v[...,1];pair=h['pair_index'][:].astype(int)
 dc=c.get('instrument_metadata',{}).get('dc_offset',[0,0]);dc=complex(*dc);echo-=dc;tx[:,c['lpi']['tx_gate'][0]:c['lpi']['tx_gate'][1]]-=dc
 times=[]
 for rep in range(a.repeats):
  start=time.perf_counter();result=estimate(echo,tx,pair,c);times.append(time.perf_counter()-start);print(json.dumps(dict(repeat=rep,seconds=times[-1])),flush=True)
 with h5py.File(a.native)as h:
  r=next(iter(h['records'].values()))['ion_line'];v=r['acf'][:];native=v[...,0]+1j*v[...,1];v=r['ground_clutter_acf'][:];ground=v[...,0]+1j*v[...,1]
  np.testing.assert_allclose(native,result[0],rtol=2e-8,atol=2e-9,equal_nan=True);np.testing.assert_allclose(ground,result[1],rtol=2e-8,atol=2e-9,equal_nan=True)
  np.testing.assert_allclose(r['variance'][:],result[2],rtol=3e-8,atol=1e-10,equal_nan=True);np.testing.assert_allclose(r['retained_fraction'][:],result[3],rtol=1e-12);np.testing.assert_allclose(r['receiver_gain'][0],result[4],rtol=1e-12)
 error=float(np.nanmax(abs(native-result[0]))/np.nanmax(abs(result[0])))
 with h5py.File(a.report,'w')as h:
  h['seconds']=times;h['reference_acf']=result[0];h['variance']=result[2];h['retained_fraction']=result[3];h.attrs['max_relative_peak_error']=error;h.attrs['finite_acf_cells']=int(np.isfinite(native).sum());h.attrs['input_sha256']=hashlib.sha256(a.input.read_bytes()).hexdigest();h.attrs['python']=platform.python_version();h.attrs['numpy']=np.__version__;h.attrs['scipy']=scipy.__version__
 print(json.dumps(dict(max_relative_peak_error=error,seconds=times)),flush=True)
if __name__=='__main__':main()
