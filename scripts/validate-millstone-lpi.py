#!/usr/bin/env python
"""Independent NumPy/SciPy WLS replay of a native bounded pulse-train product.
Checks actual signed ACF, complete variance, retained rows and receiver gain.
"""
import argparse,json
from pathlib import Path
import h5py,numpy as np
from scipy.signal import fftconvolve
p=argparse.ArgumentParser();p.add_argument('input');p.add_argument('output');p.add_argument('--report',required=True);a=p.parse_args()
with h5py.File(a.input)as h:
 c=json.loads(bytes(h['configuration_json'][()]));s=c['lpi'];v=h['echo'][()];echo=v[...,0]+1j*v[...,1];v=h['transmit'][()];tx=v[...,0]+1j*v[...,1];pair=h['pair_index'][()].astype(int)
 if c['cancel_debris']:raise ValueError('operator parity requires a cancellation-disabled control')
N=echo.shape[1];P=len(echo);half=s['filter_half'];j=np.arange(2*half);m=j-half+1e-6
fir=np.hanning(2*half)*np.sin(2*np.pi*1.2*s['passband_hz']/s['sample_rate_hz']*m)/(np.pi*m);fir/=fir.sum()
# Linear convolution, matching physical sample origin, independently via FFT.
y=np.array([fftconvolve(x,fir,mode='full')[half:half+N]for x in echo])
no=s['noise_gate'];inj=s['injection_gate'];dc=np.median(y[:,no[0]:no[1]].mean(axis=1).real)+1j*np.median(y[:,no[0]:no[1]].mean(axis=1).imag)
noise=np.median(np.mean(abs(y[:,no[0]:no[1]]-dc)**2,axis=1));injection=np.median(np.mean(abs(y[:,inj[0]:inj[1]]-dc)**2,axis=1));gain=(injection-noise)/c['injection_temperature_k']
tx[:,:s['tx_gate'][0]]=0;tx[:,s['tx_gate'][1]:]=0;tx/=np.sqrt(np.sum(abs(tx)**2,axis=1))[:,None]
dec=s['decimation'];nr=s['range_end']-s['range_start'];nl=len(s['lag_samples']);result=np.full((nr,nl),np.nan+1j*np.nan);resultg=result.copy();variance=np.full((nr,nl),np.nan);retained=[]
for li,lag in enumerate(s['lag_samples']):
 def products(z):
  z=z[:,:N-lag]*np.conj(z[:,lag:]);n=z.shape[1]//dec;return z[:,:n*dec].reshape(P,n,dec).sum(axis=2)
 e=products(y);g=products(y-y[pair]);amb=products(tx);ms=np.arange(s['measurement_start'],s['measurement_end']);M=len(ms)
 robust=np.sqrt(2*np.percentile(abs(e[:,ms])**2,34,axis=0));robustg=np.sqrt(2*np.percentile(abs(g[:,ms])**2,34,axis=0));local=np.sqrt(sum(np.roll(abs(e[:,ms])**2,j-5,axis=0)for j in range(10))/10);median=np.median(local);sigma=np.maximum(local,median)
 keep=(abs(e[:,ms])<=s['outlier_threshold']*robust)&(abs(g[:,ms])<=s['outlier_threshold']*robustg)&(sigma<=100*median)
 usable=(ms*dec>=s['echo_gate'][0]+half)&((ms+1)*dec-1+lag<s['echo_gate'][1]-half);keep&=usable
 ri=np.arange(s['range_start'],s['range_end']);indices=ms[:,None]-ri;design=amb[:,indices]
 design[:,indices<0]=0
 if s['background_column']:design=np.concatenate([design,np.ones((P,M,1),complex)],axis=2)
 rows=design[keep];obs=e[:,ms][keep];obsg=g[:,ms][keep];sig=sigma[keep];retained.append(float(keep.sum()/(P*M)))
 A=rows/sig[:,None];b=obs/sig;bg=obsg/sig;scale=np.linalg.norm(A,axis=0);active=scale>0;B=A[:,active]/scale[active]
 try:
  normal=B.conj().T@B;cov=np.linalg.inv(normal);x=np.linalg.solve(normal,B.conj().T@b)/scale[active];xg=np.linalg.solve(normal,B.conj().T@bg)/scale[active];vv=np.diag(cov).real/scale[active]**2
 except np.linalg.LinAlgError:continue
 allx=np.full(len(scale),np.nan+1j*np.nan);allg=allx.copy();allvar=np.full(len(scale),np.nan);allx[active]=x;allg[active]=xg;allvar[active]=vv
 result[:,li]=allx[:nr]/gain;resultg[:,li]=allg[:nr]/(2*gain);variance[:,li]=allvar[:nr]/gain**2
with h5py.File(a.output)as h:
 r=next(iter(h['records'].values()));actual=r['ion_line/acf'][()];actual=actual[...,0]+1j*actual[...,1];gc=r['ion_line/ground_clutter_acf'][()];gc=gc[...,0]+1j*gc[...,1]
 assert np.isfinite(actual).sum()>100,'an all-missing product is not a successful inversion'
 np.testing.assert_allclose(actual,result,rtol=2e-9,atol=2e-7,equal_nan=True);np.testing.assert_allclose(gc,resultg,rtol=2e-9,atol=2e-7,equal_nan=True);np.testing.assert_allclose(r['ion_line/variance'][()],variance,rtol=3e-9,atol=1e-6,equal_nan=True);np.testing.assert_allclose(r['ion_line/retained_fraction'][()],retained,rtol=1e-12);np.testing.assert_allclose(r['ion_line/receiver_gain'][0],gain,rtol=1e-12)
 error=float(np.nanmax(abs(actual-result))/np.nanmax(abs(result)))
with h5py.File(a.report,'w')as h:
 h['reference_acf']=result;h['reference_variance']=variance;h['retained_fraction']=retained;h.attrs['max_relative_peak_error']=error;h.attrs['pulse_count']=P;h.attrs['finite_acf_cells']=int(np.isfinite(actual).sum())
print(json.dumps(dict(pulses=P,finite_acf_cells=int(np.isfinite(actual).sum()),max_relative_peak_error=error,receiver_gain=gain,retained_fraction=retained)))
