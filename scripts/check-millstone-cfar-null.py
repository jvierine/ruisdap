#!/usr/bin/env python3
"""Run an external, unmodified CFAR detector on synthetic diffuse-only RF.

Original RF is never opened; synthetic voltages remain in memory. This probes
model-dependent false candidates, not the real-archive false-alarm rate.
"""
import argparse,hashlib,importlib.util,json,pathlib,sys,types
import h5py,numpy as np

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--detector-script',type=pathlib.Path,required=True)
 p.add_argument('--spectrum',type=pathlib.Path,required=True)
 p.add_argument('--output',type=pathlib.Path,required=True)
 p.add_argument('--trials',type=int,default=50)
 a=p.parse_args();assert not a.output.exists() and a.trials>0
 spec=importlib.util.spec_from_file_location('cfar_null_reference',a.detector_script)
 d=importlib.util.module_from_spec(spec);sys.modules[spec.name]=d;spec.loader.exec_module(d)
 with h5py.File(a.spectrum) as h:
  psd=h['ion_line_psd'][:];sr=float(h.attrs['sample_rate_hz']);physics=bytes(h['physics_json'][:]);assert sr==d.FS==1e6
 n=10000;fft=d.FFTBank(128);delay=round(550e3*2/d.C*sr)
 parameters=types.SimpleNamespace(decimation_factor=8,fft_padding=4.,max_doppler_hz=60000.,pfa=1e-10,max_echoes_per_pulse=32,receiver_delay_samples=0.)
 class InMemoryRF:
  def __init__(self,y,tx):self.values={'zenith-l':y,'tx-h':tx}
  def read_vector(self,sample,length,channel):return self.values[channel][:length]
 with h5py.File(a.output,'x') as report:
  report.attrs['complete']=False
  report.attrs['detector_sha256']=hashlib.sha256(a.detector_script.read_bytes()).hexdigest()
  report.attrs['spectrum_sha256']=hashlib.sha256(a.spectrum.read_bytes()).hexdigest()
  report.attrs['scope']='Synthetic Gaussian ion-line fluctuations in thin or Gaussian range-extended plasma; no coherent target in null cases. Approximate scattering environment, not an empirical real-archive false-alarm rate.'
  report['physics_json']=np.frombuffer(physics,np.uint8)
  report['parameters_json']=np.frombuffer(json.dumps(vars(parameters)).encode(),np.uint8)
  for family,sweep in [('coded',1),('long',300)]:
   mode=d.MODES[sweep];body=np.zeros(mode.tx1-mode.tx0,np.complex64)
   code=np.array([1,1,-1,1,-1,-1,1,-1,1,1,1,-1,-1,1,-1,-1],np.float32)
   body[:480]=1 if family=='long' else np.repeat(code,30)
   tx=np.zeros(n,np.complex64);tx[mode.tx0:mode.tx1]=body
   for width in [0,240]:
    offsets=np.array([0]) if width==0 else np.round(np.linspace(-3*width,3*width,49)).astype(int)
    weights=np.ones(1) if width==0 else np.exp(-.5*(offsets/width)**2);weights/=weights.sum()
    for kind,power in [('null',v)for v in [0.,.03,.1,.3,1.]]+[('positive',.1)]:
     rows=[];near=[];profile=np.zeros(n)
     for offset,w in zip(offsets,weights):
      start=delay+mode.tx0+int(offset);profile[start:start+len(body)]+=w*abs(body)**2
     for seed in range(1,a.trials+1):
      rng=np.random.default_rng(seed)
      y=((rng.standard_normal(n)+1j*rng.standard_normal(n))/np.sqrt(2)).astype(np.complex64)
      y[mode.tx0:mode.tx1]+=100*body
      if power:
       for offset,w in zip(offsets,weights):
        z=(rng.standard_normal(len(psd))+1j*rng.standard_normal(len(psd)))/np.sqrt(2)
        z=np.fft.ifft(z*np.sqrt(psd))*len(psd)/np.sqrt(psd.sum())
        start=delay+mode.tx0+int(offset)
        y[start:start+len(body)]+=body*z[np.arange(start,start+len(body))%len(z)]*np.sqrt(power*w)
      if kind=='positive':
       start=delay+mode.tx0;t=np.arange(start,start+len(body))
       y[start:start+len(body)]+=3*body*np.exp(2j*np.pi*7000*t/sr)
      detections=d.detect_pulse(InMemoryRF(y,tx),0,sweep,'zenith-l',parameters,fft)
      hits=[v for v in detections if abs(v['corrected_delay_samples']-delay)<max(16,3*width)]
      rows.append(detections);near.append(len(hits))
     g=report.create_group(f'{family}/width_{width}/{kind}_{power}')
     g.attrs['nominal_diffuse_power']=power
     g.attrs['maximum_expected_diffuse_receiver_ratio']=power*profile.max()
     g['near_layer_count']=near;g['detections_json']=np.frombuffer(json.dumps(rows).encode(),np.uint8)
     report.flush()
     print(f'{family} width={width} {kind} diffuse={power}: pulses detected near layer {sum(v>0 for v in near)}/{a.trials}; total candidates {sum(len(v)for v in rows)}',flush=True)
  report.attrs['complete']=True
 print('CFAR SYNTHETIC CHECK COMPLETE',a.output,flush=True)

if __name__=='__main__':main()
