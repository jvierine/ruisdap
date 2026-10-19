#!/usr/bin/env python
"""Compare native spectra to a pinned, external Kudeki–Milla implementation.
Reference implementation is downloaded separately; it is never used at runtime.
"""
import argparse,subprocess,sys,json
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import atomic_mass,m_p
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--output',required=True);a=p.parse_args()
sys.path.insert(0,a.reference)
from ISRSpectrum import Specinit
cases=[];expected=[];freq=np.linspace(-40000,40000,801)
for radar in [440.2e6,930e6]:
 for ne in [1e10,1e11,1e12]:
  for ti,te in [(500.,500.),(900.,1800.),(1200.,3600.)]:
   for fraction in [0.,0.7,1.]:
    for velocity in [-300.,0.,400.]:
     model=dict(radar_frequency_hz=radar,ions=[dict(mass_amu=16.,fraction=fraction),dict(mass_amu=30.5,fraction=1-fraction)],electron_collision_s=0.)
     state=dict(electron_density_m3=ne,ion_temperature_k=ti,electron_temperature_k=te,velocity_m_s=velocity,ion_collision_s=0.)
     # ISRSpectrum expresses ion mass in proton masses; RUISDAP uses atomic
     # mass units. Convert explicitly before comparing absolute spectral power.
     ions=[[ne*fraction,ti,velocity,1.,16.*atomic_mass/m_p,0.,0.,0.],[ne*(1-fraction),ti,velocity,1.,30.5*atomic_mass/m_p,0.,0.,0.],[ne,te,velocity,-1.,1.,0.,0.,0.]]
     ref=Specinit(centerFrequency=radar,bMag=0.,f=freq).getspec(np.array(ions),alphadeg=90.)[1]/ne
     cases.append(dict(model=model,plasma=state,frequency_hz=freq.tolist()));expected.append(ref)
output=Path(a.output);output.parent.mkdir(parents=True,exist_ok=True)
config=output.with_suffix('.input.json');config.write_text(json.dumps(cases))
proc=subprocess.run(['target/release/examples/theory_grid',str(config)],check=True,capture_output=True,text=True)
actual=np.array([r['values']for r in json.loads(proc.stdout)]);expected=np.array(expected)
relative=np.max(abs(actual-expected),axis=1)/np.max(expected,axis=1)
with h5py.File(output,'w')as h:
 h.attrs['reference_git_revision']=subprocess.check_output(['git','-C',a.reference,'rev-parse','HEAD'],text=True).strip()
 h.attrs['domain']='unmagnetized, collisionless Maxwellian ion-line, singly charged ions; no near-perpendicular magnetic validation'
 h['frequency_hz']=freq;h['rust_spectrum']=actual;h['reference_spectrum']=expected;h['relative_peak_error']=relative
assert np.max(relative)<2e-6, np.max(relative)
print(json.dumps(dict(cases=len(cases),max_relative_peak_error=float(np.max(relative)),median_relative_peak_error=float(np.median(relative)))))
