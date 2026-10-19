//! Regression for the HDF5 fit boundary and per-measurement power correction.
use guisdap_millstone::{fitting,signal,C};
use guisdap_ion::spectrum::{SpectrumModel,IonSpecies,PlasmaState};
use std::fs;
fn bytes(g:&hdf5::Group,name:&str,v:&[u8]){g.new_dataset_builder().with_data(v).create(name).unwrap();}
fn array(g:&hdf5::Group,name:&str,shape:&[usize],v:&[f64]){g.new_dataset::<f64>().shape(shape).create(name).unwrap().write_raw(v).unwrap();}
#[test]
fn physical_fit_survives_different_powers_and_pulse_normalizations(){
 let tmp=tempfile::tempdir().unwrap();let gain=6.2e11;let range:f64=300000.;let truth=PlasmaState{electron_density_m3:4e11,ion_temperature_k:900.,electron_temperature_k:1800.,velocity_m_s:80.,ion_collision_s:0.};
 let ions=vec![IonSpecies{mass_amu:16.,fraction:1.}];let model=SpectrumModel{radar_frequency_hz:440.2e6,ions:ions.clone(),electron_collision_s:0.};
 let f:Vec<_>=(-400..=400).map(|i|i as f64*250.).collect();let s=model.evaluate(truth,&f).unwrap();
 let mut impulse=vec![C::default();401];impulse[200]=C::new(1.,0.);let h=signal::lowpass(&impulse,1e6,21600.,100).unwrap();
 let response:Vec<_>=f.iter().map(|f|h.iter().enumerate().map(|(j,h)|*h*C::from_polar(1.,-2.*std::f64::consts::PI*f*j as f64/1e6)).sum::<C>().norm_sqr()).collect();
 let lags:Vec<_>=(1..=16).map(|i|i as f64*1e-5).collect();let mut paths=Vec::new();
 for (i,pt,energy)in [(0,600000.,480.),(1,1200000.,1600.)]{
  let path=tmp.path().join(format!("input-{i}.h5"));let h=hdf5::File::create(&path).unwrap();let r=h.create_group("records").unwrap().create_group("record").unwrap();
  let cfg=serde_json::json!({"station":"Millstone Hill","antenna":"ZENITH","channel":"zenith-l","mode":"test","radar_frequency_hz":440.2e6,"start_unix_s":1000.+i as f64,"end_unix_s":1001.+i as f64,"lpi":{"sample_rate_hz":1e6,"decimation":30,"code_energy_samples":[energy],"lag_samples":(10..=160).step_by(10).collect::<Vec<_>>(),"range_start":0,"range_end":1,"measurement_start":0,"measurement_end":1,"tx_gate":[0,480],"echo_gate":[500,900],"noise_gate":[900,950],"injection_gate":[950,1000],"filter_half":100,"passband_hz":18000.},"cancel_debris":true,"injection_temperature_k":1172.,"range_doppler_fft":1024,"range_doppler_samples":[],"instrument_metadata":{"transmitter_w":pt}});
  bytes(&r,"configuration_json",&serde_json::to_vec(&cfg).unwrap());array(&r,"range_m",&[1],&[range]);let g=r.create_group("ion_line").unwrap();array(&g,"lag_s",&[16],&lags);
  let values:Vec<_>=lags.iter().map(|lag|f.iter().zip(&s).zip(&response).enumerate().map(|(j,((&f,&s),&h))|C::from_polar(pt/(gain*(range/1000.).powi(2))*250.*h*truth.electron_density_m3*s*if j==0||j==800{0.5}else{1.},-2.*std::f64::consts::PI*f*lag)).sum::<C>()).flat_map(|z|[z.re,z.im]).collect();array(&g,"acf",&[1,16,2],&values);array(&g,"variance",&[1,16],&vec![(pt/1e6_f64*0.01).powi(2);16]);paths.push(path);
 }
 let manifest=serde_json::json!({"inputs":paths,"output":tmp.path().join("out"),"integration_s":180.,"magic_constant_by_antenna":{"ZENITH":gain},"ions":ions,"range_km":[290.,310.],"initial":{"electron_density_m3":2e11,"ion_temperature_k":1000.,"electron_temperature_k":1500.,"velocity_m_s":0.,"ion_collision_s":0.}});
 let m:fitting::Manifest=serde_json::from_value(manifest).unwrap();let output=tmp.path().join("fit.h5");fitting::run(&paths,&output,&m).unwrap();let h=hdf5::File::open(output).unwrap();let p=h.dataset("records/record/ion_fit/parameters").unwrap().read_raw::<f64>().unwrap();
 for (got,want)in p.iter().zip([4e11,900.,1800.,80.,0.]){assert!((got-want).abs()<want.abs().max(1.)*0.002,"{got} != {want}");}
 let e=h.dataset("records/record/ion_fit/standard_error").unwrap().read_raw::<f64>().unwrap();assert!(e[..4].iter().all(|e|e.is_finite()&&*e>0.));
 let _=fs::metadata(tmp.path()).unwrap();
}
