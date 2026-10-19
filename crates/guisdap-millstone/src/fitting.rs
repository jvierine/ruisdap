//! Local homogeneous ion fits using the common, tested ISR spectrum model.
//! Calibrations and composition are supplied explicitly; no antenna defaults.
use anyhow::{ensure,Context,Result};
use guisdap_ion::{fit::{FitProblem,FitSettings,Measurement},spectrum::{IonSpecies,PlasmaState,SpectrumModel}};
use serde::{Deserialize,Serialize};
use std::{collections::BTreeMap,fs,path::{Path,PathBuf}};
use crate::{product::Configuration,signal,C};
#[derive(Deserialize)]
pub struct Manifest {
 pub inputs:Vec<PathBuf>, pub output:PathBuf,
 #[serde(default)]pub integration_s:f64,
 pub magic_constant_by_antenna:BTreeMap<String,f64>,
 pub ions:Vec<IonSpecies>, pub range_km:[f64;2],
 pub initial:PlasmaState,
 #[serde(default)]pub settings:FitSettings,
 #[serde(default)]pub historical_settings:bool,
 #[serde(default)]pub windows:Vec<[f64;2]>,
}
#[derive(Serialize)]
struct Gate {range_m:f64,#[serde(default)]ions:Vec<IonSpecies>,result:Option<guisdap_ion::fit::FitResult>,error:Option<String>}
fn bytes(g:&hdf5::Group,name:&str,v:&[u8])->Result<()>{g.new_dataset_builder().with_data(v).create(name)?;Ok(())}
fn array(g:&hdf5::Group,name:&str,shape:&[usize],v:&[f64],units:&str)->Result<()>{let d=g.new_dataset::<f64>().shape(shape).create(name)?;d.write_raw(v)?;d.new_attr_builder().with_data(units).create("units")?;Ok(())}
pub fn jobs(m:&Manifest)->Result<Vec<Vec<PathBuf>>>{
 let mut groups:BTreeMap<(String,i64),Vec<PathBuf>>=BTreeMap::new();
 for (i,input) in m.inputs.iter().enumerate(){let h=hdf5::File::open(input)?;let rs=h.group("records")?;let r=rs.group(&rs.member_names()?[0])?;let c:Configuration=serde_json::from_slice(&r.dataset("configuration_json")?.read_raw::<u8>()?)?;
  if !m.magic_constant_by_antenna.contains_key(&c.antenna){continue;}
  let bucket=if let Some(i)=m.windows.iter().position(|w|c.start_unix_s>=w[0]&&c.start_unix_s<w[1]){i as i64}else if !m.windows.is_empty(){m.windows.len()as i64+(c.start_unix_s/m.integration_s).floor()as i64}else if m.integration_s>0.{(c.start_unix_s/m.integration_s).floor()as i64}else{i as i64};groups.entry((c.antenna,bucket)).or_default().push(input.clone());
 }Ok(groups.into_values().collect())
}
pub fn run(inputs:&[PathBuf],output:&Path,m:&Manifest)->Result<()> {
 let input=inputs.first().context("empty fit integration")?;
 ensure!(!output.exists(),"fit product already exists");
 let h=hdf5::File::open(input)?;let records=h.group("records")?;let id=records.member_names()?.into_iter().next().context("no records")?;let record=records.group(&id)?;
 let mut config:Configuration=serde_json::from_slice(&record.dataset("configuration_json")?.read_raw::<u8>()?)?;
 let gain=*m.magic_constant_by_antenna.get(&config.antenna).context("no reviewed calibration for this antenna")?;
 ensure!(gain.is_finite()&&gain>0.&&!config.lpi.code_energy_samples.is_empty(),"fit requires positive response-normalized calibration and pulse energy metadata");
 let first_ion=record.group("ion_line")?;let lags=first_ion.dataset("lag_s")?.read_raw::<f64>()?;let range=record.dataset("range_m")?.read_raw::<f64>()?;
 let count=range.len()*lags.len();let mut weight=vec![0.;count];let mut sum=vec![C::default();count];let mut source_meta=Vec::new();
 let mut begin=f64::INFINITY;let mut end:f64=f64::NEG_INFINITY;
 for input in inputs{
  let h=hdf5::File::open(input)?;let rs=h.group("records")?;let r=rs.group(&rs.member_names()?[0])?;let c:Configuration=serde_json::from_slice(&r.dataset("configuration_json")?.read_raw::<u8>()?)?;
  ensure!(c.antenna==config.antenna&&c.channel==config.channel&&c.lpi.filter_half==config.lpi.filter_half&&c.lpi.passband_hz==config.lpi.passband_hz&&c.lpi.sample_rate_hz==config.lpi.sample_rate_hz&&!c.lpi.code_energy_samples.is_empty(),"incompatible measurement responses");
  let ion=r.group("ion_line")?;ensure!(ion.dataset("lag_s")?.read_raw::<f64>()?==lags&&r.dataset("range_m")?.read_raw::<f64>()?==range,"incompatible range or lag grid");
  let pt=c.instrument_metadata["transmitter_w"].as_f64().context("missing recorded transmit power")?;ensure!(pt.is_finite()&&pt>0.,"invalid recorded power");
  let mut a=ion.dataset("acf")?.read_raw::<f64>()?;
  if m.historical_settings{let ground=ion.dataset("ground_clutter_acf")?.read_raw::<f64>()?;for (ri,&r)in range.iter().enumerate(){if r<300000.{let lo=ri*lags.len()*2;a[lo..lo+lags.len()*2].copy_from_slice(&ground[lo..lo+lags.len()*2]);}}}
  let v=ion.dataset("variance")?.read_raw::<f64>()?;
  for i in 0..count{let z=C::new(a[2*i],a[2*i+1])/pt;let variance=v[i]/pt.powi(2);if z.re.is_finite()&&z.im.is_finite()&&variance.is_finite()&&variance>0.{let w=1./variance;weight[i]+=w;sum[i]+=z*w;}}
  begin=begin.min(c.start_unix_s);end=end.max(c.end_unix_s);
  use sha2::{Digest,Sha256};source_meta.push(serde_json::json!({"path":input,"sha256":format!("{:x}",Sha256::digest(fs::read(input)?)),"start_unix_s":c.start_unix_s,"end_unix_s":c.end_unix_s,"transmitter_w":pt,"azimuth_deg":c.instrument_metadata["azimuth_deg"],"elevation_deg":c.instrument_metadata["elevation_deg"],"science_pulse_count":c.instrument_metadata["science_pulse_count"],"mode":c.mode}));
 }
 let mut acf:Vec<_>=sum.iter().zip(&weight).map(|(s,w)|if *w>0.{*s / *w}else{C::new(f64::NAN,f64::NAN)}).flat_map(|z|[z.re,z.im]).collect();
 let mut var:Vec<_>=weight.iter().map(|w|if *w>0.{1./w}else{f64::NAN}).collect();
 if m.historical_settings {
  let original=acf.clone();let original_var=var.clone();
  for (ri,&r)in range.iter().enumerate(){let radius=if r<300000.{2}else if r<700000.{3}else{5};
   let lo=ri.saturating_sub(radius);let hi=(ri+radius).min(range.len());
   for li in 0..lags.len(){let mut z=C::default();let mut total=0.;let mut v=0.;
    for j in lo..hi{let i=j*lags.len()+li;let a=C::new(original[2*i],original[2*i+1]);if a.re.is_finite()&&a.im.is_finite()&&original_var[i].is_finite()&&original_var[i]>0.{let w=range[j].powi(2);z+=a*w;total+=w;v+=w*w*original_var[i];}}
    let i=ri*lags.len()+li;if total>0.{acf[2*i]=z.re/total;acf[2*i+1]=z.im/total;var[i]=v/total.powi(2);}else{acf[2*i]=f64::NAN;acf[2*i+1]=f64::NAN;var[i]=f64::NAN;}
   }
  }
 }
 config.start_unix_s=begin;config.end_unix_s=end;
 let pt=1.; // Each measurement's recorded power was removed before averaging.
 // Model the same finite, DC-normalized receive filter used on voltages.
 // The common spectral model integrates quadrature in Hz through this kernel.
 let (frequencies,response,step)=receiver_response(&config)?;
 let model=SpectrumModel{radar_frequency_hz:config.radar_frequency_hz,ions:m.ions.clone(),electron_collision_s:0.};model.validate()?;
 let mut gates=Vec::new();let mut warm=m.initial;
 for (ri,&r)in range.iter().enumerate(){if r/1000.<m.range_km[0]||r/1000.>m.range_km[1]{continue;}
  let mut gate_model=model.clone();
  if m.historical_settings{let h=r/1000.;let f=if h<120.{1.}else{let hscale=10.-6.*((-(h-120.)/40.).min(50.)).exp();let z=(-(h-180.)/hscale).min(50.);1.-2./(1.+(1.+8.*z.exp()).sqrt())};gate_model.ions=vec![IonSpecies{mass_amu:32.,fraction:f},IonSpecies{mass_amu:16.,fraction:1.-f}];}
  let scale=pt/(gain*(r/1000.).powi(2));
  let measurements:Vec<_>=lags.iter().enumerate().filter_map(|(li,&lag)|{let i=ri*lags.len()+li;let z=C::new(acf[2*i],acf[2*i+1]);let sigma=(var[i]/2.).sqrt();if !z.re.is_finite()||!z.im.is_finite()||!sigma.is_finite()||sigma<=0.{return None;}
   let kernel=frequencies.iter().zip(&response).enumerate().map(|(j,(&f,&h))|C::from_polar(scale*h*step*if j==0||j+1==frequencies.len(){0.5}else{1.},-2.*std::f64::consts::PI*f*lag)).collect();Some(Measurement{lag_s:lag,value:z,sigma:[sigma,sigma],kernel})}).collect();
  if measurements.len()<8{gates.push(Gate{range_m:r,ions:gate_model.ions.clone(),result:None,error:Some("fewer than eight measured lags".into())});continue;}
  let mut problem=FitProblem{model:gate_model,frequency_hz:frequencies.clone(),measurements,settings:m.settings.clone()};
  let mut starts=vec![warm];if m.historical_settings{let mut negative=warm;negative.velocity_m_s=-negative.velocity_m_s;starts.push(negative);starts.push(PlasmaState{electron_density_m3:m.initial.electron_density_m3,ion_temperature_k:500.,electron_temperature_k:550.,velocity_m_s:100.,ion_collision_s:0.});}
  let mut best:Option<guisdap_ion::fit::FitResult>=None;let mut best_ions=problem.model.ions.clone();
  let fractions=if m.historical_settings&&r>700000.{vec![0.,0.05,0.15,0.3,0.5]}else{vec![-1.]};
  for fraction in fractions{if fraction>=0.{problem.model.ions=vec![IonSpecies{mass_amu:16.,fraction:1.-fraction},IonSpecies{mass_amu:1.,fraction}];}
   for start in &starts{if let Ok(f)=problem.solve(*start){if best.as_ref().is_none_or(|b|f.chi_square<b.chi_square){best_ions=problem.model.ions.clone();best=Some(f);}}}
  }
  problem.model.ions=best_ions.clone();
  let result=best.context("all initializations failed").and_then(|mut f|{problem.model_uncertainty(&mut f)?;if f.status=="converged"{warm=f.plasma;}Ok(f)});
  match result{Ok(result)=>gates.push(Gate{range_m:r,ions:best_ions.clone(),result:Some(result),error:None}),Err(e)=>gates.push(Gate{range_m:r,ions:best_ions.clone(),result:None,error:Some(e.to_string())})}
 }
 ensure!(!gates.is_empty(),"no gates in selected range");drop(first_ion);drop(record);drop(records);
 let tmp=output.with_extension(format!("{}.tmp.h5",std::process::id()));
 {let out=hdf5::File::create(&tmp)?;out.new_attr_builder().with_data("ruisdap-isr-v1").create("schema")?;
 let sources=out.create_group("sources")?.create_group("combined")?;bytes(&sources,"records_json",&serde_json::to_vec(&source_meta)?)?;out.create_group("calibration")?;
 let r=out.create_group("records")?.create_group(&id)?;
 r.new_attr_builder().with_data(serde_json::json!({"id":id,"source_id":"combined","source_key":"/","density_semantics":"calibrated_ne_m3","density_units":"m^-3","station":config.station,"antenna":config.antenna,"channel":config.channel,"mode":config.mode,"radar_frequency_hz":config.radar_frequency_hz,"start_unix_s":begin,"end_unix_s":end,"kind":"ion_fit","signal_units":"K W^-1 km^-1","source_integrations":inputs.len(),"complex_convention":"R(l)=x(t)conj(x(t+l))"}).to_string().as_str()).create("metadata")?;
 bytes(&r,"configuration_json",&serde_json::to_vec(&config)?)?;array(&r,"range_m",&[range.len()],&range,"m")?;
 let ion=r.create_group("ion_line")?;array(&ion,"acf",&[range.len(),lags.len(),2],&acf,"K W^-1 km^-1")?;array(&ion,"variance",&[range.len(),lags.len()],&var,"K^2 W^-2 km^-2")?;array(&ion,"lag_s",&[lags.len()],&lags,"s")?;
 let standard=r.create_group("fit")?;
 let mut parameters=vec![f64::NAN;range.len()*4];let mut errors=parameters.clone();
 for gate in &gates{if let Some(f)=&gate.result{let i=range.iter().position(|v|*v==gate.range_m).unwrap();let p=f.plasma;parameters[i*4..i*4+4].copy_from_slice(&[p.electron_density_m3,p.ion_temperature_k,p.electron_temperature_k,p.velocity_m_s]);for j in 0..4{errors[i*4+j]=f.standard_error[j].unwrap_or(f64::NAN);}}}
 array(&standard,"parameters",&[range.len(),4],&parameters,"m^-3,K,K,m s^-1")?;array(&standard,"standard_error",&[range.len(),4],&errors,"m^-3,K,K,m s^-1")?;
 standard.new_attr_builder().with_data("[\"Ne\",\"Ti\",\"Te\",\"vi\"]").create("parameter_names")?;
 let fit=r.create_group("ion_fit")?;
 bytes(&fit,"results_json",&serde_json::to_vec(&gates)?)?;
 bytes(&fit,"model_json",&serde_json::to_vec(&model)?)?;
 bytes(&fit,"settings_json",&serde_json::to_vec(&m.settings)?)?;
 bytes(&fit,"assumptions_json",&serde_json::to_vec(&serde_json::json!({"historical_settings":m.historical_settings,"integration_s":m.integration_s,"range_averaging":"historical r^2 spatial weights; propagated diagonal measurement variance; range covariance unavailable","composition":"Holt 32/16 amu below 700 km; profiled O+/H+ grid above 700 km when historical_settings enabled","calibration":gain,"calibration_units":"W m^-3 K^-1 km^-1","transmit_power":"removed per source integration before inverse-variance averaging","response":"native finite FIR squared frequency response; locally homogeneous plasma","uncertainty":"diagonal circular-complex lag variance; calibration, composition and shared-lag covariance excluded","velocity":"R(l)=x(t)conj(x(t+l)); signed baseband spectral velocity; physical toward-away sign unverified"}))?)?;
 let p:Vec<_>=gates.iter().flat_map(|g|if let Some(f)=&g.result{let p=f.plasma;vec![p.electron_density_m3,p.ion_temperature_k,p.electron_temperature_k,p.velocity_m_s,p.ion_collision_s]}else{vec![f64::NAN;5]}).collect();
 let e:Vec<_>=gates.iter().flat_map(|g|g.result.as_ref().map(|f|f.standard_error.map(|e|e.unwrap_or(f64::NAN))).unwrap_or([f64::NAN;5])).collect();
 array(&fit,"range_m",&[gates.len()],&gates.iter().map(|g|g.range_m).collect::<Vec<_>>(),"m")?;
 array(&fit,"parameters",&[gates.len(),5],&p,"m^-3,K,K,m s^-1,s^-1")?;array(&fit,"standard_error",&[gates.len(),5],&e,"m^-3,K,K,m s^-1,s^-1")?;
 out.flush()?;}
 fs::rename(tmp,output)?;Ok(())
}

fn receiver_response(config:&Configuration)->Result<(Vec<f64>,Vec<f64>,f64)>{
 let s=&config.lpi;let mut impulse=vec![C::default();4*s.filter_half+1];impulse[2*s.filter_half]=C::new(1.,0.);let h=signal::lowpass(&impulse,s.sample_rate_hz,1.2*s.passband_hz,s.filter_half)?;
 let step=250.;let frequencies:Vec<_>=(-400..=400).map(|i|i as f64*step).collect();
 let response:Vec<f64>=frequencies.iter().map(|f|h.iter().enumerate().map(|(j,a)|*a*C::from_polar(1.,-2.*std::f64::consts::PI*f*j as f64/s.sample_rate_hz)).sum::<C>().norm_sqr()).collect();
 Ok((frequencies,response,step))
}

/// Estimate the system constant from the actual complex ACF and finite receiver
/// response at plasma-line Ne, conditional on fitted Ti, Te, velocity and ions.
/// No zero-lag extrapolation or legacy decoder conversion enters this estimate.
#[derive(Deserialize)]
pub struct CalibrationPoint {pub fit:PathBuf,pub range_m:f64,pub frequency_hz:f64}
pub fn calibration(point:&CalibrationPoint)->Result<serde_json::Value>{
 let h=hdf5::File::open(&point.fit)?;let rs=h.group("records")?;let r=rs.group(&rs.member_names()?[0])?;
 let config:Configuration=serde_json::from_slice(&r.dataset("configuration_json")?.read_raw::<u8>()?)?;
 let ranges=r.dataset("range_m")?.read_raw::<f64>()?;
 let ri=ranges.iter().position(|v|(*v-point.range_m).abs()<0.01).context("calibration gate outside native product")?;
 let gates:serde_json::Value=serde_json::from_slice(&r.dataset("ion_fit/results_json")?.read_raw::<u8>()?)?;
 let gate=gates.as_array().context("invalid gates")?.iter().find(|g|(g["range_m"].as_f64().unwrap_or(f64::NAN)-point.range_m).abs()<0.01).context("no fitted calibration gate")?;
 let result=gate.get("result").filter(|v|!v.is_null()).context("gate has no fit")?;
 ensure!(result["status"]=="converged","calibration requires a converged ion fit");
 let mut plasma:PlasmaState=serde_json::from_value(result["plasma"].clone())?;
 ensure!(point.frequency_hz.is_finite()&&(1e6..=12.5e6).contains(&point.frequency_hz.abs()),"invalid plasma frequency");
 let k=4.*std::f64::consts::PI*config.radar_frequency_hz/299792458.;
 let omega2=(std::f64::consts::TAU*point.frequency_hz).powi(2)-3.*k*k*1.380649e-23*plasma.electron_temperature_k/9.1093837139e-31;
 ensure!(omega2>0.,"thermal correction exceeds frequency squared");
 plasma.electron_density_m3=omega2*8.8541878188e-12*9.1093837139e-31/(1.602176634e-19_f64).powi(2);
 let assumptions:serde_json::Value=serde_json::from_slice(&r.dataset("ion_fit/assumptions_json")?.read_raw::<u8>()?)?;
 let gain=assumptions["calibration"].as_f64().context("missing native calibration")?;
 let mut model:SpectrumModel=serde_json::from_slice(&r.dataset("ion_fit/model_json")?.read_raw::<u8>()?)?;
 if gate["ions"].is_array(){model.ions=serde_json::from_value(gate["ions"].clone())?;}
 let settings:FitSettings=serde_json::from_slice(&r.dataset("ion_fit/settings_json")?.read_raw::<u8>()?)?;
 let (frequencies,response,step)=receiver_response(&config)?;
 let lags=r.dataset("ion_line/lag_s")?.read_raw::<f64>()?;
 let acf=r.dataset("ion_line/acf")?.read_raw::<f64>()?;let variance=r.dataset("ion_line/variance")?.read_raw::<f64>()?;
 let scale=1./(gain*(point.range_m/1000.).powi(2));
 let measurements:Vec<_>=lags.iter().enumerate().filter_map(|(li,&lag)|{
  let i=ri*lags.len()+li;let z=C::new(acf[2*i],acf[2*i+1]);let sigma=(variance[i]/2.).sqrt();
  if !z.re.is_finite()||!z.im.is_finite()||!sigma.is_finite()||sigma<=0.{return None;}
  let kernel=frequencies.iter().zip(&response).enumerate().map(|(j,(&f,&h))|C::from_polar(scale*h*step*if j==0||j+1==frequencies.len(){0.5}else{1.},-std::f64::consts::TAU*f*lag)).collect();
  Some(Measurement{lag_s:lag,value:z,sigma:[sigma,sigma],kernel})
 }).collect();
 ensure!(measurements.len()>=8,"insufficient calibration lags");
 let problem=FitProblem{model,frequency_hz:frequencies,measurements,settings};
 let predicted=problem.predict(plasma)?;
 let mut numerator=0.;let mut denominator=0.;
 for (m,z) in problem.measurements.iter().zip(predicted){let w=1./m.sigma[0].powi(2);numerator+=(m.value*z.conj()).re*w;denominator+=z.norm_sqr()*w;}
 ensure!(numerator>0.&&denominator>0.&&numerator.is_finite()&&denominator.is_finite(),"nonpositive calibrated ACF projection");
 let constant=gain*denominator/numerator;
 Ok(serde_json::json!({"magic_constant":constant,"plasma_density_m3":plasma.electron_density_m3,"multiplier":constant,"magic_constant_units":"W m^-3 K^-1 km^-1","source_calibration":gain,"range_m":point.range_m,"lag_count":problem.measurements.len()}))
}
