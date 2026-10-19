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
}
#[derive(Serialize)]
struct Gate {range_m:f64,result:Option<guisdap_ion::fit::FitResult>,error:Option<String>}
fn bytes(g:&hdf5::Group,name:&str,v:&[u8])->Result<()>{g.new_dataset_builder().with_data(v).create(name)?;Ok(())}
fn array(g:&hdf5::Group,name:&str,shape:&[usize],v:&[f64],units:&str)->Result<()>{let d=g.new_dataset::<f64>().shape(shape).create(name)?;d.write_raw(v)?;d.new_attr_builder().with_data(units).create("units")?;Ok(())}
pub fn jobs(m:&Manifest)->Result<Vec<Vec<PathBuf>>>{
 let mut groups:BTreeMap<(String,i64),Vec<PathBuf>>=BTreeMap::new();
 for (i,input) in m.inputs.iter().enumerate(){let h=hdf5::File::open(input)?;let rs=h.group("records")?;let r=rs.group(&rs.member_names()?[0])?;let c:Configuration=serde_json::from_slice(&r.dataset("configuration_json")?.read_raw::<u8>()?)?;
  if !m.magic_constant_by_antenna.contains_key(&c.antenna){continue;}
  let bucket=if m.integration_s>0.{(c.start_unix_s/m.integration_s).floor()as i64}else{i as i64};groups.entry((c.antenna,bucket)).or_default().push(input.clone());
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
  let a=ion.dataset("acf")?.read_raw::<f64>()?;let v=ion.dataset("variance")?.read_raw::<f64>()?;
  for i in 0..count{let z=C::new(a[2*i],a[2*i+1])/pt;let variance=v[i]/pt.powi(2);if z.re.is_finite()&&z.im.is_finite()&&variance.is_finite()&&variance>0.{let w=1./variance;weight[i]+=w;sum[i]+=z*w;}}
  begin=begin.min(c.start_unix_s);end=end.max(c.end_unix_s);
  use sha2::{Digest,Sha256};source_meta.push(serde_json::json!({"path":input,"sha256":format!("{:x}",Sha256::digest(fs::read(input)?)),"start_unix_s":c.start_unix_s,"end_unix_s":c.end_unix_s,"transmitter_w":pt,"azimuth_deg":c.instrument_metadata["azimuth_deg"],"elevation_deg":c.instrument_metadata["elevation_deg"]}));
 }
 let acf:Vec<_>=sum.iter().zip(&weight).map(|(s,w)|if *w>0.{*s / *w}else{C::new(f64::NAN,f64::NAN)}).flat_map(|z|[z.re,z.im]).collect();
 let var:Vec<_>=weight.iter().map(|w|if *w>0.{1./w}else{f64::NAN}).collect();
 config.start_unix_s=begin;config.end_unix_s=end;
 let pt=1.; // Each measurement's recorded power was removed before averaging.
 // Model the same finite, DC-normalized receive filter used on voltages.
 // The common spectral model integrates quadrature in Hz through this kernel.
 let s=&config.lpi;let mut impulse=vec![C::default();4*s.filter_half+1];impulse[2*s.filter_half]=C::new(1.,0.);let h=signal::lowpass(&impulse,s.sample_rate_hz,1.2*s.passband_hz,s.filter_half)?;
 let step=250.;let frequencies:Vec<_>=(-400..=400).map(|i|i as f64*step).collect();
 let response:Vec<f64>=frequencies.iter().map(|f|h.iter().enumerate().map(|(j,a)|*a*C::from_polar(1.,-2.*std::f64::consts::PI*f*j as f64/s.sample_rate_hz)).sum::<C>().norm_sqr()).collect();
 let model=SpectrumModel{radar_frequency_hz:config.radar_frequency_hz,ions:m.ions.clone(),electron_collision_s:0.};model.validate()?;
 let mut gates=Vec::new();
 for (ri,&r)in range.iter().enumerate(){if r/1000.<m.range_km[0]||r/1000.>m.range_km[1]{continue;}
  let scale=pt/(gain*(r/1000.).powi(2));
  let measurements:Vec<_>=lags.iter().enumerate().filter_map(|(li,&lag)|{let i=ri*lags.len()+li;let z=C::new(acf[2*i],acf[2*i+1]);let sigma=(var[i]/2.).sqrt();if !z.re.is_finite()||!z.im.is_finite()||!sigma.is_finite()||sigma<=0.{return None;}
   let kernel=frequencies.iter().zip(&response).enumerate().map(|(j,(&f,&h))|C::from_polar(scale*h*step*if j==0||j+1==frequencies.len(){0.5}else{1.},-2.*std::f64::consts::PI*f*lag)).collect();Some(Measurement{lag_s:lag,value:z,sigma:[sigma,sigma],kernel})}).collect();
  if measurements.len()<8{gates.push(Gate{range_m:r,result:None,error:Some("fewer than eight measured lags".into())});continue;}
  let problem=FitProblem{model:model.clone(),frequency_hz:frequencies.clone(),measurements,settings:m.settings.clone()};
  let result=problem.solve(m.initial).and_then(|mut f|{problem.model_uncertainty(&mut f)?;Ok(f)});
  match result{Ok(result)=>gates.push(Gate{range_m:r,result:Some(result),error:None}),Err(e)=>gates.push(Gate{range_m:r,result:None,error:Some(e.to_string())})}
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
 bytes(&fit,"assumptions_json",&serde_json::to_vec(&serde_json::json!({"calibration":gain,"calibration_units":"W m^-3 K^-1 km^-1","transmit_power":"removed per source integration before inverse-variance averaging","response":"native finite FIR squared frequency response; locally homogeneous plasma","uncertainty":"diagonal circular-complex lag variance; calibration, composition and shared-lag covariance excluded","velocity":"R(l)=x(t)conj(x(t+l)); signed baseband spectral velocity; physical toward-away sign unverified"}))?)?;
 let p:Vec<_>=gates.iter().flat_map(|g|if let Some(f)=&g.result{let p=f.plasma;vec![p.electron_density_m3,p.ion_temperature_k,p.electron_temperature_k,p.velocity_m_s,p.ion_collision_s]}else{vec![f64::NAN;5]}).collect();
 let e:Vec<_>=gates.iter().flat_map(|g|g.result.as_ref().map(|f|f.standard_error.map(|e|e.unwrap_or(f64::NAN))).unwrap_or([f64::NAN;5])).collect();
 array(&fit,"range_m",&[gates.len()],&gates.iter().map(|g|g.range_m).collect::<Vec<_>>(),"m")?;
 array(&fit,"parameters",&[gates.len(),5],&p,"m^-3,K,K,m s^-1,s^-1")?;array(&fit,"standard_error",&[gates.len(),5],&e,"m^-3,K,K,m s^-1,s^-1")?;
 out.flush()?;}
 fs::rename(tmp,output)?;Ok(())
}
