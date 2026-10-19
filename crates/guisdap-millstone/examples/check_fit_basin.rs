//! Independent start sensitivity check on a stored native fit's measured ACF.
use anyhow::{Context,Result};
use guisdap_ion::{fit::{FitProblem,FitSettings,Measurement},spectrum::{SpectrumModel,PlasmaState}};
use guisdap_millstone::{signal,C};
fn json<T:serde::de::DeserializeOwned>(g:&hdf5::Group,n:&str)->Result<T>{Ok(serde_json::from_slice(&g.dataset(n)?.read_raw::<u8>()?)?)}
fn main()->Result<()>{
 let args:Vec<_>=std::env::args().collect();let h=hdf5::File::open(&args[1])?;let rs=h.group("records")?;let r=rs.group(&rs.member_names()?[0])?;
 let rg:f64=args[2].parse()?;let range=r.dataset("range_m")?.read_raw::<f64>()?;let i=range.iter().position(|r|(*r-rg*1000.).abs()<2500.).context("gate missing")?;
 let config:guisdap_millstone::product::Configuration=json(&r,"configuration_json")?;
 let gates:serde_json::Value=json(&r,"ion_fit/results_json")?;let gate=gates.as_array().unwrap().iter().find(|g|(g["range_m"].as_f64().unwrap()-range[i]).abs()<0.01).unwrap();
 let mut model:SpectrumModel=json(&r,"ion_fit/model_json")?;model.ions=serde_json::from_value(gate["ions"].clone())?;
 let settings:FitSettings=json(&r,"ion_fit/settings_json")?;let assumptions:serde_json::Value=json(&r,"ion_fit/assumptions_json")?;
 let gain=assumptions["calibration"].as_f64().unwrap();let ls=r.dataset("ion_line/lag_s")?.read_raw::<f64>()?;let acf=r.dataset("ion_line/acf")?.read_raw::<f64>()?;let var=r.dataset("ion_line/variance")?.read_raw::<f64>()?;
 let mut impulse=vec![C::default();4*config.lpi.filter_half+1];impulse[2*config.lpi.filter_half]=C::new(1.,0.);let filter=signal::lowpass(&impulse,config.lpi.sample_rate_hz,1.2*config.lpi.passband_hz,config.lpi.filter_half)?;
 let frequencies:Vec<_>=(-400..=400).map(|j|j as f64*250.).collect();let response:Vec<_>=frequencies.iter().map(|f|filter.iter().enumerate().map(|(j,a)|*a*C::from_polar(1.,-std::f64::consts::TAU*f*j as f64/config.lpi.sample_rate_hz)).sum::<C>().norm_sqr()).collect();
 let scale=1./(gain*(range[i]/1000.).powi(2));let measurements=ls.iter().enumerate().filter_map(|(li,&lag)|{let j=i*ls.len()+li;let value=C::new(acf[2*j],acf[2*j+1]);let sigma=(var[j]/2.).sqrt();if !value.re.is_finite()||!value.im.is_finite()||!sigma.is_finite()||sigma<=0.{return None;}
 let kernel=frequencies.iter().zip(&response).enumerate().map(|(j,(&f,&h))|C::from_polar(scale*h*250.*if j==0||j==800{0.5}else{1.},-std::f64::consts::TAU*f*lag)).collect();Some(Measurement{lag_s:lag,value,sigma:[sigma,sigma],kernel})}).collect();
 let problem=FitProblem{model,frequency_hz:frequencies,measurements,settings};println!("stored {}",gate["result"]);
 for ne in [1e11,5e11,1e12]{for (ti,te)in [(500.,550.),(1000.,1100.),(1200.,1300.),(1000.,2000.)]{let initial=PlasmaState{electron_density_m3:ne,ion_temperature_k:ti,electron_temperature_k:te,velocity_m_s:0.,ion_collision_s:0.};let f=problem.solve(initial)?;println!("{}",serde_json::json!({"start":initial,"result":f.plasma,"chi":f.reduced_chi_square,"status":f.status,"iterations":f.iterations}));}}
 Ok(())
}
