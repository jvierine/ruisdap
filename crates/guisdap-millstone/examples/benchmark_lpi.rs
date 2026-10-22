//! Stage timings on immutable real pulse packs. Usage: output.h5 repeats input...
//! Compilation and MPI startup are excluded; CLI process timings measure those separately.
use anyhow::{ensure,Result};
use guisdap_millstone::{product,lpi};
use std::{path::{Path,PathBuf},time::Instant};
fn main()->Result<()>{
 let a:Vec<_>=std::env::args().collect();ensure!(a.len()>3,"output.h5 repeats input...");
 let output=PathBuf::from(&a[1]);ensure!(!output.exists(),"immutable benchmark output exists");
 let repeats:usize=a[2].parse()?;let h=hdf5::File::create(&output)?;
 h.new_attr_builder().with_data("read_prepare_lpi_write; seconds; first repeat retained, no cache flush").create("timing_definition")?;
 for (pi,input)in a[3..].iter().enumerate(){
  let group=h.create_group(&format!("pack_{pi}"))?;let mut times=Vec::new();let mut baseline=None;
  for rep in 0..repeats{
   let start=Instant::now();let p=product::read(Path::new(input))?;let read=start.elapsed().as_secs_f64();
   let start=Instant::now();let(pulses,diagnostics)=product::prepare(&p)?;let prepare=start.elapsed().as_secs_f64();
   let start=Instant::now();let estimate=lpi::estimate(&pulses,&p.config.lpi,p.config.injection_temperature_k)?;let lpi=start.elapsed().as_secs_f64();
   let finite=estimate.acf.iter().filter(|z|z.re.is_finite()&&z.im.is_finite()).count();ensure!(finite>100,"missing inversion");
   let bits:Vec<_>=estimate.acf.iter().flat_map(|z|[z.re.to_bits(),z.im.to_bits()]).collect();
   if let Some(old)=&baseline{ensure!(*old==bits,"repeat changed numerical ACF");}else{baseline=Some(bits);}
   let dest=output.with_file_name(format!("bench-pack-{pi}-repeat-{rep}.h5"));
   let start=Instant::now();product::write(&dest,Path::new(input),&p,&pulses,&diagnostics,Some(&estimate))?;let write=start.elapsed().as_secs_f64();
   times.extend([read,prepare,lpi,write]);
   if rep==0{
    group.new_attr_builder().with_data(input.as_str()).create("input")?;
    group.new_attr_builder().with_data(p.hash.as_str()).create("input_sha256")?;
    group.new_dataset_builder().with_data(&serde_json::to_vec(&p.config)?).create("configuration_json")?;
    group.new_dataset_builder().with_data(&[p.samples.len()as u64,p.config.lpi.science_pulse.iter().filter(|&&x|x).count()as u64,p.config.lpi.lag_samples.len()as u64,(p.config.lpi.range_end-p.config.lpi.range_start)as u64,finite as u64,p.detections.iter().map(|d|d.len()).sum::<usize>()as u64]).create("counts_stored_science_lags_ranges_finite_seeds")?;
   }
   println!("pack={pi} repeat={rep} read={read:.6} prepare={prepare:.6} lpi={lpi:.6} write={write:.6} finite={finite}");
  }
  let d=group.new_dataset::<f64>().shape([repeats,4]).create("seconds")?;d.write_raw(&times)?;
 }
 Ok(())
}
