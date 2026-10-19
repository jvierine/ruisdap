//! Cancellation-only real-echo timing; excludes disk read and catalogue search.
use anyhow::{ensure,Result};
use guisdap_millstone::product;
use std::{path::Path,time::Instant};
fn main()->Result<()>{let a:Vec<_>=std::env::args().collect();ensure!(a.len()>3,"output.h5 repeats input...");ensure!(!Path::new(&a[1]).exists(),"output exists");let h=hdf5::File::create(&a[1])?;let n:usize=a[2].parse()?;
for(i,path)in a[3..].iter().enumerate(){let p=product::read(Path::new(path))?;ensure!(p.detections.iter().any(|d|!d.is_empty()),"no CFAR seed");let g=h.create_group(&format!("case_{i}"))?;let mut times=Vec::new();for rep in 0..n{let start=Instant::now();let(_,d)=product::prepare(&p)?;let t=start.elapsed().as_secs_f64();times.push(t);if rep==0{g.new_dataset_builder().with_data(&serde_json::to_vec(&d)?).create("diagnostics_json")?;g.new_attr_builder().with_data(p.hash.as_str()).create("input_sha256")?;}println!("case={i} repeat={rep} prepare_s={t:.6}");}g.new_dataset_builder().with_data(&times).create("seconds")?;}Ok(())}
