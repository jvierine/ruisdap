//! Export the measured fourteen-slice power as HDF5 plus bounded browser tiles.
use anyhow::{ensure, Context, Result};
use guisdap_core::{self as core, FrameMeta};
use std::{fs, path::{Path, PathBuf}};
const NR: usize = 431;
const TILE: usize = 256;
fn binary(path: &Path, v: &[f32]) -> Result<()> {
    let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
    fs::write(path, bytes)?; Ok(())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().collect();
    ensure!(args.len()==3,"usage: export_beata_power INPUT_HDF5_DIRECTORY NEW_OUTPUT_DIRECTORY");
    let input=PathBuf::from(&args[1]); let out=PathBuf::from(&args[2]);
    ensure!(!out.exists(),"choose a new output version");
    let mut frames=Vec::new();
    for entry in fs::read_dir(input)? {
        let p=entry?.path(); if p.extension().is_none_or(|e|e!="h5"){continue;}
        let h=hdf5::File::open(&p)?;
        let m:FrameMeta=serde_json::from_slice(&h.dataset("metadata_json")?.read_raw::<u8>()?)?;
        if h.link_exists("ion_short_slices"){frames.push((m,p));}
    }
    frames.sort_by(|a,b|a.0.start_unix.total_cmp(&b.0.start_unix));
    ensure!(!frames.is_empty(),"no fourteen-slice records");
    let mut starts=Vec::new();let mut ends=Vec::new();let mut values=Vec::new();
    let mut sources=Vec::new();let mut altitudes=Vec::new();let mut native=Vec::new();let mut native_altitudes=Vec::new();
    for (m,p) in frames {
        ensure!(ends.last().is_none_or(|t:&f64|*t<=m.start_unix+0.01),"overlapping recordings");
        let h=hdf5::File::open(&p)?;
        let slices=h.dataset("ion_short_slices")?.read_raw::<f32>()?;
        let sum=h.dataset("ion_short_acf")?.read_raw::<f32>()?;
        for i in 0..NR {
            let v=(0..14).map(|j|slices[(j*NR+i)*2] as f64).sum::<f64>();
            ensure!((v-sum[2*i] as f64).abs() <= 2e-6*(0..14).map(|j|slices[(j*NR+i)*2].abs() as f64).sum::<f64>().max(1.),"slice sum mismatch: {} gate {i}",m.id);
        }
        let gain=h.dataset("ion_calibration_reference")?.read_raw::<f64>().or_else(|_|h.dataset("ion_calibration_reference")?.read_raw::<f32>().map(|v|v.into_iter().map(f64::from).collect()))?[0]
            -h.dataset("noise_reference")?.read_raw::<f64>().or_else(|_|h.dataset("noise_reference")?.read_raw::<f32>().map(|v|v.into_iter().map(f64::from).collect()))?[0];
        let par=h.dataset("d_parbl")?.read_raw::<f64>()?;
        values.extend(core::ion_short_slice_raw_ne(&slices,gain,par[20],m.transmitter_w).with_context(||m.id.clone())?);
        let dt=(m.end_unix-m.start_unix)/14.;
        for j in 0..14 {starts.push((m.start_unix+j as f64*dt).max(ends.last().copied().unwrap_or(f64::NEG_INFINITY)));ends.push(m.start_unix+(j+1)as f64*dt);}
        let a:Vec<_>=(0..NR).map(|i|core::altitude_km(core::ion_short_range_km(i),m.elevation_deg)).collect();
        // Resample each native profile to a fixed altitude grid using the
        // nearest measured gate; no interpolation of power or time samples.
        if altitudes.is_empty(){altitudes=a.clone();}
        let base=values.len()-14*NR;native.extend_from_slice(&values[base..]);native_altitudes.extend_from_slice(&a);
        let source=values[base..].to_vec();
        for (i,&h) in altitudes.iter().enumerate(){
            let k=a.partition_point(|v|*v<h).min(NR-1);
            let k=if k>0&&(a[k-1]-h).abs()<(a[k]-h).abs(){k-1}else{k};
            for j in 0..14{values[base+j*NR+i]=if (a[k]-h).abs()<1.0{source[j*NR+k]}else{f32::NAN};}
        }
        sources.push(m);
    }
    fs::create_dir_all(&out)?;
    let h=hdf5::File::create(out.join("power.h5"))?;
    h.new_dataset_builder().with_data(&starts).create("start_unix_s")?;
    h.new_dataset_builder().with_data(&ends).create("end_unix_s")?;
    h.new_dataset_builder().with_data(&altitudes).create("altitude_km")?;
    h.new_dataset::<f32>().shape([starts.len(),NR]).chunk([256,NR]).deflate(4).create("raw_ne_m3")?.write_raw(&values)?;
    h.new_dataset_builder().with_data(serde_json::to_vec(&sources)?.as_slice()).create("source_metadata_json")?;
    h.new_dataset::<f32>().shape([starts.len(),NR]).chunk([256,NR]).deflate(4).create("native_raw_ne_m3")?.write_raw(&native)?;
    h.new_dataset::<f64>().shape([sources.len(),NR]).create("native_altitude_km")?.write_raw(&native_altitudes)?;
    let definition=serde_json::json!({"units":"m^-3","lag_s":1e-5,"Te_Ti":1,"base_magic_constant":0.7,"timing":"Fourteen equal-duration subintervals within recorded start/end bounds; pulse-level timing not retained in these products. Boundary overlaps below 10 ms are clipped. Native pointing grids and densities are retained separately; display uses nearest altitude gate, maximum displacement 1 km.","aggregation":"Arithmetic temporal mean of calibrated signed lag power, before logarithmic colour normalization; never a mean of log values."});
    h.new_dataset_builder().with_data(serde_json::to_vec(&definition)?.as_slice()).create("definition_json")?;
    let mut levels=Vec::new();
    for factor in [1usize,4,16,64,256] {
        let dir=out.join(format!("level-{factor}"));fs::create_dir(&dir)?;
        let mut bins:Vec<(usize,usize)>=Vec::new();let mut first=0;
        for i in 1..=starts.len() {
            if i==starts.len() || i-first==factor || starts[i]-ends[i-1]>0.001 {
                bins.push((first,i));first=i;
            }
        }
        let mut tiles=Vec::new();
        for (k,chunk) in bins.chunks(TILE).enumerate() {
            let mut edges=Vec::new();let mut cols:Vec<Option<(usize,usize)>>=Vec::new();
            for &(a,b) in chunk {
                if let Some(&prev)=edges.last(){if starts[a]-prev>0.001 {cols.push(None);edges.push(starts[a]);}}
                else {edges.push(starts[a]);}
                cols.push(Some((a,b)));edges.push(ends[b-1]);
            }
            let mut v=vec![f32::NAN;cols.len()*NR];
            for (j,col) in cols.iter().enumerate(){if let Some((a,b))=*col{for r in 0..NR{let mut n=0;let mut sum=0.;for t in a..b{let x=values[t*NR+r];if x.is_finite(){sum+=x as f64;n+=1;}}if n>0{v[r*cols.len()+j]=(sum/n as f64)as f32;}}}}
            let name=format!("level-{factor}/{k:05}.f32");binary(&out.join(&name),&v)?;
            tiles.push(serde_json::json!({"url":name,"x_edges":edges,"columns":cols.len(),"start":edges[0],"end":edges[edges.len()-1]}));
        }
        levels.push(serde_json::json!({"factor":factor,"cadence_s":factor as f64*5./14.,"tiles":tiles}));
    }
    fs::write(out.join("index.json"),serde_json::to_vec(&serde_json::json!({"schema":"beata-power-tiles-v1","profiles":starts.len(),"range_gates":NR,"start":starts[0],"end":ends[ends.len()-1],"altitude_km":altitudes,"definition":definition,"levels":levels}))?)?;
    println!("{} native profiles; HDF5 and five tiled resolutions published to {}",starts.len(),out.display());Ok(())
}
