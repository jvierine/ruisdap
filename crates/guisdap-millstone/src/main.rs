use anyhow::{ensure, Result};
use clap::{Parser, Subcommand};
use guisdap_millstone::{lpi, product};
use serde::Deserialize;
use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};
#[derive(Parser)]
#[command(about = "Native Millstone Hill voltage cancellation, spectra and weighted LPI")]
struct Args {
    #[command(subcommand)]
    command: Command,
}
#[derive(Subcommand)]
enum Command {
    Process { manifest: PathBuf },
    ProcessArchive { manifest: PathBuf },
    Fit { manifest: PathBuf },
    EstimateCalibration { manifest: PathBuf },
    Cancel { input: PathBuf, output: PathBuf },
}
#[derive(Deserialize)]
struct Manifest {
    inputs: Vec<PathBuf>,
    output: PathBuf,
    #[serde(default)]
    live_directory: Option<PathBuf>,
    #[serde(default)]
    lpi_range: Option<[usize;2]>,
    #[serde(default)]
    measurement_range: Option<[usize;2]>,
}
fn atomic_json(path: &Path, value: &serde_json::Value) -> Result<()> {
    let tmp = path.with_extension(format!("{}.tmp", std::process::id()));
    fs::write(&tmp, serde_json::to_vec(value)?)?;
    fs::rename(tmp, path)?;
    Ok(())
}
fn progress(
    dir: Option<&Path>,
    rank: usize,
    ranks: usize,
    phase: &str,
    completed: usize,
    total: usize,
    start: Instant,
) -> Result<()> {
    if let Some(dir) = dir {
        fs::create_dir_all(dir)?;
        atomic_json(
            &dir.join(format!("rank-{rank}.json")),
            &serde_json::json!({"host":std::env::var("HOSTNAME").unwrap_or_else(|_| std::fs::read_to_string("/etc/hostname").unwrap_or_else(|_|"unknown".into()).trim().into()),"rank":rank,"ranks":ranks,"phase":phase,"completed":completed,"total":total,"elapsed_s":start.elapsed().as_secs_f64(),"updated_unix_s":product::now()}),
        )?;
    }
    Ok(())
}
fn main() -> Result<()> {
    #[cfg(feature = "mpi")]
    let universe = mpi::initialize();
    #[cfg(feature = "mpi")]
    let (rank, ranks) = {
        use mpi::traits::*;
        let world = mpi::topology::SimpleCommunicator::world();
        (world.rank() as usize, world.size() as usize)
    };
    #[cfg(not(feature = "mpi"))]
    let (rank, ranks) = {
        ensure!(
            std::env::var("OMPI_COMM_WORLD_SIZE")
                .ok()
                .is_none_or(|s| s == "1")
                && std::env::var("PMI_SIZE").ok().is_none_or(|s| s == "1"),
            "MPI launcher detected; rebuild with --features mpi"
        );
        (0, 1)
    };
    let _ = rank;
    match Args::parse().command {
        Command::Cancel { input, output } => {
            ensure!(ranks == 1, "cancellation inspection takes one rank");
            let p = product::read(&input)?;
            let (pulses, diagnostics) = product::prepare(&p)?;
            product::write(&output, &input, &p, &pulses, &diagnostics, None)?;
            println!("{}", serde_json::to_string_pretty(&diagnostics)?);
        }
        Command::EstimateCalibration {manifest} => {
            ensure!(ranks==1,"calibration estimate takes one rank");
            let points:Vec<guisdap_millstone::fitting::CalibrationPoint>=serde_json::from_slice(&fs::read(&manifest)?)?;
            let estimates:Vec<_>=points.iter().map(guisdap_millstone::fitting::calibration).collect::<Result<_>>()?;
            println!("{}",serde_json::to_string(&estimates)?);
        }
        Command::Fit {manifest} => {
            let m:guisdap_millstone::fitting::Manifest=serde_json::from_slice(&fs::read(&manifest)?)?;
            fs::create_dir_all(&m.output)?;
            let mut failed=0;
            use sha2::{Digest,Sha256};
            let mut checkpoint=Sha256::new();checkpoint.update(fs::read(&manifest)?);checkpoint.update(fs::read(std::env::current_exe()?)?);let fingerprint=format!("{:x}",checkpoint.finalize());
            let jobs=guisdap_millstone::fitting::jobs(&m)?;
            ensure!(!jobs.is_empty(),"no calibrated antenna integrations");
            for (job,inputs) in jobs.iter().enumerate().filter(|(i,_)|i%ranks==rank){
                let output=m.output.join(format!("fit-{job:06}.h5"));
                let result=(||->Result<()>{
                    if output.exists(){
                        use sha2::{Digest,Sha256};
                        let h=hdf5::File::open(&output)?;ensure!(h.dataset("fit_checkpoint_sha256")?.read_raw::<u8>()?==fingerprint.as_bytes(),"fit checkpoint manifest or kernel changed");let records=h.group("records")?;let r=records.group(&records.member_names()?[0])?;
                        r.dataset("ion_fit/results_json")?;
                        let sources:Vec<serde_json::Value>=serde_json::from_slice(&h.dataset("sources/combined/records_json")?.read_raw::<u8>()?)?;
                        let old:std::collections::BTreeSet<_>=sources.iter().map(|s|s["path"].as_str().unwrap_or("").to_owned()).collect();
                        let new:std::collections::BTreeSet<_>=inputs.iter().map(|p|p.to_string_lossy().into_owned()).collect();ensure!(old==new,"fit checkpoint input set changed");
                        for s in sources{ensure!(format!("{:x}",Sha256::digest(fs::read(s["path"].as_str().unwrap())?))==s["sha256"].as_str().unwrap_or(""),"fit checkpoint input changed");}
                        return Ok(());
                    }
                    guisdap_millstone::fitting::run(inputs,&output,&m)?;
                    let h=hdf5::File::open_rw(&output)?;h.new_dataset_builder().with_data(fingerprint.as_bytes()).create("fit_checkpoint_sha256")?;h.flush()?;Ok(())
                })();
                if let Err(e)=result{failed+=1;eprintln!("rank {rank} fit {job}: {e:#}");}
            }
            ensure!(failed==0,"rank {rank}: {failed} fit jobs failed");
        }
        Command::ProcessArchive {manifest} => {
            let m:guisdap_millstone::archive::Manifest=serde_json::from_slice(&fs::read(&manifest)?)?;
            ensure!(m.end_unix_s>m.start_unix_s && (m.end_unix_s-m.start_unix_s)%10==0,"invalid archive integration window");
            let mut archive=guisdap_millstone::archive::Archive::open(&m)?;
            let total=((m.end_unix_s-m.start_unix_s)/10)as usize;
            // Allocate small contiguous blocks dynamically. Workers keep the RF
            // cache warm within each block; heavy CFAR passes cannot pin one rank.
            fs::create_dir_all(&m.output)?;let counter=m.output.join(format!("work-counter-{}-{}",m.start_unix_s,m.end_unix_s));
            if rank==0{fs::write(&counter,0u64.to_le_bytes())?;}
            #[cfg(feature="mpi")]{use mpi::traits::*;mpi::topology::SimpleCommunicator::world().barrier();}
            use std::io::{Read as _,Write as _,Seek as _};
            let mut queue=std::fs::OpenOptions::new().read(true).write(true).open(counter)?;
            let start=Instant::now();let mut failed=0;let mut done=0;
            loop{
                queue.lock()?;queue.rewind()?;let mut bytes=[0u8;8];queue.read_exact(&mut bytes)?;let lo=u64::from_le_bytes(bytes)as usize;
                if lo>=total{queue.unlock()?;break;}let hi=(lo+1).min(total);queue.rewind()?;queue.write_all(&(hi as u64).to_le_bytes())?;queue.unlock()?;
                for job in lo..hi{
                    let begin=m.start_unix_s+10*job as u64;
                    progress(m.live_directory.as_deref(),rank,ranks,"reading DigitalRF and reducing",done,total.div_ceil(ranks),start)?;
                    if let Err(e)=guisdap_millstone::archive::run_window(&mut archive,&m,begin){failed+=1;eprintln!("rank {rank} archive window {begin}: {e:#}");}done+=1;
                }
            }
            progress(m.live_directory.as_deref(),rank,ranks,if failed==0{"complete"}else{"failed"},done,done,start)?;
            ensure!(failed==0,"rank {rank}: {failed} archive windows failed");
        }
        Command::Process { manifest } => {
            let m: Manifest = serde_json::from_slice(&fs::read(&manifest)?)?;
            ensure!(!m.inputs.is_empty(), "empty manifest");
            fs::create_dir_all(&m.output)?;
            let jobs: Vec<_> = m
                .inputs
                .iter()
                .enumerate()
                .filter(|(i, _)| i % ranks == rank)
                .collect();
            let start = Instant::now();
            let mut failed = 0;
            for (done, (job, input)) in jobs.iter().enumerate() {
                progress(
                    m.live_directory.as_deref(),
                    rank,
                    ranks,
                    "reading",
                    done,
                    jobs.len(),
                    start,
                )?;
                let result = (|| -> Result<()> {
                    let mut p = product::read(input)?;
                    if let Some([lo,hi])=m.lpi_range{p.config.lpi.range_start=lo;p.config.lpi.range_end=hi;}
                    if let Some([lo,hi])=m.measurement_range{p.config.lpi.measurement_start=lo;p.config.lpi.measurement_end=hi;}
                    ensure!(
                        p.samples
                            .iter()
                            .collect::<std::collections::BTreeSet<_>>()
                            .len()
                            == p.samples.len(),
                        "duplicate pulse sample keys cannot enter LPI"
                    );
                    ensure!(
                        p.pair.iter().enumerate().all(|(i, j)| i != *j),
                        "LPI requires distinct ground-clutter partners"
                    );
                    let (pulses, diagnostics) = product::prepare(&p)?;
                    progress(
                        m.live_directory.as_deref(),
                        rank,
                        ranks,
                        "inverting_lags",
                        done,
                        jobs.len(),
                        start,
                    )?;
                    let estimate = lpi::estimate_with_progress(
                        &pulses,
                        &p.config.lpi,
                        p.config.injection_temperature_k,
                        |lag, total_lags, estimate| {
                            progress(
                                m.live_directory.as_deref(),
                                rank,
                                ranks,
                                &format!("lag {lag}/{total_lags}"),
                                done,
                                jobs.len(),
                                start,
                            )?;
                            if let Some(dir) = &m.live_directory {
                                let lease = fs::read(dir.join("lease.json")).ok().and_then(|b| {
                                    serde_json::from_slice::<serde_json::Value>(&b).ok()
                                });
                                if lease.is_some_and(|v| {
                                    v["expires_unix_s"]
                                        .as_f64()
                                        .is_some_and(|t| t > product::now())
                                }) {
                                    let nr = p.config.lpi.range_end - p.config.lpi.range_start;
                                    let nl = p.config.lpi.lag_samples.len();
                                    let preview: Vec<_> = (0..nr)
                                        .map(|r| {
                                            let z = estimate.acf[r * nl];
                                            if z.re.is_finite() {
                                                Some(z.re)
                                            } else {
                                                None
                                            }
                                        })
                                        .collect();
                                    atomic_json(
                                        &dir.join(format!("preview-{rank}.json")),
                                        &serde_json::json!({"rank":rank,"updated_unix_s":product::now(),"station":p.config.station,"antenna":p.config.antenna,"mode":p.config.mode,"first_lag_real":preview,"retained_fraction":estimate.retained_fraction,"cancelled":diagnostics.iter().filter(|d|d.accepted&&d.applied).count(),"failed_removals":diagnostics.iter().filter(|d|!d.accepted).count()}),
                                    )?;
                                }
                            }
                            Ok(())
                        },
                    )?;
                    let out = m.output.join(format!("{job:06}-{}.h5", &p.hash[..16]));
                    product::write(&out, input, &p, &pulses, &diagnostics, Some(&estimate))?;
                    Ok(())
                })();
                if let Err(e) = result {
                    failed += 1;
                    eprintln!("rank {rank} job {job} failed: {e:#}");
                }
                progress(
                    m.live_directory.as_deref(),
                    rank,
                    ranks,
                    "processing",
                    done + 1,
                    jobs.len(),
                    start,
                )?;
            }
            progress(
                m.live_directory.as_deref(),
                rank,
                ranks,
                if failed == 0 { "complete" } else { "failed" },
                jobs.len(),
                jobs.len(),
                start,
            )?;
            ensure!(failed == 0, "rank {rank}: {failed} jobs failed");
        }
    }
    #[cfg(feature = "mpi")]
    drop(universe);
    Ok(())
}
