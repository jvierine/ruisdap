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
    Cancel { input: PathBuf, output: PathBuf },
}
#[derive(Deserialize)]
struct Manifest {
    inputs: Vec<PathBuf>,
    output: PathBuf,
    #[serde(default)]
    live_directory: Option<PathBuf>,
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
        Command::Process { manifest } => {
            let m: Manifest = serde_json::from_slice(&fs::read(manifest)?)?;
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
                    let p = product::read(input)?;
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
