//! Cancellation-only reproducibility/timing check. Optional original archive
//! input is read directly; no voltage samples or RF pulse packs are written.
use anyhow::{ensure, Result};
use guisdap_millstone::{
    archive::{Archive, Manifest},
    debris::{self, Detection, Settings},
    C,
};
use sha2::{Digest, Sha256};
use std::{path::Path, time::Instant};
fn benchmark(
    group: &hdf5::Group,
    y: &[C],
    tx: &[C],
    detections: &[Detection],
    noise: [usize; 2],
    settings: &Settings,
    repeats: usize,
) -> Result<()> {
    let mut times = Vec::new();
    for rep in 0..repeats {
        let before = Instant::now();
        let result = debris::cancel(y, tx, detections, 1e6, noise, settings)?;
        times.push(before.elapsed().as_secs_f64());
        if rep == 0 {
            let mut hash = Sha256::new();
            for z in &result.cleaned {
                hash.update(z.re.to_le_bytes());
                hash.update(z.im.to_le_bytes());
            }
            for &yes in &result.usable {
                hash.update([u8::from(yes)]);
            }
            group
                .new_attr_builder()
                .with_data(format!("{:x}", hash.finalize()).as_str())
                .create("cleaned_and_mask_sha256")?;
            group
                .new_dataset_builder()
                .with_data(&serde_json::to_vec(&result.diagnostic)?)
                .create("diagnostics_json")?;
        }
    }
    group
        .new_dataset_builder()
        .with_data(&times)
        .create("seconds")?;
    println!("{} seconds={:?}", group.name(), times);
    Ok(())
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() >= 3 && !Path::new(&a[1]).exists(),
        "new report.h5 repeats [archive_manifest.json]"
    );
    let repeats: usize = a[2].parse()?;
    let h = hdf5::File::create(&a[1])?;
    let n = 10000;
    let mut tx = vec![C::default(); n];
    for i in 76..624 {
        tx[i] = C::new(if ((i - 76) / 30) % 2 == 0 { 1. } else { -1. }, 0.);
    }
    let mut state = 71u64;
    let y: Vec<_> = (0..n)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            C::new(
                ((state >> 32) as i32 as f64) / i32::MAX as f64,
                ((state & 0xffffffff) as i32 as f64) / i32::MAX as f64,
            )
        })
        .collect();
    let ds: Vec<_> = (0..4)
        .map(|i| Detection {
            delay_samples: 1700. + i as f64 * 1200. + 0.25,
            doppler_hz: -3000. + i as f64 * 1500.,
        })
        .collect();
    let mut signal = y;
    for d in &ds {
        for (y, t) in signal.iter_mut().zip(debris::template(
            &tx,
            d.delay_samples,
            d.doppler_hz,
            1e6,
            n,
            16,
        )) {
            *y += C::new(20., 10.) * t;
        }
    }
    benchmark(
        &h.create_group("synthetic")?,
        &signal,
        &tx,
        &ds,
        [9000, 9800],
        &Settings::default(),
        repeats,
    )?;
    if let Some(path) = a.get(3) {
        let m: Manifest = serde_json::from_slice(&std::fs::read(path)?)?;
        let mut archive = Archive::open(&m)?;
        let (observations, _, _) = archive.observations(m.start_unix_s)?;
        for p in observations {
            let family = if p.config.mode.starts_with("alternating") {
                "coded"
            } else {
                "long"
            };
            let i = p
                .detections
                .iter()
                .position(|d| d.len() == 1)
                .or_else(|| p.detections.iter().position(|d| !d.is_empty()))
                .ok_or_else(|| anyhow::anyhow!("no real echo in {family}"))?;
            let dc = C::new(-0.212, -0.221);
            let y: Vec<_> = p.echo[i].iter().map(|z| *z - dc).collect();
            let mut tx = p.transmit[i].clone();
            for z in &mut tx[p.config.lpi.tx_gate[0]..p.config.lpi.tx_gate[1]] {
                *z -= dc;
            }
            let group = h.create_group(family)?;
            group
                .new_attr::<u64>()
                .create("pulse_sample")?
                .write_scalar(&p.samples[i])?;
            benchmark(
                &group,
                &y,
                &tx,
                &p.detections[i],
                p.config.lpi.noise_gate,
                &p.config.cancellation,
                repeats,
            )?;
        }
    }
    Ok(())
}
