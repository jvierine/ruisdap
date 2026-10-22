//! Targeted residual-seed experiment on original RF; no voltage outputs.
//! This does not alter production detections, reduced products or fit settings.
use anyhow::{ensure, Context, Result};
use guisdap_millstone::{
    archive::{Archive, Manifest},
    debris::{self, Detection},
    C,
};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{collections::BTreeMap, path::Path};

fn bytes(g: &hdf5::Group, name: &str, value: &impl serde::Serialize) -> Result<()> {
    g.new_dataset_builder()
        .with_data(&serde_json::to_vec(value)?)
        .create(name)?;
    Ok(())
}
fn read(g: &hdf5::Group, name: &str) -> Result<Value> {
    Ok(serde_json::from_slice(&g.dataset(name)?.read_raw::<u8>()?)?)
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() == 4 && !Path::new(&a[3]).exists(),
        "archive_manifest.json residual_check.h5 new_report.h5"
    );
    let m: Manifest = serde_json::from_slice(&std::fs::read(&a[1])?)?;
    let source = hdf5::File::open(&a[2])?;
    ensure!(
        source.attr("complete")?.read_scalar::<bool>()?,
        "Incomplete residual check"
    );
    let mut selected: BTreeMap<u64, Vec<(u64, Vec<Detection>, Value)>> = BTreeMap::new();
    for product in source.member_names()? {
        let begin: u64 = product
            .split('-')
            .next()
            .context("No time prefix")?
            .parse()?;
        for id in source.group(&product)?.member_names()? {
            let group = source.group(&format!("{product}/{id}"))?;
            for name in group.member_names()? {
                let Ok(g) = group.group(&name) else { continue };
                let old = read(&g, "native_diagnostic_json")?;
                if old["accepted"] != true
                    || old["confirmed"] != true
                    || !g.attr("residual_layer_count")?.read_scalar::<i64>()?.gt(&0)
                {
                    continue;
                }
                let mut extra = Vec::new();
                for v in read(&g, "after_json")?
                    .as_array()
                    .context("No detections array")?
                {
                    let range = v["range_km"].as_f64().context("No range")?;
                    if !(525. ..=575.).contains(&range) {
                        continue;
                    }
                    let d = Detection {
                        delay_samples: v["corrected_delay_samples"].as_f64().unwrap(),
                        doppler_hz: v["doppler_hz"].as_f64().unwrap(),
                    };
                    // Avoid trivially duplicating the already fitted component.
                    if old["detections"].as_array().unwrap().iter().all(|o| {
                        (o["delay_samples"].as_f64().unwrap() - d.delay_samples).abs() > 2.
                            || (o["doppler_hz"].as_f64().unwrap() - d.doppler_hz).abs() > 300.
                    }) {
                        extra.push(d);
                    }
                }
                if !extra.is_empty() {
                    selected
                        .entry(begin)
                        .or_default()
                        .push((name.parse()?, extra, old));
                }
            }
        }
    }
    let report = hdf5::File::create(&a[3])?;
    report
        .new_attr::<bool>()
        .create("complete")?
        .write_scalar(&false)?;
    bytes(
        &report,
        "source_sha256_json",
        &json!({"manifest":format!("{:x}",Sha256::digest(std::fs::read(&a[1])?)),"residual_check":format!("{:x}",Sha256::digest(std::fs::read(&a[2])?))}),
    )?;
    let mut archive = Archive::open(&m)?;
    let mut count = 0;
    for (begin, rows) in selected {
        let (observations, _, _) = archive.observations(begin)?;
        for p in observations {
            for (sample, extra, old) in &rows {
                let Some(i) = p.samples.iter().position(|s| s == sample) else {
                    continue;
                };
                let dc = &p.config.instrument_metadata["dc_offset"];
                let dc = C::new(dc[0].as_f64().unwrap(), dc[1].as_f64().unwrap());
                let y: Vec<_> = p.echo[i].iter().map(|z| *z - dc).collect();
                let mut tx = p.transmit[i].clone();
                for z in &mut tx[p.config.lpi.tx_gate[0]..p.config.lpi.tx_gate[1]] {
                    *z -= dc;
                }
                let before = debris::cancel(
                    &y,
                    &tx,
                    &p.detections[i],
                    p.config.lpi.sample_rate_hz,
                    p.config.lpi.noise_gate,
                    &p.config.cancellation,
                )?;
                ensure!(
                    (before.diagnostic.residual_rms_over_noise
                        - old["residual_rms_over_noise"].as_f64().unwrap())
                    .abs()
                        < 1e-9,
                    "Original native cancellation changed for {sample}"
                );
                let mut joint = p.detections[i].clone();
                joint.extend(extra.iter().cloned());
                let g = report.create_group(&sample.to_string())?;
                bytes(&g, "configuration_json", &p.config)?;
                bytes(&g, "original_diagnostic_json", &before.diagnostic)?;
                bytes(&g, "added_seeds_json", extra)?;
                match debris::cancel(
                    &y,
                    &tx,
                    &joint,
                    p.config.lpi.sample_rate_hz,
                    p.config.lpi.noise_gate,
                    &p.config.cancellation,
                ) {
                    Ok(after) => {
                        bytes(&g, "augmented_diagnostic_json", &after.diagnostic)?;
                    }
                    Err(e) => {
                        bytes(&g, "augmentation_error_json", &e.to_string())?;
                    }
                }
                count += 1;
                report.flush()?;
            }
        }
    }
    report.attr("complete")?.write_scalar(&true)?;
    println!("{count} original-RF residual-candidate experiments completed; production unchanged");
    Ok(())
}
