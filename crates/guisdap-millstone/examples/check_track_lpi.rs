//! Original-RF, in-memory target injection through the complete native LPI.
//! No RF is written. Predictions use other pulses only; this diagnostic does
//! not change production subtraction, source checkpoints or published fits.
use anyhow::{ensure, Result};
use guisdap_millstone::{
    archive::{Archive, Manifest},
    debris, lpi,
    product::{self, Pack},
    signal::percentile,
    track::{self, Event, Forecaster, Settings},
    C,
};
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{path::Path, time::Instant};

fn bytes(g: &hdf5::Group, name: &str, value: &impl serde::Serialize) -> Result<()> {
    g.new_dataset_builder()
        .with_data(&serde_json::to_vec(value)?)
        .create(name)?;
    Ok(())
}
fn estimate(g: &hdf5::Group, name: &str, e: &lpi::Estimate, nr: usize, nl: usize) -> Result<()> {
    let h = g.create_group(name)?;
    for (name, values) in [
        ("acf", &e.acf),
        ("ground_clutter_acf", &e.ground_clutter_acf),
    ] {
        let data: Vec<_> = values.iter().flat_map(|z| [z.re, z.im]).collect();
        let dataset = h.new_dataset::<f64>().shape([nr, nl, 2]).create(name)?;
        dataset.write_raw(&data)?;
        dataset
            .new_attr_builder()
            .with_data("K km^-1")
            .create("units")?;
        dataset
            .new_attr_builder()
            .with_data("range,lag,complex_component")
            .create("dimensions")?;
    }
    h.new_dataset::<f64>()
        .shape([nr, nl])
        .create("variance")?
        .write_raw(&e.variance)?;
    h.new_dataset_builder()
        .with_data(&e.retained_fraction)
        .create("retained_fraction")?;
    for (name, value) in [
        ("noise_power", e.noise_power),
        ("injection_power", e.injection_power),
        ("receiver_gain", e.receiver_gain),
        ("system_temperature_k", e.system_temperature_k),
    ] {
        h.new_attr::<f64>().create(name)?.write_scalar(&value)?;
    }
    Ok(())
}
fn clone_pulses(v: &[lpi::Pulse]) -> Vec<lpi::Pulse> {
    v.iter()
        .map(|p| lpi::Pulse {
            echo: p.echo.clone(),
            transmit: p.transmit.clone(),
            pair: p.pair,
            usable: p.usable.clone(),
        })
        .collect()
}
fn main() -> Result<()> {
    let a: Vec<_> = std::env::args().collect();
    ensure!(
        a.len() == 5 && !Path::new(&a[4]).exists(),
        "archive.json begin_unix_s injected_slant_range_km new_report.h5"
    );
    let m: Manifest = serde_json::from_slice(&std::fs::read(&a[1])?)?;
    let begin: u64 = a[2].parse()?;
    let range: f64 = a[3].parse()?;
    ensure!(
        (250. ..=575.).contains(&range),
        "Control injection must avoid receiver calibration gates"
    );
    let mut archive = Archive::open(&m)?;
    let (observations, gaps, _) = archive.observations(begin)?;
    ensure!(!observations.is_empty(), "No source observations");
    let report = hdf5::File::create(&a[4])?;
    report
        .new_attr::<bool>()
        .create("complete")?
        .write_scalar(&false)?;
    report.new_attr_builder().with_data("Original RF plus declared synthetic moving coherent target; candidate-controlled coarse seeds, not a CFAR false-alarm test. No raw or cleaned voltages stored. Production unchanged.").create("scope")?;
    bytes(
        &report,
        "input_manifest_json",
        &serde_json::from_slice::<serde_json::Value>(&std::fs::read(&a[1])?)?,
    )?;
    bytes(&report, "source_gaps_json", &gaps)?;
    bytes(
        &report,
        "kernel_sha256_json",
        &format!(
            "{:x}",
            Sha256::digest(std::fs::read(std::env::current_exe()?)?)
        ),
    )?;
    let start = Instant::now();
    let mut forecaster = Forecaster::default();
    for p in observations {
        ensure!(
            !p.config.lpi.code_energy_samples.is_empty(),
            "Control requires explicit physical pulse/range normalization"
        );
        ensure!(
            p.detections.iter().all(Vec::is_empty),
            "Select an existing no-catalogue-echo source window for this control"
        );
        let family = if p.config.mode.starts_with("alternating") {
            "coded"
        } else {
            "long"
        };
        let g = report.create_group(family)?;
        bytes(&g, "source_configuration_json", &p.config)?;
        bytes(&g, "source_sha256_json", &p.hash)?;
        g.new_dataset_builder()
            .with_data(&p.samples)
            .create("pulse_sample")?;
        let (original, original_qa) = product::prepare(&p)?;
        let noise = percentile(
            &original_qa
                .iter()
                .map(|q| q.noise_rms * q.noise_rms)
                .collect::<Vec<_>>(),
            0.5,
        );
        let energy = percentile(
            &original
                .iter()
                .map(|q| q.transmit.iter().map(|z| z.norm_sqr()).sum::<f64>())
                .collect::<Vec<_>>(),
            0.5,
        );
        let amplitude = (10f64.powf(3.7) * noise / energy).sqrt();
        ensure!(
            amplitude.is_finite() && amplitude > 0.,
            "No finite injection scaling"
        );
        bytes(
            &g,
            "injection_json",
            &json!({"slant_range_km":range,"relative_amplitude_at_start":amplitude,"nominal_matched_snr_db":37.,"doppler_at_start_hz":1000.,"doppler_slope_hz_s":-250.,"amplitude_slope_per_s":0.02,"phase_at_start_rad":0.3,"coarse_seed_delay_quantum_samples":8.,"coarse_seed_doppler_quantum_hz":250.,"begin_unix_s":begin}),
        )?;
        let mut injected = Pack {
            config: p.config.clone(),
            echo: p.echo.clone(),
            transmit: p.transmit.clone(),
            pair: p.pair.clone(),
            detections: vec![vec![]; p.echo.len()],
            samples: p.samples.clone(),
            hash: String::new(),
        };
        let mut uncorrected = clone_pulses(&original);
        let mut input_hash = Sha256::new();
        for i in 0..p.samples.len() {
            let t = (p.samples[i] as i128 - begin as i128 * 1_000_000) as f64 / 1e6;
            let f = 1000. - 250. * t;
            let d = range * 2e3 / 299792458. * 1e6 - 1e6 / 440.2e6 * (1000. * t - 125. * t * t);
            let a = C::from_polar(
                amplitude * (1. + 0.02 * t),
                std::f64::consts::TAU * (1000. * t - 125. * t * t) + 0.3,
            );
            let target = debris::template(
                &original[i].transmit,
                d,
                f,
                1e6,
                original[i].echo.len(),
                p.config.cancellation.interpolation_half,
            );
            for (j, z) in target.iter().enumerate() {
                injected.echo[i][j] += a * *z;
                uncorrected[i].echo[j] += a * *z;
            }
            injected.detections[i].push(debris::Detection {
                delay_samples: (d / 8.).round() * 8.,
                doppler_hz: (f / 250.).round() * 250.,
            });
            for z in &injected.echo[i] {
                input_hash.update(z.re.to_le_bytes());
                input_hash.update(z.im.to_le_bytes());
            }
        }
        injected.hash = format!("{:x}", input_hash.finalize());
        bytes(&g, "injected_voltage_sha256_json", &injected.hash)?;
        let (legacy, legacy_qa) = product::prepare(&injected)?;
        let mut events = Vec::new();
        for (i, q) in legacy_qa.iter().enumerate() {
            if !p.config.lpi.science_pulse[i]
                || [
                    "clipped_adc",
                    "invalid_background_noise",
                    "background_gate_contaminated",
                ]
                .contains(&q.reason.as_str())
            {
                continue;
            }
            for ((d, a), snr) in q
                .detections
                .iter()
                .zip(&q.amplitudes)
                .zip(&q.amplitude_snr_db)
            {
                events.push(Event {
                    sample: p.samples[i],
                    delay_samples: d.delay_samples,
                    doppler_hz: d.doppler_hz,
                    amplitude: C::new(a[0], a[1]),
                    snr_db: *snr,
                });
            }
        }
        let settings = Settings {
            radius_s: if family == "long" { 0.5 } else { 0.35 },
            ..Settings::default()
        };
        let mut independent = clone_pulses(&uncorrected);
        let mut forecasts = Vec::new();
        let mut prediction_failures = Vec::new();
        let mut predicted_count = 0;
        for i in 0..p.samples.len() {
            let mut models = Vec::new();
            let mut failures = Vec::new();
            for group in track::groups(&events, p.samples[i], &settings)? {
                match forecaster.forecast(
                    &group,
                    p.samples[i],
                    (p.config.lpi.tx_gate[0] + p.config.lpi.tx_gate[1]) as f64 / 2.,
                    &settings,
                ) {
                    Ok(model) => {
                        if model.explained_fraction < 0.9 {
                            failures.push(format!(
                                "training_explained_fraction={}",
                                model.explained_fraction
                            ));
                            continue;
                        }
                        ensure!(
                            !model.training_samples.contains(&p.samples[i]),
                            "Held pulse leaked into training"
                        );
                        let template = debris::template(
                            &original[i].transmit,
                            model.delay_samples,
                            model.doppler_hz,
                            1e6,
                            original[i].echo.len(),
                            p.config.cancellation.interpolation_half,
                        );
                        for (j, z) in template.iter().enumerate() {
                            independent[i].echo[j] -= model.amplitude * *z;
                        }
                        models.push(model);
                    }
                    Err(error) => failures.push(error.to_string()),
                }
            }
            predicted_count += usize::from(!models.is_empty());
            forecasts.push(models);
            prediction_failures.push(failures);
        }
        bytes(&g, "independent_forecasts_json", &forecasts)?;
        bytes(&g, "legacy_cancellation_diagnostics_json", &legacy_qa)?;
        bytes(&g, "prediction_failures_json", &prediction_failures)?;
        g.new_attr::<u64>()
            .create("predicted_pulses")?
            .write_scalar(&(predicted_count as u64))?;
        let nr = p.config.lpi.range_end - p.config.lpi.range_start;
        let nl = p.config.lpi.lag_samples.len();
        g.new_dataset_builder()
            .with_data(
                &p.config
                    .lpi
                    .lag_samples
                    .iter()
                    .map(|l| *l as u64)
                    .collect::<Vec<_>>(),
            )
            .create("lag_samples")?;
        g.new_dataset_builder()
            .with_data(
                &(p.config.lpi.range_start..p.config.lpi.range_end)
                    .map(|i| i as f64 * p.config.lpi.decimation as f64 / 1e6 * 299792.458 / 2.)
                    .collect::<Vec<_>>(),
            )
            .create("slant_range_km")?;
        for (name, pulses) in [
            ("original", &original),
            ("injected", &uncorrected),
            ("same_pulse", &legacy),
            ("independent", &independent),
        ] {
            let e = lpi::estimate(pulses, &p.config.lpi, p.config.injection_temperature_k)?;
            estimate(&g, name, &e, nr, nl)?;
            println!(
                "{family} {name} gain={} elapsed={:.3}s",
                e.receiver_gain,
                start.elapsed().as_secs_f64()
            );
        }
        report.flush()?;
    }
    report.attr("complete")?.write_scalar(&true)?;
    println!("Native LPI control complete; production unchanged");
    Ok(())
}
