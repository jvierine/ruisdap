//! Native HDF5 pulse-pack boundary and common ISR product output.
use crate::{
    debris::{self, Detection, Diagnostic},
    lpi, signal, C,
};
use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::Path,
    time::{SystemTime, UNIX_EPOCH},
};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Configuration {
    #[serde(default)]
    pub instrument_metadata: serde_json::Value,
    pub station: String,
    pub antenna: String,
    pub channel: String,
    pub mode: String,
    pub radar_frequency_hz: f64,
    pub start_unix_s: f64,
    pub end_unix_s: f64,
    pub lpi: lpi::Settings,
    #[serde(default)]
    pub cancellation: debris::Settings,
    #[serde(default)]
    pub cancel_debris: bool,
    pub injection_temperature_k: f64,
    pub range_doppler_fft: usize,
    pub range_doppler_samples: Vec<usize>,
}
pub struct Pack {
    pub config: Configuration,
    pub echo: Vec<Vec<C>>,
    pub transmit: Vec<Vec<C>>,
    pub pair: Vec<usize>,
    pub detections: Vec<Vec<Detection>>,
    pub samples: Vec<u64>,
    pub hash: String,
}
fn complex(h: &hdf5::File, name: &str) -> Result<Vec<Vec<C>>> {
    let d = h.dataset(name)?;
    let shape = d.shape();
    ensure!(
        shape.len() == 3 && shape[2] == 2,
        "{name}: expected pulse,sample,complex_component"
    );
    let v = d.read_raw::<f64>()?;
    Ok(v.chunks_exact(shape[1] * 2)
        .map(|p| p.chunks_exact(2).map(|z| C::new(z[0], z[1])).collect())
        .collect())
}
pub fn read(path: &Path) -> Result<Pack> {
    let h = hdf5::File::open(path)?;
    let schema = h.dataset("schema_utf8")?.read_raw::<u8>()?;
    ensure!(schema == b"ruisdap-pulses-v1", "unsupported pulse pack");
    let config: Configuration =
        serde_json::from_slice(&h.dataset("configuration_json")?.read_raw::<u8>()?)?;
    ensure!(
        config.station == "Millstone Hill"
            && config.radar_frequency_hz > 0.
            && config.start_unix_s.is_finite()
            && config.end_unix_s > config.start_unix_s,
        "invalid observation identity"
    );
    let echo = complex(&h, "echo")?;
    let transmit = complex(&h, "transmit")?;
    let pair: Vec<_> = h
        .dataset("pair_index")?
        .read_raw::<u64>()?
        .iter()
        .map(|v| *v as usize)
        .collect();
    let samples = h.dataset("pulse_sample")?.read_raw::<u64>()?;
    ensure!(
        echo.len() == transmit.len() && echo.len() == pair.len() && echo.len() == samples.len(),
        "pulse metadata length mismatch"
    );
    ensure!(
        pair.iter().enumerate().all(|(_, &p)| p < pair.len()),
        "invalid/self ground-clutter pair"
    );
    let mut detections = vec![Vec::new(); echo.len()];
    if h.link_exists("detections") {
        let d = h.dataset("detections")?;
        ensure!(
            d.shape().len() == 2 && d.shape()[1] == 3,
            "detections: pulse index, delay samples, Doppler Hz required"
        );
        for row in d.read_raw::<f64>()?.chunks_exact(3) {
            ensure!(
                row[0].is_finite()
                    && row[0] >= 0.
                    && row[0].fract() == 0.
                    && (row[0] as usize) < echo.len(),
                "invalid detection pulse index"
            );
            detections[row[0] as usize].push(Detection {
                delay_samples: row[1],
                doppler_hz: row[2],
            });
        }
    }
    let hash = format!("{:x}", Sha256::digest(fs::read(path)?));
    Ok(Pack {
        config,
        echo,
        transmit,
        pair,
        detections,
        samples,
        hash,
    })
}
fn array(
    h: &hdf5::Group,
    name: &str,
    shape: &[usize],
    data: &[f64],
    units: &str,
    dimensions: &str,
) -> Result<()> {
    let d = h.new_dataset::<f64>().shape(shape).create(name)?;
    d.write_raw(data)?;
    d.new_attr_builder().with_data(units).create("units")?;
    d.new_attr_builder()
        .with_data(dimensions)
        .create("dimensions")?;
    Ok(())
}
fn bytes(h: &hdf5::Group, name: &str, data: &[u8]) -> Result<()> {
    h.new_dataset_builder().with_data(data).create(name)?;
    Ok(())
}
pub fn prepare(p: &Pack) -> Result<(Vec<lpi::Pulse>, Vec<Diagnostic>)> {
    let mut pulses = Vec::new();
    let mut diagnostics = Vec::new();
    for i in 0..p.echo.len() {
        let dc=&p.config.instrument_metadata["dc_offset"];
        let dc=C::new(dc[0].as_f64().unwrap_or(0.),dc[1].as_f64().unwrap_or(0.));
        let echo:Vec<_>=p.echo[i].iter().map(|z|*z-dc).collect();
        let mut transmit=p.transmit[i].clone();
        for z in &mut transmit[p.config.lpi.tx_gate[0]..p.config.lpi.tx_gate[1]]{*z-=dc;}
        let result = if p.config.cancel_debris {
            Some(
                debris::cancel(
                    &echo,
                    &transmit,
                    &p.detections[i],
                    p.config.lpi.sample_rate_hz,
                    p.config.lpi.noise_gate,
                    &p.config.cancellation,
                )
                .with_context(|| format!("cancelling pulse sample {}", p.samples[i]))?,
            )
        } else {
            None
        };
        let (cleaned, usable) = if let Some(c) = result {
            diagnostics.push(c.diagnostic);
            (c.cleaned, c.usable)
        } else {
            (echo, vec![true; p.echo[i].len()])
        };
        pulses.push(lpi::Pulse {
            echo: cleaned,
            transmit,
            pair: p.pair[i],
            usable,
        });
    }
    Ok((pulses, diagnostics))
}
pub fn write(
    path: &Path,
    source: &Path,
    p: &Pack,
    pulses: &[lpi::Pulse],
    diagnostics: &[Diagnostic],
    fit: Option<&lpi::Estimate>,
) -> Result<()> {
    ensure!(!path.exists(), "output exists; use a new product version");
    let tmp = path.with_extension(format!("{}.tmp", std::process::id()));
    let result = (|| -> Result<()> {
        let h = hdf5::File::create(&tmp)?;
        h.new_attr_builder()
            .with_data("ruisdap-isr-v1")
            .create("schema")?;
        h.new_attr_builder()
            .with_data("UTC Unix seconds; half-open [start,end)")
            .create("time_convention")?;
        h.create_group("calibration")?;
        let sources = h.create_group("sources")?;
        let src = sources.create_group(&p.hash[..20])?;
        src.new_attr_builder()
            .with_data(p.hash.as_str())
            .create("sha256")?;
        src.new_attr_builder()
            .with_data(source.to_string_lossy().as_ref())
            .create("path")?;
        let records = h.create_group("records")?;
        let r = records.create_group(&p.hash[..24])?;
        let c = &p.config;
        let metadata = serde_json::json!({"id":&p.hash[..24],"source_id":&p.hash[..20],"source_key":"/","station":c.station,"antenna":c.antenna,"channel":c.channel,"mode":c.mode,"radar_frequency_hz":c.radar_frequency_hz,"start_unix_s":c.start_unix_s,"end_unix_s":c.end_unix_s,"kind":if fit.is_some(){"ion_acf"}else{"voltage_cancellation"},"signal_units":if c.lpi.code_energy_samples.is_empty(){"K; legacy unit-energy pulse normalization"}else{"K km^-1; physical unit-peak pulse energy and range-cell measure removed per pulse"},"complex_convention":"R(l)=x(t)conj(x(t+l)); original ADC data retained in input","cancellation":c.cancel_debris});
        r.new_attr_builder()
            .with_data(metadata.to_string().as_str())
            .create("metadata")?;
        bytes(&r, "configuration_json", &serde_json::to_vec(c)?)?;
        r.new_dataset_builder()
            .with_data(&p.samples)
            .create("pulse_sample")?;
        let qa = r.create_group("debris")?;
        bytes(&qa, "diagnostics_json", &serde_json::to_vec(diagnostics)?)?;
        let ns = p.echo[0].len();
        let cleaned: Vec<_> = pulses
            .iter()
            .flat_map(|p| p.echo.iter().flat_map(|z| [z.re, z.im]))
            .collect();
        array(
            &qa,
            "cleaned_voltage",
            &[pulses.len(), ns, 2],
            &cleaned,
            "ADC",
            "pulse,sample,complex_component",
        )?;
        let mask: Vec<u8> = pulses
            .iter()
            .flat_map(|p| p.usable.iter().map(|x| u8::from(*x)))
            .collect();
        let d = qa
            .new_dataset::<u8>()
            .shape([pulses.len(), ns])
            .create("usable")?;
        d.write_raw(&mask)?;
        if let Some(f) = fit {
            let nr = c.lpi.range_end - c.lpi.range_start;
            let nl = c.lpi.lag_samples.len();
            let ion = r.create_group("ion_line")?;
            if !c.lpi.code_energy_samples.is_empty(){
                array(&ion,"code_energy_samples",&[pulses.len()],&c.lpi.code_energy_samples,"samples at unit peak power","pulse")?;
                array(&ion,"range_cell_km",&[1],&[c.lpi.decimation as f64/c.lpi.sample_rate_hz*299792.458/2.],"km","scalar")?;
            }
            let range: Vec<_> = (c.lpi.range_start..c.lpi.range_end)
                .map(|i| {
                    i as f64 * c.lpi.decimation as f64 / c.lpi.sample_rate_hz * 299792458. / 2.
                })
                .collect();
            array(&r, "range_m", &[nr], &range, "m", "range")?;
            array(
                &ion,
                "lag_s",
                &[nl],
                &c.lpi
                    .lag_samples
                    .iter()
                    .map(|l| *l as f64 / c.lpi.sample_rate_hz)
                    .collect::<Vec<_>>(),
                "s",
                "lag",
            )?;
            for (name, data) in [
                ("acf", &f.acf),
                ("ground_clutter_acf", &f.ground_clutter_acf),
            ] {
                array(
                    &ion,
                    name,
                    &[nr, nl, 2],
                    &data.iter().flat_map(|z| [z.re, z.im]).collect::<Vec<_>>(),
                    if c.lpi.code_energy_samples.is_empty(){"K"}else{"K km^-1; physical pulse-energy and range-cell normalization"},
                    "range,lag,complex_component",
                )?;
            }
            array(
                &ion,
                "variance",
                &[nr, nl],
                &f.variance,
                if c.lpi.code_energy_samples.is_empty(){"K^2; circular complex variance approximation"}else{"K^2 km^-2; circular complex variance approximation"},
                "range,lag",
            )?;
            array(
                &ion,
                "retained_fraction",
                &[nl],
                &f.retained_fraction,
                "1",
                "lag",
            )?;
            array(
                &ion,
                "cholesky_diagonal_ratio",
                &[nl],
                &f.condition_proxy,
                "1; diagnostic, not full condition number",
                "lag",
            )?;
            array(
                &ion,
                "receiver_gain",
                &[1],
                &[f.receiver_gain],
                "ADC^2 K^-1",
                "scalar",
            )?;
            array(
                &ion,
                "system_temperature_k",
                &[1],
                &[f.system_temperature_k],
                "K",
                "scalar",
            )?;
        }
        // Range–Doppler products use accepted cancellation masks. A failed removal
        // cannot silently contribute zero power or a supposedly clean spectrum.
        if !c.range_doppler_samples.is_empty() {
            let tx = &p.transmit[0][c.lpi.tx_gate[0]..c.lpi.tx_gate[1]];
            let mut sum = vec![0.; c.range_doppler_samples.len() * c.range_doppler_fft];
            let mut count = vec![0.; sum.len()];
            for pulse in pulses {
                let a = signal::range_doppler(
                    &pulse.echo,
                    tx,
                    &c.range_doppler_samples,
                    c.range_doppler_fft,
                )?;
                for (ri, &rg) in c.range_doppler_samples.iter().enumerate() {
                    if pulse.usable[rg..rg + tx.len()].iter().all(|x| *x) {
                        for fi in 0..c.range_doppler_fft {
                            let k = ri * c.range_doppler_fft + fi;
                            sum[k] += a[k];
                            count[k] += 1.;
                        }
                    }
                }
            }
            for (v, n) in sum.iter_mut().zip(&count) {
                *v = if *n > 0. { *v / n } else { f64::NAN }
            }
            let rd = r.create_group("range_doppler")?;
            array(
                &rd,
                "spectrum",
                &[c.range_doppler_samples.len(), c.range_doppler_fft],
                &sum,
                "ADC^2; Hann matched-filter power; uncalibrated",
                "range,frequency",
            )?;
            array(
                &rd,
                "count",
                &[c.range_doppler_samples.len(), c.range_doppler_fft],
                &count,
                "1",
                "range,frequency",
            )?;
            let freq: Vec<_> = (0..c.range_doppler_fft)
                .map(|i| {
                    (i as f64 - (c.range_doppler_fft / 2) as f64) * c.lpi.sample_rate_hz
                        / c.range_doppler_fft as f64
                })
                .collect();
            array(&rd, "frequency_hz", &[freq.len()], &freq, "Hz", "frequency")?;
            let range: Vec<_> = c
                .range_doppler_samples
                .iter()
                .map(|x| *x as f64 / c.lpi.sample_rate_hz * 299792458. / 2.)
                .collect();
            array(
                &rd,
                "range_m",
                &[range.len()],
                &range,
                "m; delay relative to transmit template origin",
                "range",
            )?;
        }
        h.flush()?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&tmp);
    }
    result?;
    fs::rename(tmp, path)?;
    Ok(())
}
pub fn now() -> f64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs_f64()
}
