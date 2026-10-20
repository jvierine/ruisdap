//! Read the original multi-directory DigitalRF archive directly into bounded
//! in-memory LPI observations. No RF intermediate files are ever written.
use crate::{
    debris::{Detection, Settings as Cancellation},
    lpi,
    product::{self, Configuration, Pack},
    C,
};
use anyhow::{ensure, Context, Result};
use guisdap_digital_rf::{CacheOptions, MissingSamples, Reader};
use hdf5::{File, Group};
use serde::Deserialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fs,
    path::PathBuf,
};
#[derive(Deserialize)]
pub struct Manifest {
    pub roots: Vec<PathBuf>,
    #[serde(default)]
    pub metadata_roots: Vec<PathBuf>,
    pub start_unix_s: u64,
    pub end_unix_s: u64,
    pub output: PathBuf,
    pub live_directory: Option<PathBuf>,
    #[serde(default="default_cancellation")]
    pub cancel_debris: bool,
}
fn default_cancellation()->bool{true}
struct Stream {
    directories: Vec<PathBuf>,
    cadence: u64,
    subdir: u64,
    prefix: String,
    numerator: u64,
    denominator: u64,
}
struct Metadata {
    streams: BTreeMap<String, Stream>,
    files: VecDeque<(PathBuf, File)>,
}
fn scalar(g: &Group, n: &str) -> Result<f64> {
    let a = g.dataset(n)?.read_raw::<f64>()?;
    ensure!(a.len() == 1, "non-scalar {n}");
    Ok(a[0])
}
fn text_dataset(g: &Group, n: &str) -> Result<String> {
    Ok(g.dataset(n)?.read_raw::<hdf5::types::FixedAscii<256>>()?[0]
        .as_str()
        .to_owned())
}
fn text_attr(h: &File, n: &str) -> Result<String> {
    let a = h.attr(n)?;
    match a.dtype()?.to_descriptor()? {
        hdf5::types::TypeDescriptor::Unsigned(hdf5::types::IntSize::U1) => Ok(String::from_utf8(a.read_raw::<u8>()?)?),
        hdf5::types::TypeDescriptor::VarLenUnicode => Ok(a
            .read_raw::<hdf5::types::VarLenUnicode>()?[0]
            .as_str()
            .to_owned()),
        hdf5::types::TypeDescriptor::VarLenAscii => Ok(a.read_raw::<hdf5::types::VarLenAscii>()?
            [0]
        .as_str()
        .to_owned()),
        _ => Ok(a.read_raw::<hdf5::types::FixedAscii<256>>()?[0]
            .as_str()
            .to_owned()),
    }
}
fn attr(h: &File, n: &str) -> Result<u64> {
    let a = h.attr(n)?.read_raw::<u64>()?;
    ensure!(a.len() == 1, "non-scalar property");
    Ok(a[0])
}
impl Metadata {
    fn open(roots: &[PathBuf]) -> Result<Self> {
        let mut streams: BTreeMap<String, Stream> = BTreeMap::new();
        for root in roots {
            for name in [
                "antenna_control_metadata",
                "powermeter",
                "id_metadata",
                "satellite_detections",
            ] {
                let path = root.join(name);
                if !path.join("dmd_properties.h5").is_file() {
                    continue;
                }
                let h = File::open(path.join("dmd_properties.h5"))?;
                let properties = (
                    attr(&h, "file_cadence_secs")?,
                    attr(&h, "subdir_cadence_secs")?,
                    text_attr(&h, "file_name")?,
                    attr(&h, "sample_rate_numerator")?,
                    attr(&h, "sample_rate_denominator")?,
                );
                ensure!(
                    properties.0 > 0 && properties.1 > 0 && properties.3 > 0 && properties.4 > 0,
                    "invalid DigitalMetadata properties"
                );
                let s = streams.entry(name.into()).or_insert_with(|| Stream {
                    directories: Vec::new(),
                    cadence: properties.0,
                    subdir: properties.1,
                    prefix: properties.2.clone(),
                    numerator: properties.3,
                    denominator: properties.4,
                });
                ensure!(
                    (
                        s.cadence,
                        s.subdir,
                        s.prefix.clone(),
                        s.numerator,
                        s.denominator
                    ) == properties,
                    "incompatible metadata stream {name}"
                );
                let path = path.canonicalize()?;
                if !s.directories.contains(&path) {
                    s.directories.push(path);
                }
            }
        }
        ensure!(
            streams.len() == 4,
            "missing required Millstone DigitalMetadata streams"
        );
        Ok(Self {
            streams,
            files: VecDeque::new(),
        })
    }
    fn files(&mut self, name: &str, begin: u64, end: u64) -> Result<Vec<File>> {
        let s = self.streams.get(name).context("unknown metadata stream")?;
        let mut paths = Vec::new();
        for root in &s.directories {
            for second in (begin / s.cadence * s.cadence..=end / s.cadence * s.cadence)
                .step_by(s.cadence as usize)
            {
                let dir =
                    chrono::DateTime::from_timestamp((second / s.subdir * s.subdir) as i64, 0)
                        .context("invalid date")?;
                paths.push(
                    root.join(dir.format("%Y-%m-%dT%H-%M-%S").to_string())
                        .join(format!("{}@{second}.h5", s.prefix)),
                );
            }
        }
        let mut files = Vec::new();
        for path in paths {
            if let Some(i) = self.files.iter().position(|(p, _)| p == &path) {
                let v = self.files.remove(i).unwrap();
                files.push(v.1.clone());
                self.files.push_back(v);
            } else if path.is_file() {
                let h = File::open(&path)?;
                files.push(h.clone());
                self.files.push_back((path, h));
                while self.files.len() > 64 {
                    self.files.pop_front();
                }
            }
        }
        Ok(files)
    }
    fn padded_rows(&mut self, name: &str, begin: u64, end: u64) -> Result<BTreeMap<u64, Group>> {
        let mut rows = self.rows(name, begin, end)?;
        let stream = self.streams.get(name).context("unknown metadata stream")?;
        let target = u64::try_from(
            u128::from(begin) * u128::from(stream.numerator) / u128::from(stream.denominator),
        )?;
        if rows.range(..=target).next_back().is_some() {
            return Ok(rows);
        }
        // Instrument state is event-driven. Search the preceding actual metadata
        // file, rather than inventing an expiry time for unchanged antenna state.
        let mut candidates = Vec::new();
        let prefix = format!("{}@", stream.prefix);
        let directory_limit =
            chrono::DateTime::from_timestamp((begin / stream.subdir * stream.subdir) as i64, 0)
                .context("invalid date")?
                .format("%Y-%m-%dT%H-%M-%S")
                .to_string();
        for root in &stream.directories {
            for d in fs::read_dir(root)? {
                let d = d?.path();
                if !d.is_dir()
                    || d.file_name().unwrap().to_string_lossy().as_ref() > directory_limit.as_str()
                {
                    continue;
                }
                for f in fs::read_dir(d)? {
                    let p = f?.path();
                    let n = p.file_name().unwrap().to_string_lossy();
                    if let Some(t) = n
                        .strip_prefix(&prefix)
                        .and_then(|v| v.strip_suffix(".h5"))
                        .and_then(|v| v.parse::<u64>().ok())
                    {
                        if t <= begin {
                            candidates.push((t, p));
                        }
                    }
                }
            }
        }
        candidates.sort_by(|a, b| b.0.cmp(&a.0));
        for (_, path) in candidates {
            let f = File::open(&path)?;
            let mut found = None;
            for name in f.member_names()? {
                if let Ok(k) = name.parse::<u64>() {
                    if k <= target && found.as_ref().is_none_or(|(old, _)| k > *old) {
                        found = Some((k, f.group(&name)?));
                    }
                }
            }
            if let Some((k, g)) = found {
                rows.insert(k, g);
                break;
            }
        }
        Ok(rows)
    }
    fn rows(&mut self, name: &str, begin: u64, end: u64) -> Result<BTreeMap<u64, Group>> {
        let mut out = BTreeMap::new();
        for f in self.files(name, begin, end)? {
            for key in f.member_names()? {
                if let Ok(k) = key.parse::<u64>() {
                    ensure!(
                        out.insert(k, f.group(&key)?).is_none(),
                        "overlapping metadata sample {k}"
                    );
                }
            }
        }
        Ok(out)
    }
}
fn power(time: f64, rows: &[(f64, f64)]) -> Option<f64> {
    if rows.is_empty() {
        return None;
    }
    let i = rows.partition_point(|&(t, _)| t < time);
    let nearest = if i == 0 {
        (rows[0].0 - time).abs()
    } else if i == rows.len() {
        (rows[i - 1].0 - time).abs()
    } else {
        (rows[i].0 - time).abs().min((rows[i - 1].0 - time).abs())
    };
    if nearest > 30. {
        return None;
    }
    Some(if i == 0 {
        rows[0].1
    } else if i == rows.len() {
        rows[i - 1].1
    } else {
        let (a, x) = rows[i - 1];
        let (b, y) = rows[i];
        x + (y - x) * (time - a) / (b - a)
    })
}
fn channel_delay(echo: &[C], tx: &[C], gate: [usize; 2]) -> f64 {
    let a = &echo[gate[0] - 40..gate[1] + 80];
    let b = &tx[gate[0] - 40..gate[1] + 80];
    let n = a.len() as isize;
    let scores: Vec<_> = (-40isize..=40)
        .map(|lag| {
            let lo = 0.max(-lag);
            let hi = n.min(n - lag);
            (lo..hi)
                .map(|i| b[i as usize].conj() * a[(i + lag) as usize])
                .sum::<C>()
                .norm()
        })
        .collect();
    let best = (0..scores.len())
        .max_by(|&i, &j| scores[i].total_cmp(&scores[j]))
        .unwrap();
    let delta = if best > 0 && best + 1 < scores.len() {
        let d = scores[best - 1] - 2. * scores[best] + scores[best + 1];
        if d != 0. {
            0.5 * (scores[best - 1] - scores[best + 1]) / d
        } else {
            0.
        }
    } else {
        0.
    };
    best as f64 - 40. + delta
}
pub struct Archive {
    rf: Reader,
    metadata: Metadata,
    roots: Vec<PathBuf>,
    kernel_sha256:String,
    inherited_kernel_sha256:Option<String>,
    cancel_debris:bool,
}
impl Archive {
    pub fn open(m: &Manifest) -> Result<Self> {
        ensure!(!m.roots.is_empty(), "no DigitalRF roots");
        let mut rf = Reader::open(&m.roots)?;
        rf.set_cache_options(CacheOptions {
            minimum_bytes: 64 * 1024 * 1024,
            max_open_files: 32,
        })?;
        for name in ["zenith-l", "tx-h"] {
            let p = rf.properties(name)?;
            ensure!(
                p.sample_rate_numerator == 1000000
                    && p.sample_rate_denominator == 1
                    && p.num_subchannels == 1
                    && p.is_complex,
                "unsupported Millstone RF layout"
            );
        }
        let metadata_roots = if m.metadata_roots.is_empty() {
            m.roots.iter().map(|r| r.join("metadata")).collect()
        } else {
            m.metadata_roots.clone()
        };
        let kernel_sha256=format!("{:x}",Sha256::digest(fs::read(std::env::current_exe()?)?));
        let inherited_file=m.output.parent().unwrap_or(&m.output).join("kernel.h5");
        let inherited_kernel_sha256=if inherited_file.is_file(){Some(text_attr(&File::open(inherited_file)?,"binary_sha256")?)}else{None};
        Ok(Self {
            kernel_sha256,inherited_kernel_sha256,cancel_debris:m.cancel_debris,
            rf,
            metadata: Metadata::open(&metadata_roots)?,
            roots: m.roots.clone(),
        })
    }
    pub fn observations(&mut self, begin: u64) -> Result<(Vec<Pack>, Vec<Value>, Vec<u64>)> {
        let ants = self
            .metadata
            .padded_rows("antenna_control_metadata", begin, begin + 11)?;
        let mut antenna = BTreeMap::new();
        for (k, g) in ants {
            antenna.insert(
                k,
                text_dataset(&g, "rx_antenna")? == "ZENITH"
                    && text_dataset(&g, "tx_antenna")? == "ZENITH",
            );
        }
        let powers = self
            .metadata
            .rows("powermeter", begin.saturating_sub(30), begin + 41)?;
        let power_rows: Vec<_> = powers
            .iter()
            .map(|(&k, g)| Ok((k as f64, scalar(g, "zenith_power")?)))
            .collect::<Result<_>>()?;
        let ids = self
            .metadata
            .rows("id_metadata", begin.saturating_sub(1), begin + 11)?;
        let ids: BTreeMap<_, _> = ids
            .iter()
            .map(|(&k, g)| Ok((k, scalar(g, "sweepid")? as u64)))
            .collect::<Result<_>>()?;
        let catalogues =
            self.metadata
                .files("satellite_detections", begin.saturating_sub(1), begin + 11)?;
        let mut observations = Vec::new();
        let mut gaps = Vec::new();
        let mut candidates = Vec::new();
        for family in ["coded", "long"] {
            let gate = if family == "coded" {
                [76, 624]
            } else {
                [76, 645]
            };
            let mut source = BTreeMap::new();
            for (&k, &code) in &ids {
                if !(if family == "coded" {
                    (1..=32).contains(&code)
                } else {
                    code == 300
                }) {
                    continue;
                }
                if !antenna.range(..=k).next_back().is_some_and(|(_, v)| *v) {
                    continue;
                }
                let Some(watts) = power(k as f64 / 1e6, &power_rows) else {
                    if k >= begin * 1000000 && k < (begin + 10) * 1000000 {
                        gaps.push(json!({"pulse_sample":k,"reason":"no recorded transmit power within 30 seconds"}));
                    }
                    continue;
                };
                if watts < 400000. {
                    continue;
                }
                if k >= begin * 1000000 && k < (begin + 10) * 1000000 {
                    candidates.push(k);
                }
                let mut data = Vec::new();
                let mut missing = false;
                for channel in ["zenith-l", "tx-h"] {
                    match self.rf.read_vector(k, 10000, channel, Some(0)) {
                        Ok(v) => data.push(v.into_complex64()),
                        Err(e) => {
                            if let Some(g) = e.downcast_ref::<MissingSamples>() {
                                if k >= begin * 1000000 && k < (begin + 10) * 1000000 {
                                    gaps.push(json!({"pulse_sample":k,"channel":channel,"missing_samples":g.gaps.iter().map(|r|[r.start,r.end]).collect::<Vec<_>>() }));
                                }
                                missing = true;
                                break;
                            } else {
                                return Err(e);
                            }
                        }
                    }
                }
                if !missing {
                    source.insert(k, (data.remove(0), data.remove(0), watts));
                }
            }
            let keys: Vec<_> = source.keys().copied().collect();
            let mut pair = BTreeMap::new();
            for (i, &k) in keys.iter().enumerate() {
                for j in [i.checked_sub(1), Some(i + 1).filter(|&j| j < keys.len())]
                    .into_iter()
                    .flatten()
                {
                    let other = keys[j];
                    if ids[&k] == ids[&other] && k.abs_diff(other) < 40000 {
                        pair.insert(k, other);
                        break;
                    }
                }
            }
            let science: Vec<_> = keys
                .iter()
                .copied()
                .filter(|&k| {
                    k >= begin * 1000000 && k < (begin + 10) * 1000000 && pair.contains_key(&k)
                })
                .collect();
            for &k in &keys {
                if k >= begin * 1000000 && k < (begin + 10) * 1000000 && !pair.contains_key(&k) {
                    gaps.push(
                        json!({"pulse_sample":k,"reason":"no real same-code partner within 40 ms"}),
                    );
                }
            }
            if science.len() < 16 {
                if !science.is_empty() {
                    gaps.push(json!({"pulse_samples":science,"reason":"fewer than sixteen paired science pulses in integration"}));
                }
                continue;
            }
            let mut stored: BTreeSet<_> = science.iter().copied().collect();
            for k in &science {
                stored.insert(pair[k]);
            }
            let keys: Vec<_> = stored.into_iter().collect();
            for k in &keys {
                if !science.contains(k) {
                    let partner = *science
                        .iter()
                        .find(|&&p| pair[&p] == *k)
                        .context("missing real partner")?;
                    pair.insert(*k, partner);
                }
            }
            let index: BTreeMap<_, _> = keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
            let mut echo = Vec::new();
            let mut transmit = Vec::new();
            let mut delays = Vec::new();
            let mut detections = Vec::new();
            let mut cfar = Vec::new();
            let mut ptx = Vec::new();
            for &k in &keys {
                let (y, tx, w) = source.remove(&k).unwrap();
                let delay = channel_delay(&y, &tx, gate);
                let mut template = vec![C::default(); 10000];
                template[gate[0]..gate[1]].copy_from_slice(&y[gate[0]..gate[1]]);
                let mut found = Vec::new();
                for cat in &catalogues {
                    if !cat.link_exists(&k.to_string()) {
                        continue;
                    }
                    let g = cat.group(&k.to_string())?;
                    if text_dataset(&g, "channel")? != "zenith-l" {
                        continue;
                    }
                    let raw = g.dataset("raw_delay_sample")?.read_raw::<f64>()?;
                    let dop = g.dataset("doppler_hz")?.read_raw::<f64>()?;
                    ensure!(raw.len() == dop.len(), "invalid CFAR detection dimensions");
                    for (&r, &f) in raw.iter().zip(&dop) {
                        found.push(Detection {
                            delay_samples: r - gate[0] as f64 - delay,
                            doppler_hz: f,
                        });
                        cfar.push(json!({"pulse_sample":k,"raw_delay_sample":r,"doppler_hz":f}));
                    }
                }
                echo.push(y);
                transmit.push(template);
                delays.push(delay);
                detections.push(found);
                ptx.push(w);
            }
            let selected_power =
                science.iter().map(|k| ptx[index[k]]).sum::<f64>() / science.len() as f64;
            let lags: Vec<_> = (10..=450).step_by(10).collect();
            let minimum = lags
                .iter()
                .map(|&l| ((800. - 111. - 480. * 0.2 + l as f64) / 30.).round() as usize)
                .collect();
            let config = Configuration {
                station: "Millstone Hill".into(),
                antenna: "ZENITH".into(),
                channel: "zenith-l".into(),
                mode: if family == "coded" {
                    "alternating-code-16"
                } else {
                    "long-pulse-480us"
                }
                .into(),
                radar_frequency_hz: 440.2e6,
                start_unix_s: science[0] as f64 / 1e6,
                end_unix_s: (science[science.len() - 1] + 10000) as f64 / 1e6,
                lpi: lpi::Settings {
                    code_energy_samples: vec![480.; keys.len()],
                    science_pulse: keys.iter().map(|k| science.contains(k)).collect(),
                    sample_rate_hz: 1e6,
                    decimation: 30,
                    lag_samples: lags,
                    minimum_range_index: minimum,
                    range_start: 0,
                    range_end: 240,
                    measurement_start: 27,
                    measurement_end: 273,
                    tx_gate: gate,
                    echo_gate: if family == "coded" {
                        [1000, 8200]
                    } else {
                        [1000, 7700]
                    },
                    noise_gate: if family == "coded" {
                        [7700, 8200]
                    } else {
                        [7200, 7700]
                    },
                    injection_gate: if family == "coded" {
                        [8400, 8850]
                    } else {
                        [7800, 8371]
                    },
                    filter_half: 100,
                    passband_hz: 18000.,
                    outlier_threshold: 10.,
                    background_column: true,
                },
                cancellation: Cancellation::default(),
                cancel_debris: self.cancel_debris,
                injection_temperature_k: 1172.,
                range_doppler_fft: 1024,
                range_doppler_samples: Vec::new(),
                instrument_metadata: json!({"transmitter_w":selected_power,"pulse_transmitter_w":ptx,"bit_length_us":if family=="coded"{30.}else{480.},"code_bits":if family=="coded"{16}else{1},"dc_offset":[-0.212,-0.221],"tx_to_receiver_delay_samples":delays,"cfar":cfar,"azimuth_deg":0.,"elevation_deg":90.,"science_pulse_count":science.len(),"partner_pulse_count":keys.len()-science.len(),"digital_rf_roots":self.roots,"source_gaps":gaps,"source_hash_kind":"SHA256 of configuration, exact pulse keys and original numerical receive samples","settings_reference":"historical 10 s decode; 30 us range basis; 45 measured lags; full 7200 us delay grid; actual waveform-specific gates"}),
            };
            let mut hash = Sha256::new();
            hash.update(serde_json::to_vec(&config)?);
            for k in &keys {
                hash.update(k.to_le_bytes());
            }
            for row in &echo {
                for z in row {
                    hash.update(z.re.to_le_bytes());
                    hash.update(z.im.to_le_bytes());
                }
            }
            observations.push(Pack {
                config,
                echo,
                transmit,
                pair: keys.iter().map(|k| index[&pair[k]]).collect(),
                detections,
                samples: keys,
                hash: format!("{:x}", hash.finalize()),
            });
        }
        candidates.sort_unstable();
        ensure!(
            candidates.windows(2).all(|w| w[0] != w[1]),
            "duplicate candidate science pulse"
        );
        Ok((observations, gaps, candidates))
    }
}
/// Atomic per-window audit is the resume checkpoint. Only reduced ACF products,
/// masks and QA are stored. Raw/cleaned voltage samples remain in DigitalRF.
pub fn run_window(archive: &mut Archive, m: &Manifest, begin: u64) -> Result<()> {
    fs::create_dir_all(m.output.join("audit"))?;
    fs::create_dir_all(m.output.join("products"))?;
    let audit = m.output.join("audit").join(format!("{begin}.h5"));
    if audit.exists() {
        let h = File::open(&audit)?;
        let manifest: Value =
            serde_json::from_slice(&h.dataset("outputs_json")?.read_raw::<u8>()?)?;
        for row in manifest.as_array().context("invalid checkpoint")? {
            let path = PathBuf::from(row["path"].as_str().context("invalid output path")?);
            ensure!(
                format!("{:x}", Sha256::digest(fs::read(path)?))
                    == row["sha256"].as_str().unwrap_or(""),
                "reduced output changed"
            );
        }
        return Ok(());
    }
    let (observations, mut gaps, candidates) = archive.observations(begin)?;
    let mut outputs = Vec::new();
    let mut science = Vec::new();
    let mut calibration_rejected_science = Vec::new();
    for p in observations {
        let family = if p.config.mode.starts_with("alternating") {
            "coded"
        } else {
            "long"
        };
        let path = m
            .output
            .join("products")
            .join(format!("{begin}-{family}.h5"));
        if path.exists() {
            let h = File::open(&path)?;
            ensure!(
                h.group("records")?.link_exists(&p.hash[..24]),
                "existing output differs from source observation"
            );
        } else {
            let (pulses, diagnostics) = product::prepare(&p)?;
            if pulses.iter().filter(|p|p.usable.iter().any(|&x|x)).count()<16 {
                let selected:Vec<_>=p.samples.iter().zip(&p.config.lpi.science_pulse).filter(|(_,selected)|**selected).map(|(&k,_)|k).collect();
                calibration_rejected_science.extend_from_slice(&selected);
                let calibration:Vec<_>=p.samples.iter().zip(&pulses).filter(|(_,pulse)|pulse.usable.iter().any(|&x|x)).map(|(&k,_)|k).collect();
                gaps.push(json!({"reason":"insufficient_receiver_calibration","mode":p.config.mode,"source_sha256":p.hash,"science_pulse_samples":selected,"pulse_samples":p.samples,"calibration_pulse_samples":calibration,"configuration":p.config,"cancellation_diagnostics":diagnostics}));
                continue;
            }
            let fit = match lpi::estimate(&pulses, &p.config.lpi, p.config.injection_temperature_k) {
                Ok(fit)=>fit,
                Err(error)=>{
                    if let Some(q)=error.downcast_ref::<lpi::InvalidReceiverGain>() {
                        let selected:Vec<_>=p.samples.iter().zip(&p.config.lpi.science_pulse).filter(|(_,selected)|**selected).map(|(&k,_)|k).collect();
                        calibration_rejected_science.extend_from_slice(&selected);
                        let calibration:Vec<_>=p.samples.iter().zip(&pulses).filter(|(_,pulse)|pulse.usable.iter().any(|&x|x)).map(|(&k,_)|k).collect();
                        gaps.push(json!({"reason":"invalid_receiver_gain","mode":p.config.mode,"source_sha256":p.hash,"science_pulse_samples":selected,"pulse_samples":p.samples,"calibration_pulse_samples":calibration,"noise_power":q.noise_power,"injection_power":q.injection_power,"receiver_gain":q.receiver_gain,"injection_temperature_k":q.injection_temperature_k,"configuration":p.config,"cancellation_diagnostics":diagnostics}));
                        continue;
                    }
                    return Err(error);
                }
            };
            product::write_options(
                &path,
                &m.roots[0],
                &p,
                &pulses,
                &diagnostics,
                Some(&fit),
                false,
                Some(&archive.kernel_sha256),
            )?;
        }
        for (&k, &selected) in p.samples.iter().zip(&p.config.lpi.science_pulse) {
            if selected {
                science.push(k);
            }
        }
        let product_kernel={let h=File::open(&path)?;if h.attr_names()?.iter().any(|n|n=="processing_kernel_sha256"){Some(text_attr(&h,"processing_kernel_sha256")?)}else{archive.inherited_kernel_sha256.clone()}};
        outputs
            .push(json!({"path":path,"sha256":format!("{:x}",Sha256::digest(fs::read(&path)?)),"kernel_sha256":product_kernel}));
    }
    science.sort_unstable();
    calibration_rejected_science.sort_unstable();
    ensure!(calibration_rejected_science.windows(2).all(|w|w[0]!=w[1]) && calibration_rejected_science.iter().all(|k|science.binary_search(k).is_err() && candidates.contains(k)),"invalid receiver-calibration exclusion ledger");
    ensure!(
        science.windows(2).all(|w| w[0] != w[1]),
        "duplicate science pulse"
    );
    let tmp = audit.with_extension(format!("{}.tmp", std::process::id()));
    {
        let h = File::create(&tmp)?;
        h.new_attr_builder().with_data(archive.kernel_sha256.as_str()).create("processing_kernel_sha256")?;
        h.new_dataset_builder()
            .with_data(&candidates)
            .create("candidate_science_pulse_sample")?;
        h.new_dataset_builder()
            .with_data(&science)
            .create("science_pulse_sample")?;
        h.new_dataset_builder().with_data(&calibration_rejected_science).create("calibration_rejected_science_pulse_sample")?;
        h.new_dataset_builder()
            .with_data(&serde_json::to_vec(&outputs)?)
            .create("outputs_json")?;
        h.new_dataset_builder()
            .with_data(&serde_json::to_vec(&gaps)?)
            .create("gaps_json")?;
        h.new_attr_builder()
            .with_data("direct multi-directory DigitalRF; no RF staging")
            .create("input_mode")?;
        h.flush()?;
    }
    fs::rename(tmp, audit)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cancellation_controls_are_opt_out_and_keep_source_selection() {
        let base=json!({"roots":["original-rf"],"start_unix_s":100,"end_unix_s":110,"output":"reduced"});
        let production:Manifest=serde_json::from_value(base.clone()).unwrap();
        assert!(production.cancel_debris);
        let mut control=base;
        control["cancel_debris"]=json!(false);
        let control:Manifest=serde_json::from_value(control).unwrap();
        assert!(!control.cancel_debris);
        assert_eq!(production.roots,control.roots);
        assert_eq!(production.start_unix_s,control.start_unix_s);
        assert_eq!(production.end_unix_s,control.end_unix_s);
    }
    #[test]
    fn unchanged_instrument_state_has_no_arbitrary_age_limit() {
        let tmp = tempfile::tempdir().unwrap();
        let directory = tmp.path().join("1970-01-01T00-00-00");
        fs::create_dir(&directory).unwrap();
        {
            let h = File::create(directory.join("ant@0.h5")).unwrap();
            h.create_group("5000000").unwrap();
        }
        {
            let h = File::create(directory.join("ant@600.h5")).unwrap();
            h.create_group("600000001").unwrap();
        }
        let mut streams = BTreeMap::new();
        streams.insert(
            "antenna".into(),
            Stream {
                directories: vec![tmp.path().to_owned()],
                cadence: 60,
                subdir: 3600,
                prefix: "ant".into(),
                numerator: 1000000,
                denominator: 1,
            },
        );
        let mut m = Metadata {
            streams,
            files: VecDeque::new(),
        };
        let rows = m.padded_rows("antenna", 500, 511).unwrap();
        assert_eq!(rows.range(..=500000000).next_back().unwrap().0, &5000000);
        let rows = m.padded_rows("antenna", 615, 626).unwrap();
        assert_eq!(rows.range(..=615000000).next_back().unwrap().0, &600000001);
    }
}
