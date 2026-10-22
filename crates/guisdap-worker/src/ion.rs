//! BEATA ion fitting orchestration. Legacy numerical responses are read as data,
//! never executed. Each observation owns its calibration and each job its state.
mod atmosphere;
mod stream;
use super::{load, metadata};
use anyhow::{ensure, Context, Result};
use guisdap_core::{altitude_km, FrameMeta};
use guisdap_ion::{
    fit::{FitProblem, FitResult, FitSettings, Measurement},
    spectrum::{IonSpecies, PlasmaState, SpectrumModel},
    Complex as C,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    #[serde(default)]
    pub atmosphere: Option<atmosphere::AtmosphereConfig>,
    pub quicklook_frames: PathBuf,
    pub experiment_model: PathBuf,
    pub output: PathBuf,
    pub export: PathBuf,
    pub integration_s: f64,
    pub altitude_edges_km: Vec<f64>,
    pub composition_transition_km: f64,
    pub composition_scale_km: f64,
    pub tie_temperature_below_km: f64,
    pub fit_collisions_below_km: f64,
    pub antenna_efficiency: f64,
    pub calibration_multiplier: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub calibration_provenance: Option<serde_json::Value>,
    pub threads: usize,
    #[serde(default)]
    pub streaming: bool,
    #[serde(default)]
    pub native_range_below_km: Option<f64>,
    #[serde(default)]
    pub model_selection_below_km: Option<f64>,
    #[serde(default)]
    pub warm_start_index: Option<PathBuf>,
    #[serde(default)]
    pub native_slant_edges_km: Vec<f64>,
    #[serde(default)]
    pub live_export: Option<PathBuf>,
}
#[derive(Clone)]
struct Group {
    lag_s: f64,
    width_km: f64,
    range_km: f64,
    step_km: f64,
    address: usize,
    count: usize,
    stride: usize,
    nd: f64,
    peak: f64,
    kernel: Vec<C>,
}
struct Experiment {
    groups: Vec<Group>,
    frequency_hz: Vec<f64>,
    radar_hz: f64,
    gain: f64,
    ap0: f64,
    dt_s: f64,
    hash: String,
}
struct Record {
    meta: FrameMeta,
    calibration_k_per_count: f64,
    radar_factor: f64,
    data: Vec<C>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Gate {
    #[serde(default)]
    pub prescribed_atmosphere: Option<atmosphere::Prescribed>,
    pub altitude_km: f64,
    pub lower_km: f64,
    pub upper_km: f64,
    pub response_width_max_km: f64,
    pub first_lag_signal_to_noise: Option<f64>,
    pub oxygen_fraction: f64,
    pub model: SpectrumModel,
    pub settings: FitSettings,
    pub result: Option<FitResult>,
    pub status: String,
    pub lag_us: Vec<f64>,
    pub observed: Vec<[f64; 2]>,
    pub sigma: Vec<[f64; 2]>,
    pub addresses: Vec<Vec<usize>>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Profile {
    pub start_unix: f64,
    pub end_unix: f64,
    pub source_count: usize,
    pub sources: Vec<FrameMeta>,
    pub gates: Vec<Gate>,
    #[serde(default)]
    pub processing_state: String,
    #[serde(default)]
    pub detail_url: Option<String>,
    #[serde(default)]
    pub raw_ne: Vec<Option<[f64; 2]>>,
    #[serde(default)]
    pub notes: Vec<String>,
}
#[derive(Serialize)]
struct Output<'a> {
    schema: &'static str,
    config: &'a Config,
    experiment_sha256: &'a str,
    uncertainty: &'static str,
    limitations: Vec<&'static str>,
    profiles: &'a [Profile],
}
fn numbers<'a>(m: &'a matfile::MatFile, name: &str) -> Result<&'a [f64]> {
    let a = m
        .find_by_name(name)
        .with_context(|| format!("missing model field {name}"))?;
    if let matfile::NumericData::Double { real, .. } = a.data() {
        Ok(real)
    } else {
        anyhow::bail!("expected double model field {name}")
    }
}
impl Experiment {
    fn read(path: &Path) -> Result<Self> {
        let bytes = fs::read(path)?;
        let hash = format!("{:x}", Sha256::digest(&bytes));
        ensure!(
            hash == "27180ca657512d2bca042824431f688e8b3fa66980bd0916c77aab9c9391f4ec",
            "unsupported BEATA response model SHA-256: {hash}"
        );
        let m = matfile::MatFile::parse(bytes.as_slice())?;
        let dt_s = numbers(&m, "p_dtau")?[0] * 1e-6;
        let radar_hz = numbers(&m, "ch_fradar")?[0];
        let om0 = numbers(&m, "p_om0")?[0];
        let frequency_hz: Vec<_> = numbers(&m, "p_om")?
            .iter()
            .map(|w| w * om0 / (2. * std::f64::consts::PI))
            .collect();
        let nf = frequency_hz.len();
        let df: Vec<_> = (0..nf)
            .map(|i| (frequency_hz[(i + 1).min(nf - 1)] - frequency_hz[i.saturating_sub(1)]) / 2.)
            .collect();
        let ra = numbers(&m, "lpg_ra")?;
        let nt = numbers(&m, "lpg_nt")?;
        let ri = numbers(&m, "lpg_ri")?;
        let nd = numbers(&m, "lpg_ND")?;
        let h = numbers(&m, "lpg_h")?;
        let dt = numbers(&m, "lpg_dt")?;
        let width = numbers(&m, "lpg_w")?;
        let lag = numbers(&m, "lpg_lag")?;
        let kind = numbers(&m, "lpg_bcs")?;
        let bac = numbers(&m, "lpg_bac")?;
        let (wr, wi) = match m
            .find_by_name("lpg_wom")
            .context("missing response")?
            .data()
        {
            matfile::NumericData::Double {
                real,
                imag: Some(imag),
            } => (real, imag),
            _ => anyhow::bail!("expected complex spectral response"),
        };
        let ng = ra.len();
        ensure!(wr.len() == ng * nf, "wrong response shape");
        let scale = dt_s * 299792458. / 2000.;
        let mut coverage = vec![false; 25531 - 480];
        let mut groups = vec![];
        for i in 0..ng {
            if kind[i] != 115. || !(480. ..25531.).contains(&ra[i]) {
                continue;
            }
            ensure!(
                bac[i] == 0. && nd[i] > 0.,
                "unsupported ion lag background/calibration"
            );
            let peak = (0..nf)
                .map(|f| wr[i + ng * f] / nd[i])
                .fold(f64::NEG_INFINITY, f64::max);
            ensure!(peak > 0., "invalid spectral response peak");
            let kernel = (0..nf)
                .map(|f| C::new(wr[i + ng * f], wi[i + ng * f]) * (df[f] / (nd[i] * peak)))
                .collect();
            for j in 0..nt[i] as usize {
                let a = ra[i] as usize + j * ri[i] as usize;
                ensure!(a < 25531 && !coverage[a - 480], "ion layout overlap");
                coverage[a - 480] = true;
            }
            groups.push(Group {
                lag_s: lag[i] * dt_s,
                width_km: width[i] * scale,
                range_km: h[i] * scale,
                step_km: dt[i] * scale,
                address: ra[i] as usize,
                count: nt[i] as usize,
                stride: ri[i] as usize,
                nd: nd[i],
                peak,
                kernel,
            });
        }
        ensure!(coverage.iter().all(|x| *x), "incomplete ion lag coverage");
        Ok(Self {
            groups,
            frequency_hz,
            radar_hz,
            gain: numbers(&m, "ch_gain")?[0],
            ap0: numbers(&m, "vcg_Ap")?[0],
            dt_s,
            hash,
        })
    }
    fn radar_factor(&self, tx: f64, cfg: &Config) -> f64 {
        let lambda = 299792458. / self.radar_hz;
        (2.8179403262e-15_f64).powi(2)
            * self.gain
            * lambda
            * lambda
            * 0.460
            * cfg.antenna_efficiency
            * (299792458. * self.dt_s / 2.)
            * tx
            / (cfg.calibration_multiplier
                * 4.
                * std::f64::consts::PI
                * 1.380649e-23
                * (self.ap0 / self.dt_s))
    }
}
fn read_record(meta: &FrameMeta, exp: &Experiment, cfg: &Config) -> Result<Record> {
    let path = Path::new(&meta.source);
    let (bytes, m) = load(path)?;
    let checked = metadata(path, &bytes, &m)?;
    ensure!(
        checked.source_sha256 == meta.source_sha256,
        "source changed since quick look"
    );
    let par = &m["d_parbl"].real;
    let data = m["d_data"].complex();
    ensure!(data.len() == 196493, "unexpected raw layout");
    let background = data[..16].iter().map(|x| x.re).sum::<f64>() / 16. / 32.;
    let cal = data[16..32].iter().map(|x| x.re).sum::<f64>() / 16. / 32.;
    ensure!(
        cal > background && par[20] > 0. && par[7] > 0.,
        "invalid calibration"
    );
    Ok(Record {
        meta: checked,
        calibration_k_per_count: par[20] / (cal - background),
        radar_factor: exp.radar_factor(par[7], cfg),
        data,
    })
}
fn fit_gate(
    records: &[Record],
    exp: &Experiment,
    cfg: &Config,
    lo: f64,
    hi: f64,
    seeds: &[PlasmaState],
) -> Result<Gate> {
    let native_bounds = cfg
        .altitude_edges_km
        .iter()
        .position(|h| (*h - lo).abs() < 1e-8)
        .and_then(|i| cfg.native_slant_edges_km.get(i..i + 2));
    let (lo, hi) = if let Some(r) = native_bounds {
        let elevation =
            records.iter().map(|r| r.meta.elevation_deg).sum::<f64>() / records.len() as f64;
        (altitude_km(r[0], elevation), altitude_km(r[1], elevation))
    } else {
        (lo, hi)
    };
    let altitude = (lo + hi) / 2.;
    let atmosphere = cfg
        .atmosphere
        .as_ref()
        .context("prescribed MSIS/IRI atmosphere configuration required")?;
    ensure!(
        cfg.model_selection_below_km.is_none(),
        "model selection must not override prescribed collisions/composition"
    );
    let unix = records
        .iter()
        .map(|r| (r.meta.start_unix + r.meta.end_unix) / 2.)
        .sum::<f64>()
        / records.len() as f64;
    let prescribed = atmosphere.evaluate(unix, altitude)?;
    let oxygen = 1. - prescribed.molecular_fraction;
    let model = SpectrumModel {
        radar_frequency_hz: exp.radar_hz,
        ions: vec![
            IonSpecies {
                mass_amu: 16.,
                fraction: oxygen,
            },
            IonSpecies {
                mass_amu: 30.5,
                fraction: 1. - oxygen,
            },
        ],
        electron_collision_s: 0.,
    };
    let settings = FitSettings {
        max_iterations: 160,
        tie_electron_ion_temperature: altitude < cfg.tie_temperature_below_km,
        fit_collision_rate: false,
        fixed_collision_s: prescribed.collision_s,
        ..Default::default()
    };
    let mut measurements = vec![];
    let mut addresses = vec![];
    let mut samples = vec![vec![]; records.len()];
    let mut maxwidth: f64 = 0.;
    for g in &exp.groups {
        // Each point is selected by its measured response centroid. The response width
        // is retained; this is a local homogeneous fit, not a global range inversion.
        let selected: Vec<_> = (0..g.count)
            .filter(|&j| {
                let r = g.range_km + j as f64 * g.step_km;
                if let Some(bounds) = native_bounds {
                    r >= bounds[0] - 1e-8 && r < bounds[1] - 1e-8
                } else {
                    records.iter().all(|rec| {
                        let h = altitude_km(r, rec.meta.elevation_deg);
                        h >= lo && h < hi
                    })
                }
            })
            .collect();
        ensure!(
            native_bounds.is_none() || selected.len() <= 1,
            "native range selection combined {} samples at lag={} us, step={} km, bin={lo}..{hi}",
            selected.len(),
            g.lag_s * 1e6,
            g.step_km
        );
        if selected.is_empty() {
            continue;
        }
        let mut values = vec![];
        for rec in records {
            let v = selected
                .iter()
                .map(|&j| {
                    let r = (g.range_km + j as f64 * g.step_km) * 1000.;
                    rec.data[g.address + j * g.stride]
                        * (rec.calibration_k_per_count * r * r / (g.nd * g.peak * rec.radar_factor))
                })
                .sum::<C>()
                / selected.len() as f64;
            values.push(v);
        }
        let n = values.len() as f64;
        let mean = values.iter().copied().sum::<C>() / n;
        let vr = values.iter().map(|v| (v.re - mean.re).powi(2)).sum::<f64>() / (n * (n - 1.));
        let vi = values.iter().map(|v| (v.im - mean.im).powi(2)).sum::<f64>() / (n * (n - 1.));
        ensure!(
            vr > 0. && vi > 0. && vr.is_finite() && vi.is_finite(),
            "undefined empirical variance"
        );
        for (sample, v) in samples.iter_mut().zip(values) {
            sample.push(v);
        }
        measurements.push(Measurement {
            lag_s: g.lag_s,
            value: mean,
            sigma: [vr.sqrt(), vi.sqrt()],
            kernel: g.kernel.clone(),
        });
        addresses.push(selected.iter().map(|&j| g.address + j * g.stride).collect());
        maxwidth = maxwidth.max(g.width_km);
    }
    let mut gate = Gate {
        prescribed_atmosphere: Some(prescribed),
        altitude_km: altitude,
        lower_km: lo,
        upper_km: hi,
        response_width_max_km: maxwidth,
        first_lag_signal_to_noise: None,
        oxygen_fraction: oxygen,
        model: model.clone(),
        settings: settings.clone(),
        result: None,
        status: "insufficient_lags".into(),
        lag_us: measurements.iter().map(|m| m.lag_s * 1e6).collect(),
        observed: measurements
            .iter()
            .map(|m| [m.value.re, m.value.im])
            .collect(),
        sigma: measurements.iter().map(|m| m.sigma).collect(),
        addresses,
    };
    if measurements.len() < 8 {
        return Ok(gate);
    }
    let first: Vec<_> = measurements
        .iter()
        .enumerate()
        .filter(|(_, m)| m.lag_s < 15e-6)
        .map(|(i, _)| i)
        .collect();
    if !first.is_empty() {
        let x: Vec<_> = samples
            .iter()
            .map(|s| first.iter().map(|&i| s[i].re).sum::<f64>() / first.len() as f64)
            .collect();
        let mean = x.iter().sum::<f64>() / x.len() as f64;
        let se = (x.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
            / (x.len() * (x.len() - 1)) as f64)
            .sqrt();
        gate.first_lag_signal_to_noise = Some(mean / se);
    }
    let estimate = measurements
        .iter()
        .filter(|m| m.lag_s < 15e-6)
        .map(|m| m.value.re * 2.)
        .sum::<f64>()
        / measurements
            .iter()
            .filter(|m| m.lag_s < 15e-6)
            .count()
            .max(1) as f64;
    let problem = FitProblem {
        model,
        frequency_hz: exp.frequency_hz.clone(),
        measurements,
        settings,
    };
    let initial = PlasmaState {
        electron_density_m3: estimate.clamp(1e9, 1e13),
        ion_temperature_k: if altitude < 130. { 300. } else { 1000. },
        electron_temperature_k: if altitude < 130. { 300. } else { 1800. },
        velocity_m_s: 0.,
        ion_collision_s: if altitude < 110. { 3000. } else { 300. },
    };
    let select_model = cfg.model_selection_below_km.is_some_and(|h| altitude < h);
    let mut variants = vec![problem.clone()];
    if select_model {
        if problem.settings.fit_collision_rate {
            let mut zero = problem.clone();
            zero.settings.fit_collision_rate = false;
            zero.settings.fixed_collision_s = 0.;
            variants.push(zero);
        }
        if !problem.settings.tie_electron_ion_temperature {
            let tied: Vec<_> = variants
                .iter()
                .map(|p| {
                    let mut p = p.clone();
                    p.settings.tie_electron_ion_temperature = true;
                    p
                })
                .collect();
            variants.extend(tied);
        }
    }
    let mut candidates = vec![];
    for candidate in variants {
        let mut starts = Vec::new();
        for ratio in [0.7, 1.5] {
            let mut p = initial;
            p.ion_temperature_k *= ratio;
            p.electron_temperature_k *= ratio;
            starts.push(p);
        }
        starts.extend_from_slice(seeds);
        for p in starts {
            if let Ok(fit) = candidate.solve(p) {
                let penalty = if select_model {
                    fit.free_parameter_indices.len() as f64
                        * (2. * candidate.measurements.len() as f64).ln()
                } else {
                    0.
                };
                if fit.chi_square.is_finite() {
                    candidates.push((fit.chi_square + penalty, candidate.clone(), fit));
                }
            }
        }
    }
    if !select_model {
        let suspicious = candidates
            .iter()
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .is_none_or(|(_, _, f)| {
                f.status != "converged"
                    || f.reduced_chi_square > 9.
                    || (altitude >= 120. && f.plasma.ion_temperature_k < 250.)
            });
        if suspicious {
            let density = candidates
                .iter()
                .min_by(|a, b| a.0.total_cmp(&b.0))
                .map(|(_, _, f)| f.plasma.electron_density_m3)
                .unwrap_or(initial.electron_density_m3);
            for (ti, te, velocity, collision) in [
                (350., 700., -250., 30.),
                (700., 1400., 250., 3000.),
                (1500., 3000., 0., 100.),
                (3000., 6000., 0., 10000.),
            ] {
                let seed = PlasmaState {
                    electron_density_m3: density,
                    ion_temperature_k: ti,
                    electron_temperature_k: te,
                    velocity_m_s: velocity,
                    ion_collision_s: collision,
                };
                if let Ok(f) = problem.solve(seed) {
                    if f.chi_square.is_finite() {
                        candidates.push((f.chi_square, problem.clone(), f));
                    }
                }
            }
        }
    }
    if let Some((_, problem, mut best)) = candidates.into_iter().min_by(|a, b| a.0.total_cmp(&b.0))
    {
        gate.settings = problem.settings.clone();
        problem.temporal_uncertainty(&mut best, &samples)?;
        problem.model_uncertainty(&mut best)?;
        if best.status == "converged" && gate.first_lag_signal_to_noise.is_none_or(|s| s < 3.) {
            best.status = "weak_signal".into();
        }
        if best.status == "converged" && best.reduced_chi_square > 9. {
            best.status = "model_mismatch".into();
        }
        gate.status = best.status.clone();
        gate.result = Some(best);
    } else {
        gate.status = "numerical_failure".into();
    }
    Ok(gate)
}
pub fn run(path: &Path) -> Result<()> {
    let cfg: Config = serde_json::from_slice(&fs::read(path)?)?;
    ensure!(
        cfg.atmosphere.is_some(),
        "atmosphere configuration is required: prescribed Rust MSIS collisions and IRI composition"
    );
    ensure!(
        cfg.fit_collisions_below_km == 0. && cfg.model_selection_below_km.is_none(),
        "free collision fitting and model selection are disabled for prescribed atmosphere runs"
    );
    if cfg.streaming {
        return stream::run(&cfg);
    }
    ensure!(
        cfg.integration_s >= 15.
            && cfg.integration_s.is_finite()
            && cfg.composition_scale_km > 0.
            && cfg.antenna_efficiency > 0.
            && cfg.calibration_multiplier > 0.,
        "invalid fit configuration"
    );
    ensure!(
        cfg.altitude_edges_km.len() > 1
            && cfg.altitude_edges_km.iter().all(|h| h.is_finite())
            && cfg.altitude_edges_km.windows(2).all(|h| h[1] > h[0]),
        "invalid altitude bins"
    );
    let exp = Experiment::read(&cfg.experiment_model)?;
    let mut frames = vec![];
    for entry in fs::read_dir(&cfg.quicklook_frames)? {
        let path = entry?.path();
        if path.extension().is_some_and(|e| e == "h5") {
            frames.push(super::read_meta(&hdf5::File::open(path)?)?);
        }
    }
    frames.sort_by(|a, b| a.start_unix.total_cmp(&b.start_unix));
    ensure!(!frames.is_empty(), "no input frames");
    ensure!(
        frames
            .windows(2)
            .all(|w| w[1].start_unix >= w[0].end_unix - 0.01),
        "overlapping input integrations"
    );
    let mut windows: BTreeMap<i64, Vec<FrameMeta>> = BTreeMap::new();
    for m in frames {
        let middle = (m.start_unix + m.end_unix) / 2.;
        windows
            .entry((middle / cfg.integration_s).floor() as i64)
            .or_default()
            .push(m);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cfg.threads)
        .build()?;
    let mut profiles = vec![];
    for (window, metas) in windows {
        ensure!(
            metas.len() >= 3,
            "need at least three independent records per window"
        );
        ensure!(
            metas
                .iter()
                .all(|m| (m.elevation_deg - metas[0].elevation_deg).abs() < 0.1
                    && (m.azimuth_deg - metas[0].azimuth_deg).abs() < 0.1),
            "a fit window mixes beam directions"
        );
        let records = pool.install(|| {
            metas
                .par_iter()
                .map(|m| read_record(m, &exp, &cfg))
                .collect::<Result<Vec<_>>>()
        })?;
        let gates = pool.install(|| {
            cfg.altitude_edges_km
                .par_windows(2)
                .map(|h| fit_gate(&records, &exp, &cfg, h[0], h[1], &[]))
                .collect::<Result<Vec<_>>>()
        })?;
        let good = gates.iter().filter(|g| g.status == "converged").count();
        eprintln!(
            "UTC {}: {good}/{} converged, {} source integrations",
            window as f64 * cfg.integration_s,
            gates.len(),
            records.len()
        );
        profiles.push(Profile {
            start_unix: records
                .iter()
                .map(|r| r.meta.start_unix)
                .fold(f64::INFINITY, f64::min),
            end_unix: records
                .iter()
                .map(|r| r.meta.end_unix)
                .fold(f64::NEG_INFINITY, f64::max),
            source_count: records.len(),
            sources: records.iter().map(|r| r.meta.clone()).collect(),
            gates,
            processing_state: "complete".into(),
            detail_url: None,
            raw_ne: vec![],
            notes: vec![],
        });
    }
    let output=Output{schema:"beata22-ionfit-v1",config:&cfg,experiment_sha256:&exp.hash,uncertainty:"Temporal cluster sandwich over five-second records with conditional variance floor; approximate, includes variability; no absolute-calibration or composition systematics",limitations:vec!["Local homogeneous Maxwellian model with fixed altitude-dependent O+/30.5 amu composition; composition is assumed, not measured", "Number-conserving BGK collision model at low altitude; not a full magnetized ion-neutral collision model", "Only decoded signal lags 10..410 us are fitted; broad zero lag is retained separately and requires a joint range model", "Lag slices repeat the first lag and are not counted as independent measurements", "Residual score uses diagonal weighting; correlated lags mean it is not a calibrated chi-square goodness-of-fit test", "Velocity sign refers to positive model spectral shift; physical LOS sign unverified", "Density retains the event's assumed radar calibration (efficiency and calibration multiplier in config)"],profiles:&profiles};
    fs::create_dir_all(&cfg.output)?;
    fs::create_dir_all(&cfg.export)?;
    let json = serde_json::to_vec(&output)?;
    let h5path = cfg.output.join("ion-fits.h5");
    let tmp = tempfile::NamedTempFile::new_in(&cfg.output)?;
    {
        let h = hdf5::File::create(tmp.path())?;
        crate::storage::bytes(&h, "metadata_and_diagnostics_json", &json)?;
        h.new_dataset_builder()
            .with_data(&exp.frequency_hz)
            .create("model_frequency_hz")?;
        for (i, p) in profiles.iter().enumerate() {
            let gr = h.create_group(&format!("profile_{i:03}"))?;
            gr.new_dataset_builder()
                .with_data(&[p.start_unix, p.end_unix])
                .create("time_bounds_unix_s")?;
            gr.new_dataset_builder()
                .with_data(&p.gates.iter().map(|g| g.altitude_km).collect::<Vec<_>>())
                .create("altitude_km")?;
            let values: Vec<_> = p
                .gates
                .iter()
                .flat_map(|g| {
                    g.result
                        .as_ref()
                        .map(|r| {
                            [
                                r.plasma.electron_density_m3,
                                r.plasma.ion_temperature_k,
                                r.plasma.electron_temperature_k,
                                r.plasma.velocity_m_s,
                                r.plasma.ion_collision_s,
                            ]
                        })
                        .unwrap_or([f64::NAN; 5])
                })
                .collect();
            let errors: Vec<_> = p
                .gates
                .iter()
                .flat_map(|g| {
                    g.result
                        .as_ref()
                        .map(|r| r.standard_error.map(|e| e.unwrap_or(f64::NAN)))
                        .unwrap_or([f64::NAN; 5])
                })
                .collect();
            for (name, data) in [
                ("parameters_ne_ti_te_v_nu", values),
                ("standard_error", errors),
            ] {
                let d = crate::storage::array(&gr, name, [p.gates.len(), 5], &data)?;
                d.new_attr_builder()
                    .with_data("m^-3,K,K,m s^-1,s^-1")
                    .create("column_units")?;
            }
            for (j, g) in p.gates.iter().enumerate() {
                let gg = gr.create_group(&format!("gate_{j:03}"))?;
                gg.new_dataset_builder()
                    .with_data(&g.lag_us)
                    .create("lag_us")?;
                gg.new_dataset_builder()
                    .with_data(&[g.lower_km, g.upper_km])
                    .create("altitude_bin_km")?;
                gg.new_dataset_builder()
                    .with_data(g.status.as_bytes())
                    .create("status_utf8")?;
                let mut arrays = vec![("observed", g.observed.clone()), ("sigma", g.sigma.clone())];
                if let Some(r) = &g.result {
                    arrays.push(("predicted", r.predicted.clone()));
                    arrays.push(("standardized_residual", r.standardized_residual.clone()));
                    let curve: Vec<f64> = r
                        .model_curve
                        .iter()
                        .flat_map(|p| {
                            let e = p.standard_error.unwrap_or([f64::NAN; 2]);
                            [p.lag_s, p.value[0], p.value[1], e[0], e[1]]
                        })
                        .collect();
                    let ds = crate::storage::array(
                        &gg,
                        "model_curve",
                        [r.model_curve.len(), 5],
                        &curve,
                    )?;
                    ds.new_attr_builder().with_data("lag_s,mean_real,mean_imag,sigma_real,sigma_imag; s,m^-3,m^-3,m^-3,m^-3").create("columns")?;
                    let n = r.free_parameter_indices.len();
                    if r.covariance.len() == n * n {
                        let ds = gg
                            .new_dataset::<f64>()
                            .shape([n, n])
                            .create("parameter_covariance")?;
                        ds.write_raw(&r.covariance)?;
                    }
                    gg.new_dataset_builder()
                        .with_data(
                            &r.free_parameter_indices
                                .iter()
                                .map(|i| *i as u32)
                                .collect::<Vec<_>>(),
                        )
                        .create("free_parameter_indices")?;
                }
                for (name, values) in arrays {
                    if values.is_empty() {
                        continue;
                    }
                    let flat: Vec<_> = values.iter().flatten().copied().collect();
                    let ds = crate::storage::array(&gg, name, [values.len(), 2], &flat)?;
                    ds.new_attr_builder()
                        .with_data(if name == "standardized_residual" {
                            "1"
                        } else {
                            "m^-3; real,imag; response-normalized scattering strength"
                        })
                        .create("units")?;
                }
                let mut offsets = vec![0u32];
                let mut addresses = vec![];
                for row in &g.addresses {
                    addresses.extend(row.iter().map(|a| *a as u32));
                    offsets.push(addresses.len() as u32);
                }
                gg.new_dataset_builder()
                    .with_data(&offsets)
                    .create("source_address_offsets")?;
                gg.new_dataset_builder()
                    .with_data(&addresses)
                    .create("source_addresses_zero_based")?;
            }
        }
        h.flush()?;
    }
    tmp.persist(&h5path).map_err(|e| e.error)?;
    let dest = cfg.export.join("ion-fits.json");
    let mut tmp = tempfile::NamedTempFile::new_in(&cfg.export)?;
    use std::io::Write;
    let mut public: serde_json::Value = serde_json::from_slice(&json)?;
    for key in ["quicklook_frames", "experiment_model", "output", "export"] {
        public["config"].as_object_mut().unwrap().remove(key);
    }
    for profile in public["profiles"].as_array_mut().unwrap() {
        for source in profile["sources"].as_array_mut().unwrap() {
            source.as_object_mut().unwrap().remove("source");
        }
    }
    tmp.write_all(&serde_json::to_vec(&public)?)?;
    tmp.persist(dest).map_err(|e| e.error)?;
    eprintln!("Wrote {} profiles to {}", profiles.len(), h5path.display());
    Ok(())
}

pub fn compact_products(path: &Path) -> Result<()> {
    stream::compact_products(path)
}
