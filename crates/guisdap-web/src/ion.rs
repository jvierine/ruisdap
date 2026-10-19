//! Ion-fit inspection. All interaction and plotting run in Rust/WASM.
mod colors;
mod instruments;
mod measured;
mod raster;
use super::{el, error, json, text};
use serde::Deserialize;
use std::collections::{BTreeMap, BTreeSet};
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{prelude::*, JsCast};
use wasm_bindgen_futures::spawn_local;
use web_sys::{HtmlInputElement, HtmlSelectElement};
#[derive(Clone, Deserialize)]
struct Data {
    #[serde(default)]
    calibration: Option<CalibrationMetadata>,
    schema: String,
    uncertainty: String,
    limitations: Vec<String>,
    profiles: Vec<Profile>,
    #[serde(default)]
    processing: Option<Processing>,
}
#[derive(Clone, Deserialize)]
struct CalibrationMetadata {
    magic_constant: f64,
    antenna_efficiency: f64,
}
#[derive(Clone, Deserialize)]
struct Profile {
    start_unix: f64,
    end_unix: f64,
    source_count: usize,
    gates: Vec<Gate>,
    #[serde(default)]
    detail_url: Option<String>,
    #[serde(default)]
    processing_state: String,
    #[serde(default)]
    raw_ne: Vec<Option<[f64; 2]>>,
    #[serde(default)]
    notes: Vec<String>,
}
#[derive(Clone, Deserialize)]
struct Gate {
    #[serde(default)]
    prescribed_atmosphere: Option<PrescribedAtmosphere>,
    first_lag_signal_to_noise: Option<f64>,
    settings: Settings,
    altitude_km: f64,
    lower_km: f64,
    upper_km: f64,
    response_width_max_km: f64,
    oxygen_fraction: f64,
    status: String,
    result: Option<Fit>,
    #[serde(default)]
    lag_us: Vec<f64>,
    #[serde(default)]
    observed: Vec<[f64; 2]>,
    #[serde(default)]
    sigma: Vec<[f64; 2]>,
}
#[derive(Clone, Deserialize)]
struct PrescribedAtmosphere {
    molecular_fraction: f64,
    below_iri_composition_domain: bool,
}
#[derive(Clone, Deserialize)]
struct Settings {
    tie_electron_ion_temperature: bool,
    fit_collision_rate: bool,
    #[serde(default)]
    fixed_collision_s: f64,
}
#[derive(Clone, Deserialize)]
struct Fit {
    #[serde(default)]
    predicted: Vec<[f64; 2]>,
    #[serde(default)]
    covariance: Vec<f64>,
    #[serde(default)]
    free_parameter_indices: Vec<usize>,
    plasma: Plasma,
    standard_error: [Option<f64>; 5],
    reduced_chi_square: f64,
    iterations: usize,
    #[serde(default)]
    model_curve: Vec<ModelPoint>,
    #[serde(default)]
    standardized_residual: Vec<[f64; 2]>,
}
#[derive(Clone, Deserialize)]
struct ModelPoint {
    lag_s: f64,
    value: [f64; 2],
    standard_error: Option<[f64; 2]>,
}
#[derive(Clone, Deserialize)]
struct Plasma {
    electron_density_m3: f64,
    ion_temperature_k: f64,
    electron_temperature_k: f64,
    velocity_m_s: f64,
    ion_collision_s: f64,
}
#[derive(Clone, Deserialize)]
struct Processing {
    phase: String,
    threads: usize,
    total: usize,
    completed: usize,
    #[serde(default)]
    failed: usize,
    #[serde(default)]
    insufficient: usize,
    #[serde(default)]
    source_records: usize,
    #[serde(default)]
    updated_unix: i64,
}
#[derive(Clone, Deserialize)]
struct Dataset {
    label: String,
    integration_s: f64,
    data_url: String,
}
#[derive(Deserialize)]
struct Location {
    #[serde(default)]
    instruments_url: Option<String>,
    datasets: Vec<Dataset>,
    #[serde(default)]
    study_status_url: Option<String>,
}
type MapKey = (u64, u64, usize, usize, usize, i64);
struct Viewer {
    instruments: Option<instruments::Data>,
    data: Data,
    time: usize,
    gate: usize,
    lower: usize,
    upper: usize,
    window_start: f64,
    window_end: f64,
    datasets: Vec<Dataset>,
    resolution: usize,
    details: BTreeMap<usize, Profile>,
    loading: BTreeSet<usize>,
    map_key: std::cell::Cell<Option<MapKey>>,
    statuses: Vec<Option<Processing>>,
    refresh_busy: bool,
    drag: Option<(f64, f64, f64)>,
    scales: [[f64; 2]; 5],
    quality_limit: f64,
    quality: guisdap_core::quality::Limits,
    mask_bad: bool,
    hatch_bad: bool,
}
fn rejection_reasons(g: &Gate, v: &Viewer) -> Vec<&'static str> {
    let Some(f) = &g.result else {
        return vec!["fit unavailable"];
    };
    guisdap_core::quality::reasons(
        &g.status,
        [
            f.plasma.electron_density_m3,
            f.plasma.ion_temperature_k,
            f.plasma.electron_temperature_k,
            f.plasma.velocity_m_s,
        ],
        [
            f.standard_error[0],
            f.standard_error[1],
            f.standard_error[2],
            f.standard_error[3],
        ],
        g.first_lag_signal_to_noise,
        f.reduced_chi_square,
        v.quality,
    )
}
fn rejected(g: &Gate, v: &Viewer) -> bool {
    !rejection_reasons(g, v).is_empty()
}
fn selected(v: &Viewer) -> &Profile {
    v.details.get(&v.time).unwrap_or(&v.data.profiles[v.time])
}

fn clock(t: f64) -> String {
    let s = (t.round() as i64).rem_euclid(86400);
    format!("{:02}:{:02}:{:02}", s / 3600, s / 60 % 60, s % 60)
}
fn value(p: &Plasma, i: usize) -> f64 {
    match i {
        0 => p.electron_density_m3,
        1 => p.ion_temperature_k,
        2 => p.electron_temperature_k,
        3 => p.velocity_m_s,
        _ => p.ion_collision_s,
    }
}
fn profile(v: &Viewer, param: usize) -> String {
    let p = selected(v);
    let [lo, hi] = v.scales[param];
    let log = param == 0;
    let tr = |x: f64| {
        if log {
            if x > 0. {
                x.log10()
            } else {
                f64::NAN
            }
        } else {
            x
        }
    };
    let mut good = Vec::new();
    let mut bad = Vec::new();
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut selected = Vec::new();
    let mut errors = Vec::new();
    for (i, g) in p.gates.iter().enumerate().take(v.upper + 1).skip(v.lower) {
        if v.mask_bad && rejected(g, v) {
            continue;
        }
        let Some(f) = &g.result else {
            continue;
        };
        let z = value(&f.plasma, param);
        if !z.is_finite() || (log && z <= 0.) {
            continue;
        }
        let pos = tr(z.clamp(lo, hi));
        let target = if z < lo {
            &mut left
        } else if z > hi {
            &mut right
        } else if g.status == "converged" {
            &mut good
        } else {
            &mut bad
        };
        target.extend([pos, g.altitude_km]);
        if let Some(e) = f.standard_error[param] {
            errors.extend([
                tr((z - e).max(lo)),
                g.altitude_km,
                tr((z + e).min(hi)),
                g.altitude_km,
            ]);
        }
        if i == v.gate {
            selected.extend([pos, g.altitude_km]);
        }
    }
    serde_json::json!({"title":"Selected-time plasma parameter profile","height":350,"bounds":[tr(lo),tr(hi),p.gates[v.lower].lower_km,p.gates[v.upper].upper_km],"xlabel":(["Ne / m⁻³","Ti / K","Te / K","Velocity / m s⁻¹"][param]),"xformat":if log{"log10"}else{""},"ylabel":"Altitude / km","series":[{"xy":errors,"segments":true,"color":"#7b929e"},{"xy":good,"style":0,"radius":3},{"xy":bad,"style":1,"radius":3},{"xy":left,"style":2,"radius":3},{"xy":right,"style":3,"radius":3},{"xy":selected,"style":1,"radius":6,"color":"#986327"}]}).to_string()
}

fn time_altitude(v: &Viewer, param: usize) -> Result<(), JsValue> {
    let p = &v.data.profiles[v.time];
    let bottom = p.gates[v.lower].lower_km;
    let top = p.gates[v.upper].upper_km;
    let start = v.window_start;
    let end = v.window_end;
    let [min, max] = v.scales[param];
    text(
        &format!("fit-scale-applied-{param}"),
        &format!("Applied limits: {min:.3e} to {max:.3e}"),
    );
    let nx = 1200usize;
    let ny = 240usize;
    let mut pixels = vec![f32::NAN; nx * ny];
    let mut flags = vec![0u8; nx * ny];
    for ix in 0..nx {
        let time = start + (end - start) * (ix as f64 + 0.5) / nx as f64;
        let ti = v.data.profiles.partition_point(|p| p.end_unix <= time);
        let Some(p) = v.data.profiles.get(ti).filter(|p| p.start_unix <= time) else {
            continue;
        };
        for iy in 0..ny {
            let h = bottom + (top - bottom) * (iy as f64 + 0.5) / ny as f64;
            let gi = p.gates.partition_point(|g| g.upper_km <= h);
            let Some(g) = p.gates.get(gi).filter(|g| g.lower_km <= h) else {
                continue;
            };
            let (raw, flag) = if param == 4 {
                let Some(Some(raw)) = p.raw_ne.get(gi) else {
                    continue;
                };
                (raw[0], false)
            } else {
                let Some(f) = &g.result else { continue };
                if g.status == "singular" {
                    continue;
                }
                let bad = rejected(g, v);
                if bad && v.mask_bad {
                    continue;
                }
                (value(&f.plasma, param), bad && v.hatch_bad)
            };
            if !raw.is_finite() || ((param == 0 || param == 4) && raw <= 0.) {
                continue;
            }
            pixels[iy * nx + ix] = raw as f32;
            flags[iy * nx + ix] = u8::from(flag);
        }
    }
    let spec = serde_json::json!({"title":(["Ne / m⁻³","Ti / K","Te / K","Velocity / m s⁻¹","Raw Ne / m⁻³"][param]),"bounds":[start,end,bottom,top],"limits":v.scales[param],"mode":if param==0||param==4{1}else{0},"normalized":false,"palette":if param==3{0}else{1},"external_navigation":true,"ylabel":"Altitude / km"});
    raster::draw(
        &format!("fit-raster-{param}"),
        &pixels,
        &flags,
        nx,
        ny,
        &spec,
    )?;
    Ok(())
}

fn cursor(v: &Viewer, param: usize) {
    let p = &v.data.profiles[v.time];
    let g = &p.gates[v.gate];
    let x = ((p.start_unix + p.end_unix) / 2. - v.window_start) / (v.window_end - v.window_start);
    let y = (p.gates[v.upper].upper_km - g.altitude_km)
        / (p.gates[v.upper].upper_km - p.gates[v.lower].lower_km);
    el(&format!("fit-cursor-{param}")).set_inner_html(&if (0. ..=1.).contains(&x){format!("<div class='fit-selection-cursor' style='left:calc(76px + (100% - 208px)*{x});top:22px;bottom:58px'><span style='top:{}%'></span></div>",y*100.)}else{String::new()});
}

fn acf(v: &Viewer, component: usize) -> String {
    let g = &selected(v).gates[v.gate];
    let mut magnitude = g
        .observed
        .iter()
        .zip(&g.sigma)
        .map(|(a, e)| a[component].abs() + e[component])
        .fold(1., f64::max);
    if let Some(f) = &g.result {
        for p in &f.model_curve {
            magnitude = magnitude
                .max(p.value[component].abs() + p.standard_error.map_or(0., |e| e[component]));
        }
    }
    let scale = 10_f64.powf(magnitude.log10().floor());
    let ymax = magnitude / scale * 1.1;
    let observed: Vec<_> = g
        .lag_us
        .iter()
        .zip(&g.observed)
        .flat_map(|(t, a)| [*t, a[component] / scale])
        .collect();
    let errors: Vec<_> = g
        .lag_us
        .iter()
        .zip(&g.observed)
        .zip(&g.sigma)
        .flat_map(|((t, a), e)| {
            [
                *t,
                (a[component] - e[component]) / scale,
                *t,
                (a[component] + e[component]) / scale,
            ]
        })
        .collect();
    let model: Vec<f64> = g
        .result
        .iter()
        .flat_map(|f| &f.model_curve)
        .flat_map(|p| [p.lag_s * 1e6, p.value[component] / scale])
        .collect();
    let band: Vec<f64> = g
        .result
        .iter()
        .flat_map(|f| &f.model_curve)
        .flat_map(|p| {
            [
                p.lag_s * 1e6,
                p.standard_error
                    .map_or(f64::NAN, |e| (p.value[component] - e[component]) / scale),
                p.standard_error
                    .map_or(f64::NAN, |e| (p.value[component] + e[component]) / scale),
            ]
        })
        .collect();
    serde_json::json!({"title":if component==0{"Real ACF"}else{"Imaginary ACF"},"bounds":[0.,420.,-ymax,ymax],"xlabel":"Measured lag / µs","ylabel":format!("ACF / {scale:.0e} m⁻³"),"series":[{"xy":errors,"segments":true,"color":"#7b929e"},{"xy":observed,"line":false,"radius":2.4,"color":"#647b89"},{"xy":model,"color":"#236b8e"}],"bands":[{"values":band,"color":"#236b8e4d"}]}).to_string()
}

fn render(v: &Viewer) {
    colors::sync(v);
    text("fit-uncertainty", &v.data.uncertainty);
    text("fit-assumptions", &v.data.limitations.join(". "));
    text("fit-calibration-used",&v.data.calibration.as_ref().map(|c|format!("Magic constant used: {:.4} · antenna efficiency: {:.3} · calibration recorded in this analysis product.",c.magic_constant,c.antenna_efficiency)).unwrap_or_else(||"Magic constant used: not recorded in this product.".into()));
    instruments::render(v);
    let p = selected(v);
    let g = &p.gates[v.gate];
    let altitude_slider: HtmlInputElement = el("fit-altitude").dyn_into().unwrap();
    altitude_slider.set_min(&v.lower.to_string());
    altitude_slider.set_max(&v.upper.to_string());
    altitude_slider.set_value(&v.gate.to_string());
    altitude_slider
        .set_attribute(
            "aria-valuetext",
            &format!("{:.0}–{:.0} km", g.lower_km, g.upper_km),
        )
        .unwrap();
    let selector: HtmlSelectElement = el("fit-gate").dyn_into().unwrap();
    selector.set_value(&v.gate.to_string());
    text(
        "fit-altitude-label",
        &format!("Selected altitude: {:.0}–{:.0} km", g.lower_km, g.upper_km),
    );
    text(
        "fit-range-label",
        &format!(
            "Displayed altitude: {:.0}–{:.0} km",
            p.gates[v.lower].lower_km, p.gates[v.upper].upper_km
        ),
    );
    for (id, index, altitude) in [
        ("fit-altitude-min", v.lower, p.gates[v.lower].lower_km),
        ("fit-altitude-max", v.upper, p.gates[v.upper].upper_km),
    ] {
        let slider: HtmlInputElement = el(id).dyn_into().unwrap();
        slider.set_value(&index.to_string());
        slider
            .set_attribute("aria-valuetext", &format!("{altitude:.0} km"))
            .unwrap();
    }
    text(
        "fit-time",
        &format!(
            "{}–{} UTC · {} records",
            clock(p.start_unix),
            clock(p.end_unix),
            p.source_count
        ),
    );
    text(
        "fit-status",
        &format!(
            "{} · {} / {} gates converged",
            p.processing_state,
            p.gates.iter().filter(|g| g.status == "converged").count(),
            p.gates.len()
        ),
    );
    for i in 0..4 {
        crate::plotting::shared_curve(&format!("fit-profile-{i}"), &profile(v, i));
    }
    text("fit-quality-value", &format!("{:.1}", v.quality_limit));
    let rejected_count = p.gates.iter().filter(|g| rejected(g, v)).count();
    text("fit-quality-summary", &format!("Selected interval: {rejected_count}/{} gates rejected or unavailable. Screened using convergence, bounds, uncertainty and signal significance.", p.gates.len()));
    let mut counts = BTreeMap::new();
    for g in &p.gates {
        for reason in rejection_reasons(g, v) {
            *counts.entry(reason).or_insert(0) += 1;
        }
    }
    text(
        "fit-quality-reasons",
        &counts
            .iter()
            .map(|(r, n)| format!("{r}: {n}"))
            .collect::<Vec<_>>()
            .join(" · "),
    );
    for (id, value) in [
        (
            "fit-ne-error-value",
            format!("{:.0}%", v.quality.density_relative_error * 100.),
        ),
        (
            "fit-temperature-error-value",
            format!("{:.0}%", v.quality.temperature_relative_error * 100.),
        ),
        (
            "fit-velocity-error-value",
            format!("{:.0} m/s", v.quality.velocity_error_m_s),
        ),
        (
            "fit-signal-value",
            format!("{:.1}", v.quality.minimum_signal_significance),
        ),
    ] {
        text(id, &value);
    }
    let revision = v.data.processing.as_ref().map_or(0, |p| p.updated_unix);
    let key = (
        v.window_start.to_bits(),
        v.window_end.to_bits(),
        v.lower,
        v.upper,
        v.resolution,
        revision,
    );
    if v.map_key.get() != Some(key) {
        for i in 0..5 {
            if let Err(e) = time_altitude(v, i) {
                text("fit-status", &format!("Plot error: {e:?}"));
            }
        }
        v.map_key.set(Some(key));
    }
    for i in 0..5 {
        cursor(v, i);
    }
    text(
        "fit-window",
        &format!("{}–{} UTC", clock(v.window_start), clock(v.window_end)),
    );
    let first = v
        .data
        .profiles
        .partition_point(|p| p.end_unix <= v.window_start)
        .min(v.data.profiles.len() - 1);
    let last = v
        .data
        .profiles
        .partition_point(|p| p.start_unix < v.window_end)
        .saturating_sub(1);
    let time_slider: HtmlInputElement = el("fit-slider").dyn_into().unwrap();
    time_slider.set_min(&first.to_string());
    time_slider.set_max(&last.to_string());
    time_slider.set_value(&v.time.to_string());
    for (id, t) in [
        ("fit-window-start", v.window_start),
        ("fit-window-end", v.window_end),
    ] {
        el(id)
            .dyn_into::<HtmlInputElement>()
            .unwrap()
            .set_value(&clock(t));
    }
    for i in 0..2 {
        crate::plotting::shared_curve(&format!("fit-acf-{i}"), &acf(v, i));
    }
    crate::plotting::shared_curve("fit-spectrum", &analytical_spectrum(v));
    crate::plotting::shared_curve("fit-measured-spectrum", &measured::spectrum(v));
    text(
        "fit-gate-model",
        &format!(
            "Te/Ti: {}. Collision rate: {}. First-lag signal/noise: {:.1}.",
            if g.settings.tie_electron_ion_temperature {
                "fixed to 1"
            } else {
                "fitted"
            },
            if g.settings.fit_collision_rate {
                "fitted effective BGK rate".to_string()
            } else {
                if let Some(atmosphere) = &g.prescribed_atmosphere {
                    format!(
                        "MSIS {:.3e} s⁻¹ (fixed); IRI DS-1995 molecular fraction {:.3}{}",
                        g.settings.fixed_collision_s,
                        atmosphere.molecular_fraction,
                        if atmosphere.below_iri_composition_domain {
                            "; below IRI domain, 75 km boundary composition"
                        } else {
                            ""
                        }
                    )
                } else {
                    format!("fixed {:.2} s⁻¹", g.settings.fixed_collision_s)
                }
            },
            g.first_lag_signal_to_noise.unwrap_or(f64::NAN)
        ),
    );
    let detail = if let Some(f) = &g.result {
        let rms = (f
            .standardized_residual
            .iter()
            .flat_map(|a| a.iter())
            .map(|x| x * x)
            .sum::<f64>()
            / (2 * f.standardized_residual.len()).max(1) as f64)
            .sqrt();
        format!("{} · Ne {:.3e} m⁻³ · Ti {:.0} K · Te {:.0} K · v {:.1} m/s · ν {:.1} s⁻¹ · weighted residual score {:.2} · residual RMS {:.2} · {} iterations",g.status,f.plasma.electron_density_m3,f.plasma.ion_temperature_k,f.plasma.electron_temperature_k,f.plasma.velocity_m_s,f.plasma.ion_collision_s,f.reduced_chi_square,rms,f.iterations)
    } else {
        g.status.clone()
    };
    let detail = if !v.details.contains_key(&v.time) && p.detail_url.is_some() {
        format!("Loading measured/model lags… {}", p.notes.join("; "))
    } else {
        format!("{detail} {}", p.notes.join("; "))
    };
    text("fit-detail", &detail);
    text("fit-support",&format!("Selected altitude bin {:.0}–{:.0} km; maximum equivalent response width {:.2} km (slant). Assumed O⁺ fraction {:.3}. {} response groups; all available decoded ion lags retained.",g.lower_km,g.upper_km,g.response_width_max_km,g.oxygen_fraction,g.lag_us.len()));
}
fn status_view(v: &Viewer) {
    let rows=v.datasets.iter().zip(&v.statuses).map(|(d,p)|{
        if let Some(p)=p {format!("{}: {} — {}/{} intervals; {} failed, {} insufficient-record windows; {} source records; {} worker threads; updated {} UTC",d.label,p.phase,p.completed,p.total,p.failed,p.insufficient,p.source_records,p.threads,clock(p.updated_unix as f64))}
        else{format!("{}: queued / status not yet available",d.label)}
    }).collect::<Vec<_>>().join("\n\n");
    text("processing-fits", &rows);
    if let Some(p) = v.statuses.get(v.resolution).and_then(Option::as_ref) {
        text(
            "fit-progress",
            &format!(
                "{} · {}/{} intervals · {} threads",
                p.phase, p.completed, p.total, p.threads
            ),
        );
    }
}
fn fit_index(data: &Data, time: f64) -> usize {
    data.profiles
        .partition_point(|p| p.end_unix <= time)
        .min(data.profiles.len() - 1)
}
fn set_window(v: &mut Viewer, start: f64, end: f64) {
    let full_start = v.data.profiles[0].start_unix;
    let full_end = v.data.profiles.last().unwrap().end_unix;
    let width = (end - start)
        .max(v.datasets[v.resolution].integration_s)
        .min(full_end - full_start);
    v.window_start = start.clamp(full_start, full_end - width);
    v.window_end = v.window_start + width;
    let selected = (v.data.profiles[v.time].start_unix + v.data.profiles[v.time].end_unix) / 2.;
    v.time = fit_index(
        &v.data,
        selected.clamp(v.window_start + 0.001, v.window_end - 0.001),
    );
}
fn load_detail(state: Rc<RefCell<Viewer>>) {
    let request = {
        let mut v = state.borrow_mut();
        let index = v.time;
        if v.details.contains_key(&index) || v.loading.contains(&index) {
            return;
        }
        let Some(url) = v.data.profiles[index].detail_url.clone() else {
            return;
        };
        let base = v.datasets[v.resolution]
            .data_url
            .rsplit_once('/')
            .map_or("", |(base, _)| base)
            .to_string();
        let resolution = v.resolution;
        v.loading.insert(index);
        (index, resolution, format!("{base}/{url}"))
    };
    spawn_local(async move {
        let result = json::<Profile>(&request.2).await;
        let mut v = state.borrow_mut();
        if v.resolution != request.1 {
            return;
        }
        v.loading.remove(&request.0);
        match result {
            Ok(p) => {
                v.details.insert(request.0, p);
                while v.details.len() > 8 {
                    let furthest = *v
                        .details
                        .keys()
                        .max_by_key(|&&i| i.abs_diff(v.time))
                        .unwrap();
                    v.details.remove(&furthest);
                }
                if v.time == request.0 {
                    render(&v);
                }
            }
            Err(e) => {
                if v.time == request.0 {
                    text(
                        "fit-detail",
                        &format!("Detail load failed; retry by selecting interval: {e:?}"),
                    );
                }
            }
        }
    });
}
async fn change_resolution(state: Rc<RefCell<Viewer>>, resolution: usize) -> Result<(), JsValue> {
    let (url, target) = {
        let v = state.borrow();
        (
            v.datasets[resolution].data_url.clone(),
            v.data.profiles[v.time].start_unix + 0.1,
        )
    };
    el("fit-resolution").set_attribute("disabled", "")?;
    text("fit-progress", "Loading resolution…");
    let loaded = json::<Data>(&url).await;
    el("fit-resolution").remove_attribute("disabled")?;
    let data = loaded?;
    if data.profiles.is_empty() {
        return Err(error("empty resolution"));
    }
    {
        let mut v = state.borrow_mut();
        let start = v.window_start;
        let end = v.window_end;
        let old = &v.data.profiles[v.time].gates;
        let low = old[v.lower].lower_km;
        let high = old[v.upper].upper_km;
        let height = old[v.gate].altitude_km;
        v.data = data;
        let gates = &v.data.profiles[0].gates;
        let lower = gates
            .partition_point(|g| g.upper_km <= low)
            .min(gates.len() - 1);
        let upper = gates
            .partition_point(|g| g.lower_km < high)
            .saturating_sub(1)
            .max(lower);
        let gate = gates
            .partition_point(|g| g.upper_km <= height)
            .clamp(lower, upper);
        let max = (gates.len() - 1).to_string();
        el("fit-resolution")
            .dyn_into::<HtmlSelectElement>()?
            .set_value(&resolution.to_string());
        el("fit-gate").set_inner_html(
            &gates
                .iter()
                .enumerate()
                .map(|(i, g)| {
                    format!(
                        "<option value='{i}'>{:.1}–{:.1} km</option>",
                        g.lower_km, g.upper_km
                    )
                })
                .collect::<String>(),
        );
        v.lower = lower;
        v.upper = upper;
        v.gate = gate;
        for id in ["fit-altitude", "fit-altitude-min", "fit-altitude-max"] {
            el(id).dyn_into::<HtmlInputElement>()?.set_max(&max);
        }
        v.resolution = resolution;
        v.time = fit_index(&v.data, target);
        v.details.clear();
        v.loading.clear();
        v.map_key.set(None);
        v.statuses[resolution] = v.data.processing.clone();
        set_window(&mut v, start, end);
        render(&v);
        status_view(&v);
    }
    load_detail(state);
    Ok(())
}
async fn refresh(state: Rc<RefCell<Viewer>>, study: Option<String>) {
    {
        let mut v = state.borrow_mut();
        if v.refresh_busy {
            return;
        }
        v.refresh_busy = true;
    }
    let datasets = state.borrow().datasets.clone();
    let original = web_sys::window()
        .unwrap()
        .location()
        .search()
        .unwrap_or_default()
        .trim_start_matches('?')
        .split('&')
        .any(|p| p == "original=1");
    if !original {
        if let Ok(current) = json::<Location>(&super::experiment_config("ion-config.json")).await {
            let changed = current
                .datasets
                .iter()
                .map(|d| &d.data_url)
                .ne(datasets.iter().map(|d| &d.data_url));
            let _ = el("fit-product-update").set_attribute(
                "style",
                if changed {
                    "height:2em;visibility:visible"
                } else {
                    "height:2em;visibility:hidden"
                },
            );
        }
    }

    for (i, d) in datasets.iter().enumerate() {
        let status_url = format!(
            "{}/status.json",
            d.data_url.rsplit_once('/').map_or("", |(a, _)| a)
        );
        if let Ok(p) = json::<Processing>(&status_url).await {
            let needs = {
                let v = state.borrow();
                v.resolution == i
                    && v.data
                        .processing
                        .as_ref()
                        .is_none_or(|old| old.completed != p.completed || old.phase != p.phase)
            };
            state.borrow_mut().statuses[i] = Some(p);
            if needs {
                if let Ok(data) = json::<Data>(&d.data_url).await {
                    let mut v = state.borrow_mut();
                    if v.resolution == i && !data.profiles.is_empty() {
                        v.data = data;
                        v.map_key.set(None);
                        render(&v);
                    }
                }
            }
        }
    }
    if let Some(url) = study {
        if let Ok(p) = json::<serde_json::Value>(&url).await {
            text(
                "processing-study",
                &format!("Dataset {} · {} · {} configured worker threads · {} source records\nUpdated {} UTC\n\n{}",p["dataset"].as_str().unwrap_or(""),p["state"].as_str().unwrap_or(""),p["configured_worker_threads"],p["source_records"],clock(p["updated_unix"].as_f64().unwrap_or(0.)),p["stages"].as_array().map(|stages|stages.iter().map(|s|format!("{}: {} — {} / {}",s["name"].as_str().unwrap_or(""),s["state"].as_str().unwrap_or(""),s["completed"],s["total"])).collect::<Vec<_>>().join("\n")).unwrap_or_default()),
            );
        }
    }
    {
        let mut v = state.borrow_mut();
        status_view(&v);
        v.refresh_busy = false;
    }
    load_detail(state);
}
fn seconds(s: &str) -> Option<f64> {
    let parts: Vec<_> = s
        .split(':')
        .map(str::parse::<u32>)
        .collect::<Result<_, _>>()
        .ok()?;
    if parts.len() < 2 || parts[0] > 23 || parts[1] > 59 || parts.get(2).is_some_and(|s| *s > 59) {
        return None;
    }
    Some((parts[0] * 3600 + parts[1] * 60 + parts.get(2).copied().unwrap_or(0)) as f64)
}
#[wasm_bindgen]
pub async fn start_ion() -> Result<(), JsValue> {
    crate::navigation::setup()?;
    let original = web_sys::window()
        .unwrap()
        .location()
        .search()?
        .trim_start_matches('?')
        .split('&')
        .any(|p| p == "original=1");
    let config_url = super::experiment_config("ion-config.json");
    let location: Location = if original {
        json("ion-config-original.json").await?
    } else if config_url == "ion-config.json" {
        match json("catalog/default-ion-config.json").await {
            Ok(v) => v,
            Err(_) => json(&config_url).await?,
        }
    } else {
        json(&config_url).await?
    };
    if location.datasets.is_empty() {
        return Err(error("no configured resolutions"));
    }
    let instruments = if let Some(url) = &location.instruments_url {
        json::<instruments::Data>(url).await.ok()
    } else {
        None
    };
    let resolution = web_sys::window()
        .unwrap()
        .location()
        .search()?
        .trim_start_matches('?')
        .split('&')
        .find_map(|p| {
            p.strip_prefix("resolution=")
                .and_then(|v| v.parse::<usize>().ok())
        })
        .unwrap_or(0)
        .min(location.datasets.len() - 1);
    let data: Data = json(&location.datasets[resolution].data_url).await?;
    if data.schema != "beata22-ionfit-v2" || data.profiles.is_empty() {
        return Err(error("unsupported or empty ion product"));
    }
    el("fit-resolution").set_inner_html(
        &location
            .datasets
            .iter()
            .enumerate()
            .map(|(i, d)| format!("<option value='{i}'>{}</option>", d.label))
            .collect::<String>(),
    );
    el("fit-resolution")
        .dyn_into::<HtmlSelectElement>()?
        .set_value(&resolution.to_string());
    el("fit-gate").set_inner_html(
        &data.profiles[0]
            .gates
            .iter()
            .enumerate()
            .map(|(i, g)| {
                format!(
                    "<option value='{i}'>{:.0}–{:.0} km</option>",
                    g.lower_km, g.upper_km
                )
            })
            .collect::<String>(),
    );
    let time = fit_index(&data, 1667408700.1);
    let upper = data.profiles[0].gates.len() - 1;
    for id in ["fit-altitude", "fit-altitude-min", "fit-altitude-max"] {
        el(id)
            .dyn_into::<HtmlInputElement>()?
            .set_max(&upper.to_string());
    }
    let gate = data.profiles[0]
        .gates
        .iter()
        .position(|g| g.altitude_km >= 200.)
        .unwrap_or(0);
    let mut statuses = vec![None; location.datasets.len()];
    statuses[resolution] = data.processing.clone();
    let state = Rc::new(RefCell::new(Viewer {
        instruments,
        window_start: data.profiles[0].start_unix,
        window_end: data.profiles.last().unwrap().end_unix,
        data,
        time,
        gate,
        lower: 0,
        upper,
        datasets: location.datasets,
        resolution,
        details: BTreeMap::new(),
        loading: BTreeSet::new(),
        map_key: std::cell::Cell::new(None),
        statuses,
        refresh_busy: false,
        drag: None,
        scales: [
            [1e8, 1e13],
            [0., 2000.],
            [0., 4000.],
            [-1000., 1000.],
            [1e8, 1e13],
        ],
        quality_limit: 10.,
        quality: guisdap_core::quality::Limits::default(),
        mask_bad: false,
        hatch_bad: false,
    }));
    render(&state.borrow());
    status_view(&state.borrow());
    load_detail(state.clone());
    settings_drag()?;
    colors::install(state.clone())?;
    for (id, open) in [("fit-settings-open", true), ("fit-settings-close", false)] {
        let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            if open {
                el("fit-settings").remove_attribute("hidden").unwrap();
            } else {
                el("fit-settings").set_attribute("hidden", "").unwrap();
            }
            el("fit-settings-open")
                .set_attribute("aria-expanded", if open { "true" } else { "false" })
                .unwrap();
        });
        el(id).add_event_listener_with_callback("click", cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    for param in 0..5 {
        for edge in 0..2 {
            let id = format!("fit-scale-{param}-{edge}");
            let state = state.clone();
            let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
                let read = |k| {
                    el(&format!("fit-scale-{param}-{k}"))
                        .dyn_into::<HtmlInputElement>()
                        .unwrap()
                        .value()
                        .parse::<f64>()
                        .ok()
                };
                if let (Some(a), Some(b)) = (read(0), read(1)) {
                    if a.is_finite()
                        && b.is_finite()
                        && a < b
                        && (!(param == 0 || param == 4) || a > 0.)
                        && (param != 3 || (a < 0. && b > 0.))
                    {
                        let mut v = state.borrow_mut();
                        v.scales[param] = [a, b];
                        v.map_key.set(None);
                        render(&v);
                        text("fit-scale-error", "");
                        return;
                    }
                }
                text("fit-scale-error", "Colour limits need finite min < max; density limits must be positive; velocity limits must straddle zero.");
            });
            if edge == 0 {
                el(&format!("fit-scale-apply-{param}"))
                    .add_event_listener_with_callback("click", cb.as_ref().unchecked_ref())?;
            }
            el(&id).add_event_listener_with_callback("input", cb.as_ref().unchecked_ref())?;
            el(&id).add_event_listener_with_callback("change", cb.as_ref().unchecked_ref())?;
            cb.forget();
        }
    }
    for (id, event) in [
        ("fit-quality", "input"),
        ("fit-ne-error", "input"),
        ("fit-temperature-error", "input"),
        ("fit-velocity-error", "input"),
        ("fit-signal", "input"),
        ("fit-residual-enabled", "change"),
        ("fit-mask-bad", "change"),
        ("fit-hatch-bad", "change"),
        ("fit-slider", "input"),
        ("fit-gate", "change"),
        ("fit-altitude", "input"),
        ("fit-altitude-min", "input"),
        ("fit-altitude-max", "input"),
        ("fit-window-start", "change"),
        ("fit-window-end", "change"),
        ("fit-zoom-in", "click"),
        ("fit-zoom-out", "click"),
        ("fit-view-all", "click"),
        ("fit-pan-left", "click"),
        ("fit-pan-right", "click"),
    ] {
        let state = state.clone();
        let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            {
                let mut v = state.borrow_mut();
                let last = v.data.profiles[v.time].gates.len() - 1;
                match id {
                    "fit-quality"
                    | "fit-mask-bad"
                    | "fit-hatch-bad"
                    | "fit-ne-error"
                    | "fit-temperature-error"
                    | "fit-velocity-error"
                    | "fit-signal"
                    | "fit-residual-enabled" => {
                        v.quality_limit = 10_f64.powf(
                            el("fit-quality")
                                .dyn_into::<HtmlInputElement>()
                                .unwrap()
                                .value()
                                .parse::<f64>()
                                .unwrap_or(1.),
                        );
                        let number = |id: &str, fallback: f64| {
                            el(id)
                                .dyn_into::<HtmlInputElement>()
                                .unwrap()
                                .value()
                                .parse::<f64>()
                                .unwrap_or(fallback)
                        };
                        v.quality.density_relative_error = number("fit-ne-error", 50.) / 100.;
                        v.quality.temperature_relative_error =
                            number("fit-temperature-error", 100.) / 100.;
                        v.quality.velocity_error_m_s = number("fit-velocity-error", 300.);
                        v.quality.minimum_signal_significance = number("fit-signal", 3.);
                        v.quality.maximum_residual = el("fit-residual-enabled")
                            .dyn_into::<HtmlInputElement>()
                            .unwrap()
                            .checked()
                            .then_some(v.quality_limit);
                        v.mask_bad = el("fit-mask-bad")
                            .dyn_into::<HtmlInputElement>()
                            .unwrap()
                            .checked();
                        v.hatch_bad = el("fit-hatch-bad")
                            .dyn_into::<HtmlInputElement>()
                            .unwrap()
                            .checked();
                        v.map_key.set(None);
                    }

                    "fit-window-start" | "fit-window-end" => {
                        let day = (v.data.profiles[0].start_unix / 86400.).floor() * 86400.;
                        let a = el("fit-window-start")
                            .dyn_into::<HtmlInputElement>()
                            .unwrap()
                            .value();
                        let b = el("fit-window-end")
                            .dyn_into::<HtmlInputElement>()
                            .unwrap()
                            .value();
                        if let (Some(a), Some(b)) = (seconds(&a), seconds(&b)) {
                            if b > a {
                                set_window(&mut v, day + a, day + b);
                            }
                        }
                    }
                    "fit-view-all" => {
                        let a = v.data.profiles[0].start_unix;
                        let b = v.data.profiles.last().unwrap().end_unix;
                        set_window(&mut v, a, b);
                    }
                    "fit-zoom-in" | "fit-zoom-out" => {
                        let center = (v.data.profiles[v.time].start_unix
                            + v.data.profiles[v.time].end_unix)
                            / 2.;
                        let width = (v.window_end - v.window_start)
                            * if id == "fit-zoom-in" { 0.5 } else { 2. };
                        set_window(&mut v, center - width / 2., center + width / 2.);
                    }
                    "fit-pan-left" | "fit-pan-right" => {
                        let shift = (v.window_end - v.window_start)
                            * if id == "fit-pan-left" { -0.5 } else { 0.5 };
                        let a = v.window_start + shift;
                        let b = v.window_end + shift;
                        set_window(&mut v, a, b);
                    }
                    _ => {
                        let raw = if id == "fit-gate" {
                            el(id).dyn_into::<HtmlSelectElement>().unwrap().value()
                        } else {
                            el(id).dyn_into::<HtmlInputElement>().unwrap().value()
                        };
                        let index = raw.parse::<usize>().unwrap_or(0);
                        match id {
                            "fit-slider" => v.time = index.min(v.data.profiles.len() - 1),
                            "fit-altitude-min" => v.lower = index.min(v.upper),
                            "fit-altitude-max" => v.upper = index.max(v.lower).min(last),
                            _ => v.gate = index.min(last),
                        }
                    }
                }
                v.gate = v.gate.clamp(v.lower, v.upper);
                render(&v);
            }
            load_detail(state.clone());
        });
        el(id).add_event_listener_with_callback(event, cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    {
        let state = state.clone();
        let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            let index = el("fit-resolution")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value()
                .parse::<usize>()
                .unwrap_or(0);
            let state = state.clone();
            spawn_local(async move {
                if let Err(e) = change_resolution(state, index).await {
                    text("fit-progress", &format!("Resolution unavailable: {e:?}"));
                }
            });
        });
        el("fit-resolution")
            .add_event_listener_with_callback("change", cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    for (id, processing) in [("show-analysis", false), ("show-processing", true)] {
        let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            if processing {
                let _ = web_sys::window()
                    .unwrap()
                    .location()
                    .set_href("catalog.html#queue");
                return;
            }
            let (show, hide) = if processing {
                ("processing-view", "analysis-view")
            } else {
                ("analysis-view", "processing-view")
            };
            el(show).remove_attribute("hidden").unwrap();
            el(hide).set_attribute("hidden", "").unwrap();
        });
        el(id).add_event_listener_with_callback("click", cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    for param in 0..5 {
        let state = state.clone();
        let reset = Closure::<dyn FnMut(web_sys::Event)>::new(move |event: web_sys::Event| {
            event.prevent_default();
            let mut v = state.borrow_mut();
            let start = v.data.profiles[0].start_unix;
            let end = v.data.profiles.last().unwrap().end_unix;
            v.lower = 0;
            v.upper = v.data.profiles[v.time].gates.len() - 1;
            v.drag = None;
            selection(None);
            set_window(&mut v, start, end);
            render(&v);
        });
        el(&format!("fit-raster-{param}"))
            .add_event_listener_with_callback("contextmenu", reset.as_ref().unchecked_ref())?;
        reset.forget();
    }
    for param in 0..5 {
        for event in [
            "pointerdown",
            "pointerup",
            "pointermove",
            "pointercancel",
            "lostpointercapture",
        ] {
            let state = state.clone();
            let id = format!("fit-raster-{param}");
            let canvas_id = id.clone();
            let cb = Closure::<dyn FnMut(web_sys::PointerEvent)>::new(
                move |event_data: web_sys::PointerEvent| {
                    let rect = el(&canvas_id).get_bounding_client_rect();
                    let x =
                        ((event_data.client_x() as f64 - rect.x()) / rect.width()).clamp(0., 1.);
                    let y =
                        ((event_data.client_y() as f64 - rect.y()) / rect.height()).clamp(0., 1.);
                    let mut v = state.borrow_mut();
                    let t = v.window_start + x * (v.window_end - v.window_start);
                    if event == "pointercancel" || event == "lostpointercapture" {
                        v.drag = None;
                        selection(None);
                        return;
                    }
                    if event == "pointerdown" {
                        if event_data.button() != 0 {
                            return;
                        }
                        event_data.prevent_default();
                        let _ = el(&canvas_id).set_pointer_capture(event_data.pointer_id());
                        v.drag = Some((x, y, t));
                        return;
                    }
                    let ti = fit_index(&v.data, t);
                    let bottom = v.data.profiles[ti].gates[v.lower].lower_km;
                    let top = v.data.profiles[ti].gates[v.upper].upper_km;
                    let h = top - y * (top - bottom);
                    let gi = v.data.profiles[ti]
                        .gates
                        .partition_point(|g| g.upper_km <= h)
                        .min(v.upper);
                    if event == "pointermove" {
                        if let Some((x0, y0, t0)) = v.drag {
                            if event_data.buttons() == 0 {
                                v.drag = None;
                                selection(None);
                            } else {
                                selection(Some((x0, x, y0, y, t0, t)));
                            }
                        }
                        let p = &v.data.profiles[ti];
                        let g = &p.gates[gi];
                        let val = if param == 4 {
                            p.raw_ne.get(gi).and_then(|a| *a).map(|a| a[0])
                        } else {
                            if v.mask_bad && rejected(g, &v) {
                                None
                            } else {
                                g.result.as_ref().map(|f| value(&f.plasma, param))
                            }
                        };
                        text(
                            "fit-hover",
                            &format!(
                                "{}–{} UTC · {:.0}–{:.0} km · {} · {}",
                                clock(p.start_unix),
                                clock(p.end_unix),
                                g.lower_km,
                                g.upper_km,
                                val.map_or("NaN / unavailable".into(), |n| format!("{n:.4e}")),
                                if param == 4 {
                                    "raw 10 µs lag".to_string()
                                } else {
                                    format!(
                                        "{} · score {} · {}",
                                        g.status,
                                        g.result.as_ref().map_or(
                                            "unavailable".into(),
                                            |f| format!("{:.2}", f.reduced_chi_square)
                                        ),
                                        if rejected(g, &v) {
                                            format!(
                                                "rejected: {}",
                                                rejection_reasons(g, &v).join(", ")
                                            )
                                        } else {
                                            "passes display screening".to_string()
                                        }
                                    )
                                }
                            ),
                        );
                        return;
                    }
                    selection(None);
                    if let Some((x0, y0, t0)) = v.drag.take() {
                        if (x - x0).abs() > 0.005 || (y - y0).abs() > 0.02 {
                            if (y - y0).abs() > 0.02 {
                                let h0 = top - y0 * (top - bottom);
                                let gates = &v.data.profiles[ti].gates;
                                let lower = gates
                                    .partition_point(|g| g.upper_km <= h.min(h0))
                                    .min(gates.len() - 1);
                                let upper = gates
                                    .partition_point(|g| g.lower_km < h.max(h0))
                                    .saturating_sub(1)
                                    .max(lower)
                                    .min(gates.len() - 1);
                                v.lower = lower;
                                v.upper = upper;
                            }
                            if (x - x0).abs() > 0.005 {
                                set_window(&mut v, t.min(t0), t.max(t0));
                            }
                        } else {
                            v.time = ti;
                            v.gate = gi;
                        }
                        render(&v);
                        drop(v);
                        load_detail(state.clone());
                    }
                },
            );
            el(&id).add_event_listener_with_callback(event, cb.as_ref().unchecked_ref())?;
            cb.forget();
        }
    }
    let poll_state = state.clone();
    let study = location.study_status_url.clone();
    let poll = Closure::<dyn FnMut()>::new(move || {
        let state = poll_state.clone();
        let study = study.clone();
        spawn_local(refresh(state, study));
    });
    web_sys::window()
        .unwrap()
        .set_interval_with_callback_and_timeout_and_arguments_0(
            poll.as_ref().unchecked_ref(),
            15000,
        )?;
    poll.forget();
    spawn_local(refresh(state, location.study_status_url));
    Ok(())
}

fn selection(region: Option<(f64, f64, f64, f64, f64, f64)>) {
    let html=region.map_or(String::new(),|(x0,x1,y0,y1,t0,t1)|format!("<div style='position:absolute;left:calc(76px + (100% - 208px)*{});top:calc(22px + (100% - 80px)*{});width:calc((100% - 208px)*{});height:calc((100% - 80px)*{});border:1px solid #236b8e;background:#236b8e22'><span style='background:#fff;color:#236b8e'>{}–{} UTC</span></div>",x0.min(x1),y0.min(y1),(x1-x0).abs(),(y1-y0).abs(),clock(t0.min(t1)),clock(t0.max(t1))));
    for i in 0..5 {
        el(&format!("fit-selection-{i}")).set_inner_html(&html);
    }
}

fn move_settings(x: f64, y: f64) {
    let panel = el("fit-settings");
    let rect = panel.get_bounding_client_rect();
    let window = web_sys::window().unwrap();
    let width = window
        .inner_width()
        .ok()
        .and_then(|v| v.as_f64())
        .unwrap_or(1024.);
    let height = window
        .inner_height()
        .ok()
        .and_then(|v| v.as_f64())
        .unwrap_or(768.);
    let x = x.clamp(0., (width - rect.width()).max(0.));
    let y = y.clamp(0., (height - rect.height()).max(0.));
    let _ = panel.set_attribute("style", &format!("left:{x}px;top:{y}px;right:auto"));
}
fn settings_drag() -> Result<(), JsValue> {
    let drag = Rc::new(RefCell::new(None::<(i32, f64, f64)>));
    for kind in ["pointerdown", "pointermove", "pointerup", "pointercancel"] {
        let drag = drag.clone();
        let cb = Closure::<dyn FnMut(web_sys::PointerEvent)>::new(
            move |event: web_sys::PointerEvent| {
                let handle = el("fit-settings-drag");
                if kind == "pointerdown" {
                    if event.button() != 0 {
                        return;
                    }
                    let r = el("fit-settings").get_bounding_client_rect();
                    *drag.borrow_mut() = Some((
                        event.pointer_id(),
                        event.client_x() as f64 - r.x(),
                        event.client_y() as f64 - r.y(),
                    ));
                    let _ = handle.set_pointer_capture(event.pointer_id());
                    event.prevent_default();
                } else if let Some((pointer, dx, dy)) = *drag.borrow() {
                    if pointer != event.pointer_id() {
                        return;
                    }
                    if kind == "pointermove" {
                        move_settings(event.client_x() as f64 - dx, event.client_y() as f64 - dy);
                    } else {
                        let _ = handle.release_pointer_capture(pointer);
                    }
                }
                if kind == "pointerup" || kind == "pointercancel" {
                    *drag.borrow_mut() = None;
                }
            },
        );
        el("fit-settings-drag")
            .add_event_listener_with_callback(kind, cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    let keys =
        Closure::<dyn FnMut(web_sys::KeyboardEvent)>::new(move |event: web_sys::KeyboardEvent| {
            let step = if event.shift_key() { 40. } else { 10. };
            let (dx, dy) = match event.key().as_str() {
                "ArrowLeft" => (-step, 0.),
                "ArrowRight" => (step, 0.),
                "ArrowUp" => (0., -step),
                "ArrowDown" => (0., step),
                _ => return,
            };
            event.prevent_default();
            let r = el("fit-settings").get_bounding_client_rect();
            move_settings(r.x() + dx, r.y() + dy);
        });
    el("fit-settings-drag")
        .add_event_listener_with_callback("keydown", keys.as_ref().unchecked_ref())?;
    keys.forget();
    let resize = Closure::<dyn FnMut()>::new(move || {
        if !el("fit-settings").has_attribute("hidden") {
            let r = el("fit-settings").get_bounding_client_rect();
            move_settings(r.x(), r.y());
        }
    });
    web_sys::window()
        .unwrap()
        .add_event_listener_with_callback("resize", resize.as_ref().unchecked_ref())?;
    resize.forget();
    Ok(())
}

/// Intrinsic spectrum from the same stateless kinetic model used by the server.
fn analytical_spectrum(v: &Viewer) -> String {
    use guisdap_ion::spectrum::{IonSpecies, PlasmaState, SpectrumModel};
    let g = &selected(v).gates[v.gate];

    let Some(f) = &g.result else {
        return empty_curve(&g.status);
    };
    let model = SpectrumModel {
        radar_frequency_hz: 930e6,
        electron_collision_s: 0.,
        ions: vec![
            IonSpecies {
                mass_amu: 16.,
                fraction: g.oxygen_fraction,
            },
            IonSpecies {
                mass_amu: 30.5,
                fraction: 1. - g.oxygen_fraction,
            },
        ],
    };
    let frequency: Vec<_> = (0..=800).map(|i| -40000. + i as f64 * 100.).collect();
    let p = PlasmaState {
        electron_density_m3: f.plasma.electron_density_m3,
        ion_temperature_k: f.plasma.ion_temperature_k,
        electron_temperature_k: f.plasma.electron_temperature_k,
        velocity_m_s: f.plasma.velocity_m_s,
        ion_collision_s: f.plasma.ion_collision_s,
    };
    let evaluate = |p: PlasmaState| {
        model.evaluate(p, &frequency).map(|values| {
            values
                .into_iter()
                .map(|y| y * p.electron_density_m3)
                .collect::<Vec<_>>()
        })
    };
    let Ok(values) = evaluate(p) else {
        return empty_curve("Spectrum unavailable for this fit");
    };
    let n = f.free_parameter_indices.len();
    let mut gradients = Vec::new();
    if n > 0 && f.covariance.len() == n * n {
        let base = [
            p.electron_density_m3,
            p.ion_temperature_k,
            p.electron_temperature_k,
            p.velocity_m_s,
            p.ion_collision_s,
        ];
        for &index in &f.free_parameter_indices {
            if index >= 5 {
                break;
            }
            let delta = (base[index].abs() * 1e-3).max(if index == 0 { 1e5 } else { 0.1 });
            let change = |amount: f64| {
                let mut a = base;
                a[index] += amount;
                if index == 1 && g.settings.tie_electron_ion_temperature {
                    a[2] = a[1];
                }
                PlasmaState {
                    electron_density_m3: a[0],
                    ion_temperature_k: a[1],
                    electron_temperature_k: a[2],
                    velocity_m_s: a[3],
                    ion_collision_s: a[4],
                }
            };
            if let (Ok(plus), Ok(minus)) = (evaluate(change(delta)), evaluate(change(-delta))) {
                gradients.push(
                    plus.iter()
                        .zip(minus)
                        .map(|(a, b)| (a - b) / (2. * delta))
                        .collect::<Vec<_>>(),
                );
            } else {
                break;
            }
        }
    }
    let errors: Option<Vec<f64>> = if gradients.len() == n && n > 0 {
        Some(
            (0..values.len())
                .map(|k| {
                    let mut variance = 0.;
                    for i in 0..n {
                        for j in 0..n {
                            variance += gradients[i][k] * f.covariance[i * n + j] * gradients[j][k];
                        }
                    }
                    variance.max(0.).sqrt()
                })
                .collect(),
        )
    } else {
        None
    };
    let peak = values
        .iter()
        .enumerate()
        .map(|(i, y)| y + errors.as_ref().map_or(0., |e| e[i]))
        .fold(0., f64::max)
        .max(1.);
    let scale = 10f64.powf(peak.log10().floor());
    let top = peak * 1.08;
    let points: Vec<_> = frequency
        .iter()
        .zip(&values)
        .flat_map(|(f, v)| [*f / 1000., *v / scale])
        .collect();
    let band: Vec<_> = errors
        .iter()
        .flat_map(|e| {
            (0..values.len()).flat_map(|i| {
                [
                    frequency[i] / 1000.,
                    (values[i] - e[i]).max(0.) / scale,
                    (values[i] + e[i]) / scale,
                ]
            })
        })
        .collect();
    serde_json::json!({"title":"Analytical ISR spectrum","bounds":[-40.,40.,0.,top/scale],"xlabel":"Frequency offset / kHz","ylabel":format!("Ne S(f) / {scale:.0e} m⁻³ Hz⁻¹"),"series":[{"xy":points}],"bands":[{"values":band,"color":"#236b8e4d"}],"caption":if errors.is_some(){"Intrinsic two-ion model · ±1σ parameter uncertainty"}else{"Intrinsic two-ion model · uncertainty unavailable"}}).to_string()
}
fn empty_curve(message: &str) -> String {
    serde_json::json!({"bounds":[0,1,0,1],"caption":message}).to_string()
}
