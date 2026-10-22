//! Explicit user-selected plasma-line calibration; no automatic rescaling.
use super::*;
use guisdap_core::calibration::{estimate, frequency_hz, Estimate};
use serde::Serialize;
#[derive(Deserialize, Clone)]
struct Settings {
    ion_index_url: String,
    #[serde(default)]
    experiment_id: Option<String>,
    radar_frequency_hz: f64,
    calibration_multiplier: f64,
    #[serde(default)]
    analysis_config_url: Option<String>,
}
#[derive(Deserialize)]
struct Response {
    frequency_hz: Vec<f64>,
    kernel_real: Vec<Vec<f64>>,
}
#[derive(Deserialize)]
struct Fits {
    profiles: Vec<Profile>,
}
#[derive(Deserialize, Clone)]
struct Profile {
    start_unix: f64,
    end_unix: f64,
    gates: Vec<Gate>,
    detail_url: Option<String>,
}
#[derive(Deserialize, Clone)]
struct Gate {
    oxygen_fraction: f64,
    altitude_km: f64,
    lower_km: f64,
    upper_km: f64,
    status: String,
    result: Option<Fit>,
}
#[derive(Deserialize, Clone)]
struct Fit {
    plasma: Plasma,
    standard_error: [Option<f64>; 5],
    reduced_chi_square: f64,
}
#[derive(Deserialize, Clone)]
struct Plasma {
    velocity_m_s: f64,
    ion_collision_s: f64,
    electron_density_m3: f64,
    electron_temperature_k: f64,
    ion_temperature_k: f64,
}
#[derive(Serialize, Deserialize, Clone)]
struct Point {
    #[serde(default)]
    raw_product_url: String,
    raw_ne_proxy_m3: f64,
    equivalent_zero_lag_raw_ne_m3: f64,
    short_lag_response: f64,
    power_time_window_s: f64,
    power_altitude_width_km: f64,
    raw_sample_count: usize,
    time_unix: f64,
    plasma_start_unix: f64,
    plasma_end_unix: f64,
    broad_gate: usize,
    frequency_hz: f64,
    frequency_bin_width_hz: f64,
    altitude_km: f64,
    altitude_bounds_km: [f64; 2],
    ion_start_unix: f64,
    ion_end_unix: f64,
    ion_product_url: String,
    radar_frequency_hz: f64,
    original_calibration_multiplier: f64,
    #[serde(default)]
    magic_constant: f64,
    te_k: f64,
    ti_k: f64,
    fitted_ne_m3: f64,
    fit_standard_error: [Option<f64>; 5],
    fit_status: String,
    residual_score: f64,
    estimate: Estimate,
    included: bool,
}
#[derive(Clone)]
struct Pending {
    time: f64,
    start: f64,
    end: f64,
    broad: usize,
    frequency: f64,
    profile: Profile,
}
pub(super) struct Calibration {
    sync: Rc<RefCell<Sync>>,
    suppress_click: bool,
    app: State,
    response: Option<Response>,
    editing: Option<usize>,
    settings: Option<Settings>,
    fits: Option<Fits>,
    pending: Option<Pending>,
    points: Vec<Point>,
    request: u64,
}
pub(super) type Cal = Rc<RefCell<Calibration>>;
const KEY: &str = "guisdap26-beata-calibration-v1";
fn input(id: &str) -> HtmlInputElement {
    el(id).dyn_into().unwrap()
}
fn save(c: &Calibration) {
    {
        let mut sync = c.sync.borrow_mut();
        sync.pending = Some(serde_json::to_value(&c.points).unwrap());
    }
    flush(c.sync.clone());
    if let Ok(Some(storage)) = web_sys::window().unwrap().local_storage() {
        if let Ok(s) = serde_json::to_string(&c.points) {
            let _ = storage.set_item(KEY, &s);
        }
    }
}
fn safe_fit(g: &Gate) -> bool {
    let Some(f) = &g.result else {
        return false;
    };
    g.status == "converged"
        && [0, 1, 2]
            .iter()
            .all(|&i| f.standard_error[i].is_some_and(|x| x.is_finite() && x >= 0.))
        && f.standard_error[0].unwrap() <= 0.5 * f.plasma.electron_density_m3
        && f.standard_error[1].unwrap() <= f.plasma.ion_temperature_k
        && f.standard_error[2].unwrap() <= f.plasma.electron_temperature_k
}
fn point(c: &Calibration) -> Result<Point, String> {
    let s = c
        .settings
        .as_ref()
        .ok_or("Calibration settings unavailable")?;
    let p = c
        .pending
        .as_ref()
        .ok_or("Click a broad-gate spectrum first")?;
    let i = el("cal-altitude")
        .dyn_into::<HtmlSelectElement>()
        .unwrap()
        .value()
        .parse::<usize>()
        .map_err(|_| "Select the source altitude")?;
    let g = p.profile.gates.get(i).ok_or("Altitude unavailable")?;
    let f = g.result.as_ref().ok_or("No fit at that altitude")?;
    let height = input("cal-height")
        .value()
        .parse::<f64>()
        .map_err(|_| "Invalid altitude")?;
    let freq = input("cal-frequency")
        .value()
        .parse::<f64>()
        .map_err(|_| "Invalid frequency")?
        * 1e6;
    // Frequency sign is retained, but the unmagnetized dispersion uses its square.
    if !(-6.15e6..=-1.275e6).contains(&freq) {
        return Err("Frequency is outside this event's recorded plasma bands".into());
    }
    let mut result = estimate(
        freq,
        s.radar_frequency_hz,
        f.plasma.electron_temperature_k,
        f.plasma.ion_temperature_k,
        f.plasma.electron_density_m3,
    )
    .map_err(|e| e.to_string())?;
    let response = c
        .response
        .as_ref()
        .ok_or("Short-lag response unavailable")?;
    let time_width = input("cal-time-width")
        .value()
        .parse::<f64>()
        .map_err(|_| "Invalid time width")?;
    let range_width = input("cal-range-width")
        .value()
        .parse::<f64>()
        .map_err(|_| "Invalid altitude width")?;
    if !time_width.is_finite() || time_width <= 0. || !range_width.is_finite() || range_width <= 0.
    {
        return Err("Positive averaging widths required".into());
    }
    let model = guisdap_ion::spectrum::SpectrumModel {
        radar_frequency_hz: s.radar_frequency_hz,
        electron_collision_s: 0.,
        ions: vec![
            guisdap_ion::spectrum::IonSpecies {
                mass_amu: 16.,
                fraction: g.oxygen_fraction,
            },
            guisdap_ion::spectrum::IonSpecies {
                mass_amu: 30.5,
                fraction: 1. - g.oxygen_fraction,
            },
        ],
    };
    let spectrum = model
        .evaluate(
            guisdap_ion::spectrum::PlasmaState {
                electron_density_m3: result.plasma_density_m3,
                electron_temperature_k: p.profile.gates[i]
                    .result
                    .as_ref()
                    .unwrap()
                    .plasma
                    .electron_temperature_k,
                ion_temperature_k: f.plasma.ion_temperature_k,
                velocity_m_s: f.plasma.velocity_m_s,
                ion_collision_s: f.plasma.ion_collision_s,
            },
            &response.frequency_hz,
        )
        .map_err(|e| e.to_string())?;
    let app = c.app.borrow();
    let n = app.m.frames.len();
    let cadence = app.m.frames[0].end_unix - app.m.frames[0].start_unix;
    if time_width + 0.01 < cadence {
        return Err(format!(
            "Choose at least {cadence:.0} s or switch the plasma viewer to a finer resolution"
        ));
    }
    let mut sum = 0.;
    let mut weight = 0.;
    let mut expected = 0.;
    let mut count = 0;
    for (t, frame) in app.m.frames.iter().enumerate() {
        let midpoint = (frame.start_unix + frame.end_unix) / 2.;
        if midpoint < p.time - time_width / 2. || midpoint >= p.time + time_width / 2. {
            continue;
        }
        for (j, r) in app.m.ion_range_km.iter().enumerate() {
            let target_height = height;
            let height = altitude(*r, frame.elevation_deg);
            if (height - target_height).abs() > range_width / 2. {
                continue;
            }
            let raw = app.power[j * n + t] as f64;
            if !raw.is_finite() {
                continue;
            }
            let w = frame.end_unix - frame.start_unix;
            let kernel = response
                .kernel_real
                .get(j)
                .ok_or("Missing native response")?;
            let response_value = kernel
                .iter()
                .zip(&spectrum)
                .map(|(k, v)| k * v)
                .sum::<f64>();
            sum += w * raw;
            expected += w * response_value;
            weight += w;
            count += 1;
        }
    }
    if count == 0 || sum <= 0. {
        return Err("No positive measured short-lag power in the selected window".into());
    }
    let raw = sum / weight;
    let response_value = expected / weight;
    if response_value <= 0. || !response_value.is_finite() {
        return Err(
            "Short-lag model response is nonpositive; conversion is ill-conditioned".into(),
        );
    }
    let zero = guisdap_core::calibration::short_lag_conversion(&mut result, raw, response_value)
        .map_err(|e| e.to_string())?;
    Ok(Point {
        raw_product_url: format!("{}{}ion-short-raw-ne.f32", app.root, app.m.overview_root),
        raw_ne_proxy_m3: raw,
        equivalent_zero_lag_raw_ne_m3: zero,
        short_lag_response: response_value,
        power_time_window_s: time_width,
        power_altitude_width_km: range_width,
        raw_sample_count: count,
        time_unix: p.time,
        plasma_start_unix: p.start,
        plasma_end_unix: p.end,
        broad_gate: p.broad + 1,
        frequency_hz: freq,
        frequency_bin_width_hz: 25000.,
        altitude_km: height,
        altitude_bounds_km: [g.lower_km, g.upper_km],
        ion_start_unix: p.profile.start_unix,
        ion_end_unix: p.profile.end_unix,
        ion_product_url: p
            .profile
            .detail_url
            .clone()
            .map(|d| format!("{}{d}", s.ion_index_url.trim_end_matches("index.json")))
            .unwrap_or_else(|| s.ion_index_url.clone()),
        radar_frequency_hz: s.radar_frequency_hz,
        original_calibration_multiplier: s.calibration_multiplier,
        magic_constant: s.calibration_multiplier * result.calibration_multiplier_ratio,
        te_k: f.plasma.electron_temperature_k,
        ti_k: f.plasma.ion_temperature_k,
        fitted_ne_m3: f.plasma.electron_density_m3,
        fit_standard_error: f.standard_error,
        fit_status: g.status.clone(),
        residual_score: f.reduced_chi_square,
        included: safe_fit(g) && result.k2_debye2 < 0.1,
        estimate: result,
    })
}
fn preview(c: &Calibration) {
    match point(c) {
        Err(e) => {
            text("cal-result", &e);
            text("cal-result-details", "");
        }
        Ok(p) => {
            let a = &p.estimate;
            if let Ok(Some(svg)) = el("cal-candidates").query_selector("svg") {
                if let Ok(Some(old)) = svg.query_selector("#cal-height-marker") {
                    old.remove();
                }
                let low = input("cal-height").min().parse::<f64>().unwrap_or(0.);
                let high = input("cal-height").max().parse::<f64>().unwrap_or(500.);
                let y = 255. - 210. * (p.altitude_km - low) / (high - low);
                let marker = doc()
                    .create_element_ns(Some("http://www.w3.org/2000/svg"), "path")
                    .unwrap();
                let _ = marker.set_attribute("id", "cal-height-marker");
                let _ = marker.set_attribute("d", &format!("M60,{y} H610"));
                let _ = marker.set_attribute("stroke", "#a9fff3");
                let _ = svg.append_child(&marker);
            }
            let ratio = p.magic_constant / p.original_calibration_multiplier;
            text("cal-result", &format!("Input-product density is {} (C = {:.4}).\nThis point estimates C = {:.4}: {} electron density by {:.1}% (×{:.3}).\n{:.1} km · {}–{} UTC · {}", if ratio >= 1. { "underestimated" } else { "overestimated" }, p.original_calibration_multiplier, p.magic_constant, if ratio >= 1. { "increase" } else { "decrease" }, (ratio - 1.).abs()*100., ratio, p.altitude_km, utc(p.ion_start_unix), utc(p.ion_end_unix), if p.included { "eligible for inclusion" } else { "excluded by default: review fit quality" }));
            text("cal-result-details",&format!("Altitude {:.1} km · one-minute fit {}–{} UTC · status {} · residual {:.2}\nTe {:.0} ± {} K; Ti {:.0} ± {} K; Te/Ti {:.3}\nCold Ne {:.4e}; thermal-corrected Ne {:.4e} m⁻³\nThermal frequency term {:.3} MHz²; Debye length {:.4} m; k²λD² {:.4}; power denominator {:.4}\nFitted Ne {:.4e} m⁻³; fitted-power denominator {:.4}\nEstimated magic constant {:.4} (dimensionless)\nInput product magic constant {:.4}; relative density correction ×{:.4}\n{}",
                p.altitude_km,utc(p.ion_start_unix),utc(p.ion_end_unix),p.fit_status,p.residual_score,p.te_k,p.fit_standard_error[2].map(|x|format!("{x:.0}")).unwrap_or("unknown".into()),p.ti_k,p.fit_standard_error[1].map(|x|format!("{x:.0}")).unwrap_or("unknown".into()),a.temperature_ratio,a.cold_density_m3,a.plasma_density_m3,a.thermal_frequency_squared_hz2/1e12,a.debye_length_m,a.k2_debye2,a.ion_power_denominator,p.fitted_ne_m3,a.fitted_ion_power_denominator,p.magic_constant,p.original_calibration_multiplier,a.calibration_multiplier_ratio,
                if p.included {"Eligible for inclusion. Confirm the altitude association before saving."} else {"Flagged fit or k²λD² ≥ 0.1: saved excluded by default. Thermal approximation / fit quality requires review."}));
            let current = el("cal-result-details").text_content().unwrap_or_default();
            text("cal-result-details",&format!("{current}\nMeasured raw 10 µs Ne proxy {:.4e}; equivalent zero-lag raw Ne {:.4e} m⁻³\nModel short-lag response {:.5}; {} samples · {:.0} s × {:.1} km requested window",p.raw_ne_proxy_m3,p.equivalent_zero_lag_raw_ne_m3,p.short_lag_response,p.raw_sample_count,p.power_time_window_s,p.power_altitude_width_km));
        }
    }
    markers(&c.app.borrow(), c);
}
fn table(state: &Cal) {
    let c = state.borrow();
    let mut html=String::from("<table><thead><tr><th>Use</th><th>UTC</th><th>Gate</th><th>MHz</th><th>km</th><th>Te/Ti</th><th>Ne / m⁻³</th><th>Magic constant</th><th>Fit</th><th></th></tr></thead><tbody>");
    for (i, p) in c.points.iter().enumerate() {
        html.push_str(&format!("<tr><td><input id='cal-use-{i}' aria-label='Include calibration point {}' type='checkbox' {}></td><td>{}</td><td>{}</td><td>{:.3}</td><td>{:.1}</td><td>{:.2}</td><td>{:.3e}</td><td>{:.4}</td><td>{}</td><td><button id='cal-edit-{i}'>Edit altitude</button><button id='cal-remove-{i}'>Remove {}</button></td></tr>",i+1,if p.included {"checked"}else{""},utc(p.time_unix),p.broad_gate,p.frequency_hz/1e6,p.altitude_km,p.estimate.temperature_ratio,p.estimate.plasma_density_m3,p.magic_constant,p.fit_status,i+1));
    }
    html.push_str("</tbody></table>");
    el("cal-points").set_inner_html(&html);
    let mut values: Vec<_> = c
        .points
        .iter()
        .filter(|p| p.included)
        .map(|p| p.magic_constant)
        .collect();
    values.sort_by(f64::total_cmp);
    text(
        "cal-summary",
        &if values.is_empty() {
            "No included calibration points.".into()
        } else {
            let n = values.len();
            let median = (values[(n - 1) / 2] + values[n / 2]) / 2.;
            format!("{} included points · median magic constant {:.4} · range {:.4}–{:.4}. Repeat clicks are not independent measurements. No processing settings have been changed.",n,median,values[0],values[n-1])
        },
    );
    let export = serde_json::json!({"schema":"guisdap26-plasma-calibration-v1","method":"Unmagnetized Bohm-Gross; finite-Debye ion power Ne/[(1+q)(1+q+Te/Ti)]. Measured 10-us raw Ne converted with the experiment spectral response and fitted shape at plasma-line density; finite-Debye equivalent power. No raw short-lag zero-lag substitution; no drift correction. Approximation, not a kinetic dispersion inversion.","points":c.points,"magic_constant_fit":guisdap_core::calibration::fit_gain(&c.points.iter().filter(|p|p.included).map(|p|(p.raw_ne_proxy_m3/p.original_calibration_multiplier,p.raw_ne_proxy_m3*p.estimate.calibration_multiplier_ratio,p.ion_start_unix.to_bits())).collect::<Vec<_>>())});
    let uri = js_sys::encode_uri_component(&serde_json::to_string_pretty(&export).unwrap());
    let _ = el("cal-export").set_attribute(
        "href",
        &format!("data:application/json;charset=utf-8,{uri}"),
    );
    let lower = c
        .points
        .iter()
        .filter(|p| p.magic_constant > 0.)
        .map(|p| p.magic_constant.log10())
        .fold(0., f64::min)
        .floor()
        - 0.2;
    let upper = c
        .points
        .iter()
        .filter(|p| p.magic_constant > 0.)
        .map(|p| p.magic_constant.log10())
        .fold(0., f64::max)
        .ceil()
        + 0.2;
    let dots = |included: bool| {
        c.points
            .iter()
            .filter(|p| p.included == included && p.magic_constant > 0.)
            .flat_map(|p| [p.magic_constant.log10(), p.altitude_km])
            .collect::<Vec<_>>()
    };
    crate::plotting::shared_curve("cal-accumulated",&serde_json::json!({"title":"Magic constant versus altitude","bounds":[lower,upper,40.,480.],"xlabel":"Magic constant","xformat":"log10","ylabel":"Altitude / km","series":[{"xy":dots(true),"style":0,"radius":4},{"xy":dots(false),"style":1,"radius":4,"color":"#986327"}]}).to_string());
    let lo = c
        .points
        .iter()
        .flat_map(|p| {
            [
                p.raw_ne_proxy_m3 / p.original_calibration_multiplier,
                p.raw_ne_proxy_m3 * p.estimate.calibration_multiplier_ratio,
            ]
        })
        .filter(|v| *v > 0.)
        .map(f64::log10)
        .fold(11., f64::min)
        .floor()
        - 0.1;
    let hi = c
        .points
        .iter()
        .flat_map(|p| {
            [
                p.raw_ne_proxy_m3 / p.original_calibration_multiplier,
                p.raw_ne_proxy_m3 * p.estimate.calibration_multiplier_ratio,
            ]
        })
        .filter(|v| *v > 0.)
        .map(f64::log10)
        .fold(11., f64::max)
        .ceil()
        + 0.1;
    let samples: Vec<_> = c
        .points
        .iter()
        .filter(|p| p.included)
        .map(|p| {
            (
                p.raw_ne_proxy_m3 / p.original_calibration_multiplier,
                p.raw_ne_proxy_m3 * p.estimate.calibration_multiplier_ratio,
                p.ion_start_unix.to_bits(),
            )
        })
        .collect();
    let gain_fit = guisdap_core::calibration::fit_gain(&samples);
    let scatter = |included: bool| {
        c.points
            .iter()
            .filter(|p| p.included == included)
            .flat_map(|p| {
                [
                    (p.raw_ne_proxy_m3 / p.original_calibration_multiplier).log10(),
                    (p.raw_ne_proxy_m3 * p.estimate.calibration_multiplier_ratio).log10(),
                ]
            })
            .collect::<Vec<_>>()
    };
    let line: Vec<_> = gain_fit
        .iter()
        .flat_map(|fit| [lo, lo + fit.gain.log10(), hi, hi + fit.gain.log10()])
        .collect();
    let band: Vec<_> = gain_fit
        .iter()
        .filter_map(|f| f.standard_error.map(|e| (f.gain, e)))
        .flat_map(|(g, e)| {
            [
                lo,
                lo + (g - e).max(1e-30).log10(),
                lo + (g + e).log10(),
                hi,
                hi + (g - e).max(1e-30).log10(),
                hi + (g + e).log10(),
            ]
        })
        .collect();
    crate::plotting::shared_curve("cal-power-scatter",&serde_json::json!({"title":"Predicted versus measured short-lag power","bounds":[lo,hi,lo,hi],"xlabel":"Measured raw Ne at C=1 / m⁻³","ylabel":"Predicted raw Ne / m⁻³","xformat":"log10","yformat":"log10","series":[{"xy":[lo,lo,hi,hi],"color":"#7b929e"},{"xy":line},{"xy":scatter(true),"style":0,"radius":4},{"xy":scatter(false),"style":1,"radius":4,"color":"#986327"}],"bands":[{"values":band,"color":"#236b8e4d"}]}).to_string());
    text(
        "cal-fit-summary",
        &match gain_fit {
            Some(f) => {
                let correction = c.settings.as_ref().map(|s| {
                let ratio=f.gain/s.calibration_multiplier;
                format!(" Compared with the input constant {:.4}, {} electron density by {:.1}% (×{:.3}); the old density is {}.", s.calibration_multiplier, if ratio>=1. {"increase"} else {"decrease"}, (ratio-1.).abs()*100., ratio, if ratio>=1. {"underestimated"} else {"overestimated"})
            }).unwrap_or_default();
                format!("Magic constant C = {:.4}{}.{correction} {} included points from {} one-minute intervals. Uncertainty is ±1σ and excludes systematic errors. Processing uses a saved snapshot; editing points does not change existing results.",f.gain,f.standard_error.map(|se|format!(" ± {se:.4}")).unwrap_or_default(),f.points,f.time_clusters)
            }
            None => "Include calibration points to estimate the absolute magic constant.".into(),
        },
    );
    markers(&c.app.borrow(), &c);
    let n = c.points.len();
    save(&c);
    drop(c);
    for i in 0..n {
        let st = state.clone();
        let _ = on(&format!("cal-edit-{i}"), "click", move |_| {
            let mut c = st.borrow_mut();
            let saved = c.points[i].clone();
            let Some(profile) = c
                .fits
                .as_ref()
                .and_then(|f| {
                    f.profiles
                        .iter()
                        .find(|p| p.start_unix == saved.ion_start_unix)
                })
                .cloned()
            else {
                text("cal-result","Select a broad-gate point first to load the one-minute fits, then edit this saved point.");
                return;
            };
            c.pending = Some(Pending {
                time: saved.time_unix,
                start: saved.plasma_start_unix,
                end: saved.plasma_end_unix,
                broad: saved.broad_gate - 1,
                frequency: saved.frequency_hz,
                profile,
            });
            c.editing = Some(i);
            input("cal-frequency").set_value(&(saved.frequency_hz / 1e6).to_string());
            input("cal-time-width").set_value(&saved.power_time_window_s.to_string());
            input("cal-range-width").set_value(&saved.power_altitude_width_km.to_string());
            let app = c.app.borrow();
            let elevation = app.m.frames[0].elevation_deg;
            drop(app);
            let bounds = [[47., 221.], [128., 302.], [209., 383.], [290., 464.]]
                [saved.broad_gate - 1]
                .map(|r| altitude(r, elevation));
            candidates(&c, bounds);
            let gate = c
                .pending
                .as_ref()
                .unwrap()
                .profile
                .gates
                .iter()
                .position(|g| (g.altitude_km - saved.altitude_km).abs() < 0.01)
                .unwrap_or(0);
            el("cal-altitude")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .set_value(&gate.to_string());
            input("cal-height").set_value(&saved.altitude_km.to_string());
            text(
                "cal-selection",
                &format!("Editing point {} · {} UTC", i + 1, utc(saved.time_unix)),
            );
            preview(&c);
        });
        let st = state.clone();
        let _ = on(&format!("cal-use-{i}"), "change", move |_| {
            st.borrow_mut().points[i].included = input(&format!("cal-use-{i}")).checked();
            table(&st);
        });
        let st = state.clone();
        let _ = on(&format!("cal-remove-{i}"), "click", move |_| {
            st.borrow_mut().points.remove(i);
            st.borrow_mut().editing = None;
            table(&st);
        });
    }
}
fn candidates(c: &Calibration, bounds: [f64; 2]) {
    let Some(p) = &c.pending else {
        return;
    };
    input("cal-height").set_min(&bounds[0].to_string());
    input("cal-height").set_max(&bounds[1].to_string());
    input("cal-height").set_value(&((bounds[0] + bounds[1]) / 2.).to_string());
    let s = c.settings.as_ref().unwrap();
    let mut options = String::from("<option value=''>Select source altitude…</option>");
    let mut good = Vec::new();
    let mut bad = Vec::new();
    for (i, g) in p.profile.gates.iter().enumerate() {
        if g.altitude_km < bounds[0] || g.altitude_km > bounds[1] {
            continue;
        }
        let Some(f) = &g.result else {
            continue;
        };
        let predicted = frequency_hz(
            f.plasma.electron_density_m3,
            f.plasma.electron_temperature_k,
            s.radar_frequency_hz,
        ) / 1e6;
        options.push_str(&format!(
            "<option value='{i}'>{:.1} km · predicted {:.3} MHz · {}</option>",
            g.altitude_km, predicted, g.status
        ));
        if (1. ..=7.).contains(&predicted) {
            if safe_fit(g) {
                good.extend([predicted, g.altitude_km]);
            } else {
                bad.extend([predicted, g.altitude_km]);
            }
        }
    }
    crate::plotting::shared_curve("cal-candidates",&serde_json::json!({"title":"Candidate plasma frequency","bounds":[1.,7.,bounds[0],bounds[1]],"xlabel":"|Frequency offset| / MHz","ylabel":"Altitude / km","series":[{"xy":good,"style":0,"radius":4},{"xy":bad,"style":1,"radius":4,"color":"#986327"},{"xy":[p.frequency.abs()/1e6,bounds[0],p.frequency.abs()/1e6,bounds[1]],"color":"#986327"}]}).to_string());
    el("cal-altitude").set_inner_html(&options);
}
pub(super) fn start(app: State) -> Result<(), JsValue> {
    let mut points: Vec<Point> = web_sys::window()
        .unwrap()
        .local_storage()
        .ok()
        .flatten()
        .and_then(|s| s.get_item(KEY).ok().flatten())
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_default();
    for p in &mut points {
        p.magic_constant =
            p.original_calibration_multiplier * p.estimate.calibration_multiplier_ratio;
    }
    let storage = web_sys::window()
        .unwrap()
        .local_storage()?
        .ok_or_else(|| error("Browser storage unavailable"))?;
    let key = if let Some(k) = storage.get_item("guisdap26-calibration-notebook")? {
        k
    } else {
        let mut bytes = [0u8; 32];
        web_sys::window()
            .unwrap()
            .crypto()?
            .get_random_values_with_u8_array(&mut bytes)?;
        let k = bytes.iter().map(|b| format!("{b:02x}")).collect::<String>();
        storage.set_item("guisdap26-calibration-notebook", &k)?;
        k
    };
    let sync = Rc::new(RefCell::new(Sync {
        key,
        revision: 0,
        ready: false,
        busy: false,
        pending: None,
    }));
    let state = Rc::new(RefCell::new(Calibration {
        sync: sync.clone(),
        suppress_click: false,
        app: app.clone(),
        response: None,
        editing: None,
        settings: None,
        fits: None,
        pending: None,
        points,
        request: 0,
    }));
    app.borrow_mut().calibration = Some(state.clone());
    gestures(&state)?;
    recalculation(&state)?;
    table(&state);
    let loaded = state.clone();
    spawn_local(async move {
        let key = loaded.borrow().sync.borrow().key.clone();
        let result: Result<(), JsValue> = async {
            let n: Saved = json(&format!("api/calibrations/{key}")).await?;
            {
                let mut c = loaded.borrow_mut();
                if n.revision > 0 {
                    c.points = n.points;
                    for p in &mut c.points {
                        p.magic_constant = p.original_calibration_multiplier
                            * p.estimate.calibration_multiplier_ratio;
                    }
                }
                let mut sync = c.sync.borrow_mut();
                sync.revision = n.revision;
                sync.ready = true;
            }
            let settings: Settings = json("calibration-config.json").await?;
            if let Some(url)=&settings.analysis_config_url {
                let summary: Result<String,JsValue> = async {
                    let cfg: serde_json::Value=json(url).await?;
                    let dataset=cfg["datasets"].as_array().and_then(|d|d.iter().find(|d|d["integration_s"]==60)).ok_or_else(||error("No one-minute product"))?;
                    let index: serde_json::Value=json(dataset["data_url"].as_str().ok_or_else(||error("Missing product URL"))?).await?;
                    let constant=index["calibration"]["magic_constant"].as_f64().ok_or_else(||error("Missing calibration provenance"))?;
                    Ok(format!("Current one-minute ion analysis: C = {constant:.5}. Calibration reference measurements below retain C = {:.4}; these are the original data used to estimate the new constant.",settings.calibration_multiplier))
                }.await;
                text("cal-analysis-used",&summary.unwrap_or_else(|_|"Current ion-analysis calibration unavailable; open Ion-line to inspect product provenance.".into()));
            }
            let fits: Fits = json(&settings.ion_index_url).await?;
            let response: Response = json("data/beata-short-lag-response.json").await?;
            {
                let mut c = loaded.borrow_mut();
                c.settings = Some(settings);
                c.fits = Some(fits);
                c.response = Some(response);
            }
            table(&loaded);
            Ok(())
        }
        .await;
        if let Err(e) = result {
            text(
                "cal-save-status",
                &format!("Server unavailable; browser copy retained. {e:?}"),
            );
        }
    });
    for id in [
        "cal-altitude",
        "cal-frequency",
        "cal-time-width",
        "cal-range-width",
        "cal-height",
    ] {
        let s = state.clone();
        on(id, "input", move |_| {
            let c = s.borrow();
            if let Some(p) = &c.pending {
                if id == "cal-altitude" {
                    if let Ok(i) = el("cal-altitude")
                        .dyn_into::<HtmlSelectElement>()
                        .unwrap()
                        .value()
                        .parse::<usize>()
                    {
                        input("cal-height").set_value(&p.profile.gates[i].altitude_km.to_string());
                    }
                }
                if id == "cal-height" {
                    let h = input("cal-height").value().parse::<f64>().unwrap_or(0.);
                    let index = p
                        .profile
                        .gates
                        .iter()
                        .enumerate()
                        .filter(|(_, g)| g.result.is_some())
                        .min_by(|(_, a), (_, b)| {
                            (a.altitude_km - h)
                                .abs()
                                .total_cmp(&(b.altitude_km - h).abs())
                        })
                        .map(|(i, _)| i);
                    if let Some(i) = index {
                        el("cal-altitude")
                            .dyn_into::<HtmlSelectElement>()
                            .unwrap()
                            .set_value(&i.to_string());
                    }
                }
            }
            drop(c);
            update_current(&s);
        })?;
    }
    for event_name in ["pointerdown", "pointermove"] {
        let st = state.clone();
        on("cal-candidates", event_name, move |event| {
            let Ok(e) = event.dyn_into::<web_sys::MouseEvent>() else {
                return;
            };
            if event_name == "pointermove" && e.buttons() != 1 {
                return;
            }
            let c = st.borrow();
            let Some(p) = &c.pending else {
                return;
            };
            let rect = el("cal-candidates")
                .query_selector("svg")
                .unwrap()
                .unwrap()
                .get_bounding_client_rect();
            let fraction = ((e.client_y() as f64 - rect.top()) / rect.height() * 300. - 45.) / 210.;
            let app = c.app.borrow();
            let elevation = app.m.frames[0].elevation_deg;
            let bounds = [[47., 221.], [128., 302.], [209., 383.], [290., 464.]][p.broad]
                .map(|r| altitude(r, elevation));
            let h = bounds[1] - fraction.clamp(0., 1.) * (bounds[1] - bounds[0]);
            let index = p
                .profile
                .gates
                .iter()
                .enumerate()
                .filter(|(_, g)| {
                    g.result.is_some() && g.altitude_km >= bounds[0] && g.altitude_km <= bounds[1]
                })
                .min_by(|(_, a), (_, b)| {
                    (a.altitude_km - h)
                        .abs()
                        .total_cmp(&(b.altitude_km - h).abs())
                })
                .map(|(i, _)| i);
            if let Some(i) = index {
                input("cal-height").set_value(&h.to_string());
                el("cal-altitude")
                    .dyn_into::<HtmlSelectElement>()
                    .unwrap()
                    .set_value(&i.to_string());
                drop(app);
                drop(c);
                update_current(&st);
            }
        })?;
    }
    let s = state.clone();
    on("cal-add", "click", move |_| {
        let p = point(&s.borrow());
        match p {
            Ok(p) => {
                let mut c = s.borrow_mut();
                if let Some(i) = c.editing.take() {
                    c.points[i] = p;
                } else {
                    c.points.push(p);
                }
                drop(c);
                table(&s);
                text(
                    "cal-result",
                    "Point saved. Click another feature or select another altitude.",
                );
            }
            Err(e) => {
                text("cal-result", &e);
                text("cal-result-details", "");
            }
        }
    })?;
    for broad in 0..4 {
        let st = state.clone();
        let app = app.clone();
        on(&format!("broad{broad}"), "click", move |event| {
            if st.borrow().suppress_click {
                st.borrow_mut().suppress_click = false;
                return;
            }
            if !input("cal-enabled").checked() {
                return;
            }
            let Ok(event) = event.dyn_into::<web_sys::MouseEvent>() else {
                return;
            };
            let rect = el(&format!("broad{broad}")).get_bounding_client_rect();
            let x = ((event.client_x() as f64 - rect.left()) / rect.width())
                .clamp(0., 1. - f64::EPSILON);
            let y = (1. - (event.client_y() as f64 - rect.top()) / rect.height())
                .clamp(0., 1. - f64::EPSILON);
            let a = app.borrow();
            let dt = a.m.frames[0].end_unix - a.m.frames[0].start_unix;
            let columns = ((a.view_end - a.view_start) / dt).ceil().max(1.);
            let t = a.view_start
                + ((x * columns).floor() + 0.5) / columns * (a.view_end - a.view_start);
            let Some(frame) =
                a.m.frames
                    .iter()
                    .find(|f| f.start_unix <= t && f.end_unix > t)
            else {
                text("cal-result", "No plasma data at this time (gap).");
                return;
            };
            let start = frame.start_unix;
            let end = frame.end_unix;
            let elevation = frame.elevation_deg;
            let centers = frame.centers_hz;
            let low = centers[0].min(centers[1]) - 1.25e6;
            let high = centers[0].max(centers[1]) + 1.225e6;
            let first = ((a.broad_frequency[0] - low) / 25000.).ceil().max(0.);
            let last = ((a.broad_frequency[1].min(high) - low) / 25000.).floor();
            let freq = low + (first + (y * (last - first + 1.)).floor()) * 25000.;
            let selected =
                a.m.frames
                    .iter()
                    .position(|f| f.start_unix == start)
                    .unwrap();
            let hit = st.borrow().points.iter().position(|p| {
                p.broad_gate == broad + 1
                    && ((p.time_unix - t) / (a.view_end - a.view_start) * rect.width()).abs() < 10.
                    && ((p.frequency_hz - freq) / (a.broad_frequency[1] - a.broad_frequency[0])
                        * rect.height())
                    .abs()
                        < 10.
            });
            if let Some(i) = hit {
                drop(a);
                {
                    let mut c = st.borrow_mut();
                    c.points.remove(i);
                    c.pending = None;
                    c.editing = None;
                }
                table(&st);
                return;
            }

            drop(a);
            input("cal-frequency").set_value(&format!("{:.6}", freq / 1e6));
            text(
                "cal-selection",
                &format!(
                    "Broad gate {} · {}–{} UTC · {:.3} MHz (25 kHz bin)",
                    broad + 1,
                    utc(start),
                    utc(end),
                    freq / 1e6
                ),
            );
            text("cal-result", "Loading matching one-minute fits…");
            st.borrow_mut().request += 1;
            let request = st.borrow().request;
            st.borrow_mut().pending = None;
            st.borrow_mut().editing = None;
            el("cal-altitude").set_inner_html("<option value=''>Loading…</option>");
            let st = st.clone();
            let app = app.clone();
            spawn_local(async move {
                let result: Result<(), JsValue> = async {
                    if st.borrow().settings.is_none() {
                        let settings: Settings = json("calibration-config.json").await?;
                        st.borrow_mut().settings = Some(settings);
                    }
                    if st.borrow().response.is_none() {
                        let response: Response = json("data/beata-short-lag-response.json").await?;
                        st.borrow_mut().response = Some(response);
                    }
                    if st.borrow().fits.is_none() {
                        let url = st.borrow().settings.as_ref().unwrap().ion_index_url.clone();
                        let fits: Fits = json(&url).await?;
                        st.borrow_mut().fits = Some(fits);
                    }
                    if st.borrow().request != request {
                        return Ok(());
                    }
                    let profile = st
                        .borrow()
                        .fits
                        .as_ref()
                        .unwrap()
                        .profiles
                        .iter()
                        .find(|p| {
                            p.start_unix <= t
                                && p.end_unix > t
                                && (p.end_unix - p.start_unix - 60.).abs() < 0.01
                        })
                        .cloned()
                        .ok_or_else(|| error("No one-minute fit covers the selected time"))?;
                    let range = [[47., 221.], [128., 302.], [209., 383.], [290., 464.]][broad];
                    let bounds = range.map(|r| altitude(r, elevation));
                    st.borrow_mut().pending = Some(Pending {
                        time: t,
                        start,
                        end,
                        broad,
                        frequency: freq,
                        profile,
                    });
                    candidates(&st.borrow(), bounds);
                    let peak = peak_height(&st.borrow(), bounds).ok_or_else(|| {
                        error("No finite positive raw Ne within this gate and time window")
                    })?;
                    set_height(&st.borrow(), peak);
                    select(app, selected).await?;
                    let candidate = point(&st.borrow());
                    if let Ok(mut p) = candidate {
                        p.included = false; // Peak association remains reviewable before inclusion.
                        let mut c = st.borrow_mut();
                        let i = c.points.len();
                        c.points.push(p);
                        c.editing = Some(i);
                    }
                    table(&st);
                    preview(&st.borrow());
                    Ok(())
                }
                .await;
                if let Err(e) = result {
                    if st.borrow().request == request {
                        text("cal-result", &format!("Calibration unavailable: {e:?}"));
                    }
                }
            });
        })?;
    }
    Ok(())
}

fn set_height(c: &Calibration, h: f64) {
    let Some(p) = &c.pending else { return };
    let lo = input("cal-height").min().parse::<f64>().unwrap_or(40.);
    let hi = input("cal-height").max().parse::<f64>().unwrap_or(470.);
    let h = h.clamp(lo, hi);
    input("cal-height").set_value(&h.to_string());
    if let Some((i, _)) = p
        .profile
        .gates
        .iter()
        .enumerate()
        .filter(|(_, g)| g.result.is_some() && g.altitude_km >= lo && g.altitude_km <= hi)
        .min_by(|(_, a), (_, b)| {
            (a.altitude_km - h)
                .abs()
                .total_cmp(&(b.altitude_km - h).abs())
        })
    {
        el("cal-altitude")
            .dyn_into::<HtmlSelectElement>()
            .unwrap()
            .set_value(&i.to_string());
    }
}
pub(super) fn markers(a: &App, c: &Calibration) {
    let mut layers = vec![String::new(); 5];
    let mut profile = String::new();
    let f = &a.m.frames[a.selected];

    for (i, saved) in c.points.iter().enumerate() {
        let p = saved;
        let color = if c.editing == Some(i) {
            "#236b8e"
        } else {
            "#986327"
        };
        let tx = (p.time_unix - a.view_start) / (a.view_end - a.view_start);
        let fy = (a.broad_frequency[1] + 12500. - p.frequency_hz)
            / (a.broad_frequency[1] - a.broad_frequency[0] + 25000.);
        let ay = (a.profile_altitude[1] - p.altitude_km)
            / (a.profile_altitude[1] - a.profile_altitude[0]);
        let label = format!("{} · {:.3} MHz", i + 1, p.frequency_hz / 1e6);
        let dot = |x: f64, y: f64| {
            format!("<circle cx='{x}' cy='{y}' r='6' fill='#ffffff' stroke='{color}' stroke-width='2'/><text x='{}' y='{}' fill='{color}' font-size='12'>{label}</text>",x+9.,y-9.)
        };
        if (0. ..=1.).contains(&tx) {
            let b = p.broad_gate - 1;
            if (0. ..=1.).contains(&fy) {
                let r = el(&format!("broad{b}")).get_bounding_client_rect();
                layers[b].push_str(&dot(tx * r.width(), fy * r.height()));
            }
            if (0. ..=1.).contains(&ay) {
                let r = el("power").get_bounding_client_rect();
                layers[4].push_str(&dot(tx * r.width(), ay * r.height()));
            }
        }
        if p.plasma_start_unix < f.end_unix
            && p.plasma_end_unix > f.start_unix
            && (0. ..=1.).contains(&ay)
        {
            let rect = el("profile")
                .query_selector("canvas")
                .ok()
                .flatten()
                .map(|c| c.get_bounding_client_rect());
            let width = rect.as_ref().map_or(380., |r| r.width());
            let height = rect.as_ref().map_or(380., |r| r.height());
            let px = 76.
                + (width - 100.) * (p.raw_ne_proxy_m3 / a.profile_power[0]).log10()
                    / (a.profile_power[1] / a.profile_power[0]).log10();
            if (76. ..=width - 24.).contains(&px) {
                profile.push_str(&format!("<span style='position:absolute;left:{px}px;top:{}px;transform:translate(-50%,-50%);color:#986327;font-size:20px'>○<small>{}</small></span>",22.+(height-80.)*ay,i+1));
                let expected = p.raw_ne_proxy_m3 * p.estimate.calibration_multiplier_ratio;
                let ex = 76.
                    + (width - 100.) * (expected / a.profile_power[0]).log10()
                        / (a.profile_power[1] / a.profile_power[0]).log10();
                if (76. ..=width - 24.).contains(&ex) {
                    profile.push_str(&format!("<span title='Predicted short-lag power' style='position:absolute;left:{ex}px;top:{}px;width:8px;height:8px;transform:translate(-50%,-50%);border:1px solid #986327'></span>",22.+(height-80.)*ay));
                }
            }
        }
    }
    for (i, layer) in layers.iter().enumerate() {
        let id = if i == 4 {
            "power".into()
        } else {
            format!("broad{i}")
        };
        let r = el(&id).get_bounding_client_rect();
        let target = if i == 4 {
            "cal-power".into()
        } else {
            format!("cal-broad-{i}")
        };
        el(&target).set_inner_html(&format!(
            "<svg width='100%' height='100%' viewBox='0 0 {} {}'>{layer}</svg>",
            r.width(),
            r.height()
        ));
    }
    let host = el("profile");
    if doc().get_element_by_id("cal-profile-marks").is_none() {
        let layer = doc().create_element("div").unwrap();
        layer.set_id("cal-profile-marks");
        let _ = layer.set_attribute("style", "position:absolute;inset:0;pointer-events:none");
        let _ = host.append_child(&layer);
    }
    el("cal-profile-marks").set_inner_html(&profile);
}
#[derive(Clone, Copy)]
struct Drag {
    x: f64,
    y: f64,
    px: i32,
    py: i32,
    point: Option<usize>,
    moved: bool,
}
fn coords(id: &str, e: &web_sys::PointerEvent) -> (f64, f64) {
    let target = if id == "profile" {
        el(id).query_selector("canvas").unwrap().unwrap()
    } else {
        el(id)
    };
    let r = target.get_bounding_client_rect();
    let mut x = (e.client_x() as f64 - r.left()) / r.width();
    let mut y = (e.client_y() as f64 - r.top()) / r.height();
    if id == "profile" {
        x = (x * r.width() - 76.) / (r.width() - 100.);
        y = (y * r.height() - 22.) / (r.height() - 80.);
    }
    (x.clamp(0., 1.), y.clamp(0., 1.))
}
fn gestures(state: &Cal) -> Result<(), JsValue> {
    for id in ["broad0", "broad1", "broad2", "broad3", "power", "profile"] {
        let _ = el(id).set_attribute(
            "style",
            "touch-action:none;cursor:crosshair;position:relative",
        );
        let drag: Rc<RefCell<Option<Drag>>> = Rc::new(RefCell::new(None));
        let st = state.clone();
        let d = drag.clone();
        on(id, "pointerdown", move |event| {
            let Ok(e) = event.dyn_into::<web_sys::PointerEvent>() else {
                return;
            };
            if e.button() != 0 {
                return;
            }
            st.borrow_mut().suppress_click = false;
            e.prevent_default();
            let (x, y) = coords(id, &e);
            let c = st.borrow();
            let a = c.app.borrow();
            let r = el(id).get_bounding_client_rect();
            let hit = c
                .points
                .iter()
                .enumerate()
                .find(|(_, p)| {
                    let (xx, yy) = if id == "profile" {
                        let f = &a.m.frames[a.selected];
                        if p.plasma_end_unix <= f.start_unix || p.plasma_start_unix >= f.end_unix {
                            return false;
                        }
                        (
                            (p.raw_ne_proxy_m3 / a.profile_power[0]).log10()
                                / (a.profile_power[1] / a.profile_power[0]).log10(),
                            (a.profile_altitude[1] - p.altitude_km)
                                / (a.profile_altitude[1] - a.profile_altitude[0]),
                        )
                    } else if id == "power" {
                        (
                            (p.time_unix - a.view_start) / (a.view_end - a.view_start),
                            (a.profile_altitude[1] - p.altitude_km)
                                / (a.profile_altitude[1] - a.profile_altitude[0]),
                        )
                    } else {
                        return false;
                    };
                    ((xx - x) * r.width()).abs() < 14. && ((yy - y) * r.height()).abs() < 14.
                })
                .map(|(i, _)| i);
            drop(a);
            drop(c);
            if let Some(i) = hit {
                el(&format!("cal-edit-{i}"))
                    .dyn_into::<web_sys::HtmlElement>()
                    .unwrap()
                    .click();
            }
            *d.borrow_mut() = Some(Drag {
                x,
                y,
                px: e.client_x(),
                py: e.client_y(),
                point: hit,
                moved: false,
            });
            let _ = el(id).set_pointer_capture(e.pointer_id());
        })?;
        let st = state.clone();
        let d = drag.clone();
        on(id, "pointermove", move |event| {
            let Ok(e) = event.dyn_into::<web_sys::PointerEvent>() else {
                return;
            };
            let Some(mut drag) = *d.borrow() else { return };
            drag.moved |= (e.client_x() - drag.px)
                .abs()
                .max((e.client_y() - drag.py).abs())
                > 5;
            *d.borrow_mut() = Some(drag);
            if !drag.moved {
                return;
            }
            let (x, y) = coords(id, &e);
            if drag.point.is_some() {
                let c = st.borrow();
                let bounds = c.app.borrow().profile_altitude;
                set_height(&c, bounds[1] - y * (bounds[1] - bounds[0]));
                let p = point(&c);
                drop(c);
                if let (Some(i), Ok(p)) = (drag.point, p) {
                    st.borrow_mut().points[i] = p;
                }
                preview(&st.borrow());
                return;
            }
            let r = el(id).get_bounding_client_rect();
            let (left, top, w, h) = if id == "profile" {
                (
                    (76. + (r.width() - 100.) * x.min(drag.x)) / r.width(),
                    (22. + (r.height() - 80.) * y.min(drag.y)) / r.height(),
                    (r.width() - 100.) * (x - drag.x).abs() / r.width(),
                    (r.height() - 80.) * (y - drag.y).abs() / r.height(),
                )
            } else {
                (
                    x.min(drag.x),
                    y.min(drag.y),
                    (x - drag.x).abs(),
                    (y - drag.y).abs(),
                )
            };
            let overlay = if let Some(v) = doc().get_element_by_id("cal-brush") {
                v
            } else {
                let v = doc().create_element("div").unwrap();
                let _ = v.set_attribute("id", "cal-brush");
                let _ = doc().body().unwrap().append_child(&v);
                v
            };
            let _=overlay.set_attribute("style",&format!("position:fixed;pointer-events:none;z-index:9999;border:1px solid #a9fff3;background:#a9fff333;left:{}px;top:{}px;width:{}px;height:{}px",r.left()+left*r.width(),r.top()+top*r.height(),w*r.width(),h*r.height()));
        })?;
        let st = state.clone();
        let d = drag.clone();
        on(id, "pointerup", move |event| {
            let Ok(e) = event.dyn_into::<web_sys::PointerEvent>() else {
                return;
            };
            let Some(drag) = d.borrow_mut().take() else {
                return;
            };
            let _ = el(id).release_pointer_capture(e.pointer_id());
            if let Some(b) = doc().get_element_by_id("cal-brush") {
                b.remove();
            }
            let (x, y) = coords(id, &e);
            if let Some(i) = drag.point {
                if drag.moved {
                    let p = point(&st.borrow());
                    if let Ok(p) = p {
                        st.borrow_mut().points[i] = p;
                    }
                } else {
                    let mut c = st.borrow_mut();
                    c.points.remove(i);
                    c.pending = None;
                    c.editing = None;
                }
                table(&st);
                return;
            }
            if !drag.moved {
                if id == "power" || id == "profile" {
                    if id == "power" {
                        let mut c = st.borrow_mut();
                        let app = c.app.clone();
                        let a = app.borrow();
                        let t = a.view_start + x * (a.view_end - a.view_start);
                        let frame =
                            a.m.frames
                                .iter()
                                .find(|f| f.start_unix <= t && f.end_unix > t);
                        let fit = c
                            .fits
                            .as_ref()
                            .and_then(|fs| {
                                fs.profiles
                                    .iter()
                                    .find(|p| p.start_unix <= t && p.end_unix > t)
                            })
                            .cloned();
                        if let (Some(f), Some(fit), Some(p)) = (frame, fit, c.pending.as_mut()) {
                            p.time = (f.start_unix + f.end_unix) / 2.;
                            p.start = f.start_unix;
                            p.end = f.end_unix;
                            p.profile = fit;
                            c.editing = None;
                        } else {
                            text("cal-result","Select a plasma-line frequency in a broad gate first; a matching one-minute fit is required.");
                            return;
                        }
                    }
                    let c = st.borrow();
                    let bounds = c.app.borrow().profile_altitude;
                    set_height(&c, bounds[1] - y * (bounds[1] - bounds[0]));
                    let p = point(&c);
                    drop(c);
                    if let Ok(p) = p {
                        let mut c = st.borrow_mut();
                        if let Some(i) = c.editing {
                            c.points[i] = p;
                        } else {
                            let i = c.points.len();
                            c.points.push(p);
                            c.editing = Some(i);
                        }
                        drop(c);
                        table(&st);
                    }
                }
                return;
            }
            st.borrow_mut().suppress_click = id.starts_with("broad");
            let app = st.borrow().app.clone();
            let mut a = app.borrow_mut();
            if (x - drag.x).abs() > 0.005 && (y - drag.y).abs() > 0.005 {
                if id == "profile" {
                    let [lo, hi] = a.profile_power;
                    a.profile_power = [
                        lo * (hi / lo).powf(x.min(drag.x)),
                        lo * (hi / lo).powf(x.max(drag.x)),
                    ];
                } else {
                    let lo = a.view_start;
                    let span = a.view_end - lo;
                    a.view_start = lo + x.min(drag.x) * span;
                    a.view_end = lo + x.max(drag.x) * span;
                }
                if id.starts_with("broad") {
                    let [lo, hi] = a.broad_frequency;
                    a.broad_frequency = [
                        hi - y.max(drag.y) * (hi - lo),
                        hi - y.min(drag.y) * (hi - lo),
                    ];
                } else {
                    let [lo, hi] = a.profile_altitude;
                    let lower = hi - y.max(drag.y) * (hi - lo);
                    let upper = hi - y.min(drag.y) * (hi - lo);
                    if upper - lower >= 3. {
                        a.profile_altitude = [lower, upper];
                    }
                }
            }
            let _ = overview(&a);
        })?;
        let d = drag.clone();
        on(id, "pointercancel", move |_| {
            *d.borrow_mut() = None;
            if let Some(b) = doc().get_element_by_id("cal-brush") {
                b.remove();
            }
        })?;
        let st = state.clone();
        on(id, "contextmenu", move |e| {
            e.prevent_default();
            let app = st.borrow().app.clone();
            let mut a = app.borrow_mut();
            if id != "profile" {
                a.view_start = a.m.frames[0].start_unix;
                a.view_end = a.m.frames.last().unwrap().end_unix;
            }
            if id.starts_with("broad") {
                a.broad_frequency = [-6.15e6, -1.275e6];
            } else {
                a.profile_altitude = [40., 400.];
                a.profile_power = [1e8, 1e13];
            }
            let _ = overview(&a);
        })?;
    }
    Ok(())
}

#[derive(Deserialize)]
struct Saved {
    revision: u64,
    points: Vec<Point>,
}
struct Sync {
    key: String,
    revision: u64,
    ready: bool,
    busy: bool,
    pending: Option<serde_json::Value>,
}
fn flush(state: Rc<RefCell<Sync>>) {
    {
        let mut s = state.borrow_mut();
        if !s.ready || s.busy {
            return;
        }
        s.busy = true;
    }
    spawn_local(async move {
        loop {
            let (key, revision, points) = {
                let mut s = state.borrow_mut();
                let Some(p) = s.pending.take() else {
                    s.busy = false;
                    break;
                };
                (s.key.clone(), s.revision, p)
            };
            text("cal-save-status", "Saving calibration to server…");
            let result: Result<Saved, JsValue> = async {
                let init = web_sys::RequestInit::new();
                init.set_method("PUT");
                let body = serde_json::json!({"revision":revision,"points":points}).to_string();
                init.set_body(&JsValue::from_str(&body));
                let headers = js_sys::Object::new();
                js_sys::Reflect::set(&headers, &"Content-Type".into(), &"application/json".into())?;
                init.set_headers(&headers);
                let response = JsFuture::from(
                    web_sys::window()
                        .unwrap()
                        .fetch_with_str_and_init(&format!("api/calibrations/{key}"), &init),
                )
                .await?
                .dyn_into::<web_sys::Response>()?;
                if !response.ok() {
                    return Err(error(format!(
                        "HTTP {} (reload before editing if another tab changed this notebook)",
                        response.status()
                    )));
                }
                let value = JsFuture::from(response.text()?).await?.as_string().unwrap();
                serde_json::from_str(&value).map_err(|e| error(e.to_string()))
            }
            .await;
            match result {
                Ok(n) => {
                    state.borrow_mut().revision = n.revision;
                    text(
                        "cal-save-status",
                        &format!(
                            "Saved on avaruus · revision {} · browser copy retained",
                            n.revision
                        ),
                    );
                }
                Err(e) => {
                    let mut s = state.borrow_mut();
                    if s.pending.is_none() {
                        s.pending = Some(points);
                    }
                    s.busy = false;
                    text(
                        "cal-save-status",
                        &format!("Not saved on server; browser copy retained. {e:?}"),
                    );
                    break;
                }
            }
        }
    });
}
fn update_current(state: &Cal) {
    let p = point(&state.borrow());
    if let Ok(p) = p {
        let mut c = state.borrow_mut();
        if let Some(i) = c.editing {
            c.points[i] = p;
        }
    }
    table(state);
    preview(&state.borrow());
}
fn peak_height(c: &Calibration, bounds: [f64; 2]) -> Option<f64> {
    let p = c.pending.as_ref()?;
    let a = c.app.borrow();
    let width = input("cal-time-width").value().parse::<f64>().ok()?;
    let n = a.m.frames.len();
    a.m.ion_range_km
        .iter()
        .enumerate()
        .filter_map(|(g, &r)| {
            let mut sum = 0.;
            let mut weight = 0.;
            let mut height = 0.;
            for (t, f) in a.m.frames.iter().enumerate() {
                let mid = (f.start_unix + f.end_unix) / 2.;
                let h = altitude(r, f.elevation_deg);
                let v = a.power[g * n + t] as f64;
                if mid >= p.time - width / 2.
                    && mid < p.time + width / 2.
                    && h >= bounds[0]
                    && h <= bounds[1]
                    && v.is_finite()
                {
                    let w = f.end_unix - f.start_unix;
                    sum += v * w;
                    height += h * w;
                    weight += w;
                }
            }
            (weight > 0. && sum > 0.).then(|| (height / weight, sum / weight))
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .map(|p| p.0)
}

/// Freeze the exact saved notebook revision into an authenticated queue job.
fn recalculation(state: &Cal) -> Result<(), JsValue> {
    let state = state.clone();
    let callback = Closure::<dyn FnMut()>::new(move || {
        let request = (|| -> Result<(String, String, u64), JsValue> {
            let c = state.borrow();
            let sync = c.sync.borrow();
            if !sync.ready || sync.busy || sync.pending.is_some() {
                return Err(error(
                    "Wait for the calibration revision to finish saving before recalculating.",
                ));
            }
            if !c.points.iter().any(|p| p.included) {
                return Err(error("Include reviewed calibration points first."));
            }
            let query = web_sys::window().unwrap().location().search()?;
            let id = query
                .trim_start_matches('?')
                .split('&')
                .find_map(|s| s.strip_prefix("experiment=").map(str::to_owned))
                .or_else(|| c.settings.as_ref().and_then(|s| s.experiment_id.clone()))
                .ok_or_else(|| error("No processing recipe is associated with this observation"))?;
            Ok((id, sync.key.clone(), sync.revision))
        })();
        let (id, key, revision) = match request {
            Ok(v) => v,
            Err(e) => {
                text("cal-job-status", &e.as_string().unwrap_or_default());
                return;
            }
        };
        let token = input("cal-processing-key").value();
        if token.trim().is_empty() {
            text(
                "cal-job-status",
                "Enter the processing key under Processing access.",
            );
            return;
        }
        text("cal-job-status", "Queuing the saved calibration revision…");
        spawn_local(async move {
            let result = async {
                let init = web_sys::RequestInit::new();
                init.set_method("POST");
                let headers = web_sys::Headers::new()?;
                headers.set("Authorization", &format!("Bearer {}", token.trim()))?;
                headers.set("Content-Type", "application/json")?;
                init.set_headers(&headers);
                init.set_body(
                    &serde_json::json!({"notebook_key":key,"revision":revision})
                        .to_string()
                        .into(),
                );
                let response: web_sys::Response = JsFuture::from(
                    web_sys::window()
                        .unwrap()
                        .fetch_with_str_and_init(&format!("api/catalog/jobs/{id}"), &init),
                )
                .await?
                .dyn_into()?;
                let body = JsFuture::from(response.text()?)
                    .await?
                    .as_string()
                    .unwrap_or_default();
                if !response.ok() {
                    return Err(error(body));
                }
                Ok(())
            }
            .await;
            match result {
                Ok(()) => {
                    el("cal-job-status").set_inner_html("Ion-line recalculation queued. <a href='live.html'>View processing status</a>.");
                }
                Err(e) => text("cal-job-status", &e.as_string().unwrap_or_default()),
            }
        });
    });
    el("cal-recalculate")
        .add_event_listener_with_callback("click", callback.as_ref().unchecked_ref())?;
    callback.forget();
    Ok(())
}
