//! Instrument-neutral observation browser over versioned HDF5 static exports.
use crate::{el, error, json, plotting, text, utc};
use serde::{Deserialize, Serialize};
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{prelude::*, JsCast};
use wasm_bindgen_futures::{spawn_local, JsFuture};
use web_sys::{HtmlInputElement, HtmlSelectElement};
#[derive(Deserialize)]
struct Config {
    indexes: Vec<String>,
}
#[derive(Deserialize)]
struct Index {
    schema: String,
    records: Vec<Entry>,
}
#[derive(Clone, Deserialize)]
struct Entry {
    id: String,
    station: String,
    antenna: String,
    mode: String,
    kind: String,
    start_unix_s: f64,
    end_unix_s: f64,
    url: String,
    #[serde(default)]
    has_spectra: bool,
    #[serde(default)]
    has_fit: bool,
    #[serde(skip)]
    root: String,
}
#[derive(Clone, Deserialize)]
struct Spectrum {
    name: String,
    url: String,
    shape: [usize; 2],
    frequency_hz: Vec<f64>,
    range_m: Vec<f64>,
    units: String,
    frequency_average_bins: usize,
}
#[derive(Clone, Deserialize)]
struct Fit {
    parameters: Vec<[Option<f64>; 4]>,
    standard_error: Vec<[Option<f64>; 4]>,
    #[serde(rename = "P_tx")]
    transmitter_w: Option<f64>,
    #[serde(rename = "T_sys")]
    system_temperature_k: Option<f64>,
    #[serde(default)]
    space_object_count: Vec<f64>,
}
#[derive(Clone, Deserialize)]
struct Detail {
    id: String,
    station: String,
    antenna: String,
    mode: String,
    radar_frequency_hz: f64,
    start_unix_s: f64,
    end_unix_s: f64,
    #[serde(default)]
    range_m: Vec<f64>,
    #[serde(default)]
    altitude_m: Vec<f64>,
    density_semantics: Option<String>,
    #[serde(default)]
    spectra: Vec<Spectrum>,
    fit: Option<Fit>,
    #[serde(default)]
    raw_density: Vec<Option<f64>>,
}
#[derive(Clone, Serialize, Deserialize)]
struct Point {
    observation_id: String,
    fit_id: String,
    station: String,
    antenna: String,
    mode: String,
    frequency_hz: f64,
    frequency_standard_error_hz: f64,
    range_m: f64,
    gate: usize,
    plasma_density_m3: f64,
    k2_debye2: f64,
    multiplier: f64,
    included: bool,
    flags: Vec<String>,
}
struct View {
    entries: Vec<Entry>,
    selected: usize,
    band: usize,
    detail: Option<Detail>,
    fit: Option<Detail>,
    values: Vec<f32>,
    points: Vec<Point>,
    request: u64,
}
type State = Rc<RefCell<View>>;
fn input(id: &str) -> f64 {
    el(id)
        .dyn_into::<HtmlInputElement>()
        .unwrap()
        .value_as_number()
}
fn selector(id: &str) -> String {
    el(id).dyn_into::<HtmlSelectElement>().unwrap().value()
}
fn label(e: &Entry) -> String {
    format!(
        "{} · {} · {} · {}–{} UTC",
        e.station,
        e.antenna,
        e.mode,
        utc(e.start_unix_s),
        utc(e.end_unix_s)
    )
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn options(id: &str, items: impl Iterator<Item = (usize, String)>) {
    el(id).set_inner_html(
        &items
            .map(|(i, s)| format!("<option value='{i}'>{}</option>", escape(&s)))
            .collect::<String>(),
    );
}
async fn load(state: State, index: usize) -> Result<(), JsValue> {
    let (entry, request) = {
        let mut v = state.borrow_mut();
        v.request += 1;
        v.selected = index;
        v.band = 0;
        (v.entries[index].clone(), v.request)
    };
    text("status", "Loading integration…");
    let detail: Detail = json(&format!("{}{}", entry.root, entry.url)).await?;
    let candidate = {
        let v = state.borrow();
        v.entries
            .iter()
            .enumerate()
            .filter(|(_, e)| {
                e.has_fit
                    && e.station == entry.station
                    && e.antenna == entry.antenna
                    && ((e.start_unix_s + e.end_unix_s - entry.start_unix_s - entry.end_unix_s)
                        / 2.)
                        .abs()
                        <= 300.
            })
            .min_by(|(_, a), (_, b)| {
                ((a.start_unix_s + a.end_unix_s - entry.start_unix_s - entry.end_unix_s) / 2.)
                    .abs()
                    .total_cmp(
                        &((b.start_unix_s + b.end_unix_s - entry.start_unix_s - entry.end_unix_s)
                            / 2.)
                            .abs(),
                    )
            })
            .map(|(i, e)| (i, e.clone()))
    };
    let counterpart = if let Some((_, e)) = &candidate {
        Some(json::<Detail>(&format!("{}{}", e.root, e.url)).await?)
    } else {
        None
    };
    {
        let mut v = state.borrow_mut();
        if request != v.request {
            return Ok(());
        }
        v.detail = Some(detail.clone());
        v.fit = counterpart;
        options(
            "observation-fit",
            std::iter::once((usize::MAX, "No counterpart".into())).chain(
                v.entries
                    .iter()
                    .enumerate()
                    .filter(|(_, e)| {
                        e.has_fit && e.station == entry.station && e.antenna == entry.antenna
                    })
                    .map(|(i, e)| (i, label(e))),
            ),
        );
        if let Some((i, _)) = candidate {
            el("observation-fit")
                .dyn_into::<HtmlSelectElement>()?
                .set_value(&i.to_string());
        }
        options(
            "observation-band",
            detail
                .spectra
                .iter()
                .enumerate()
                .map(|(i, s)| (i, s.name.clone())),
        );
        text("observation-context", &label(&entry));
        let query = web_sys::window().unwrap().location().search()?;
        let ion = query
            .trim_start_matches('?')
            .split('&')
            .any(|s| s == "view=ion");
        for (kind, id) in [
            ("plasma", "observation-plasma-tab"),
            ("ion", "observation-ion-tab"),
        ] {
            el(id).set_attribute("href", &format!("?record={}&view={kind}", entry.id))?;
            if (kind == "ion") == ion {
                el(id).set_attribute("aria-current", "page")?;
            } else {
                el(id).remove_attribute("aria-current")?;
            }
        }
        if ion {
            text("observation-heading", "Ion-line profiles");
            text("observation-help", "Select an integration to inspect its fitted density profile. A trial multiplier changes the displayed density only.");
            el("observation-band-control").set_attribute("hidden", "")?;
            el("observation-display-control").set_attribute("hidden", "")?;
            el("observation-calibration").set_attribute("hidden", "")?;
            el("observation-plasma-panel")
                .parent_element()
                .unwrap()
                .set_class_name("observation-layout ion-only");
        }
    }
    load_band(state).await
}
async fn load_band(state: State) -> Result<(), JsValue> {
    let (s, root, request) = {
        let v = state.borrow();
        (
            v.detail
                .as_ref()
                .and_then(|d| d.spectra.get(v.band))
                .cloned(),
            v.entries[v.selected].root.clone(),
            v.request,
        )
    };
    if let Some(s) = s {
        let r = crate::response(&format!("{root}{}", s.url)).await?;
        let b = JsFuture::from(r.array_buffer()?).await?;
        let bytes = js_sys::Uint8Array::new(&b).to_vec();
        if bytes.len() != s.shape[0] * s.shape[1] * 4 {
            return Err(error("spectrum shape mismatch"));
        }
        let values = bytes
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
            .collect();
        let mut v = state.borrow_mut();
        if request != v.request {
            return Ok(());
        }
        v.values = values;
    } else {
        state.borrow_mut().values.clear();
    }
    render(&state.borrow())?;
    text("status", "Observation loaded");
    Ok(())
}
fn render(v: &View) -> Result<(), JsValue> {
    if let Some(s) = v.detail.as_ref().and_then(|d| d.spectra.get(v.band)) {
        let mut values = v.values.clone();
        if selector("observation-baseline") == "contrast" {
            for row in values.chunks_exact_mut(s.shape[1]) {
                let mut sorted: Vec<_> = row.iter().copied().filter(|v| v.is_finite()).collect();
                sorted.sort_by(f32::total_cmp);
                let median = sorted.get(sorted.len() / 2).copied().unwrap_or(0.);
                for v in row {
                    *v -= median
                }
            }
        }
        let limit = crate::percentile(&values, 0.995);
        plotting::rectilinear("observation-spectrum", &values, &s.frequency_hz, &s.range_m)?;
        text(
            "observation-spectrum-units",
            &format!(
                "{} · {} source bins averaged",
                s.units, s.frequency_average_bins
            ),
        );
        text(
            "observation-range-top",
            &format!(
                "{:.1}–{:.1}",
                s.range_m[0] / 1000.,
                s.range_m.last().unwrap() / 1000.
            ),
        );
        text(
            "observation-frequency-low",
            &format!("{:.3}", s.frequency_hz[0] / 1e6),
        );
        text(
            "observation-frequency-high",
            &format!("{:.3}", s.frequency_hz.last().unwrap() / 1e6),
        );
        text("observation-scale-low", &format!("−{limit:.3e}"));
        text("observation-scale-high", &format!("+{limit:.3e}"));
    } else {
        text(
            "observation-spectrum-units",
            "No plasma spectrum in this record",
        );
        plotting::shared_field(
            "observation-spectrum",
            &[f32::NAN],
            1,
            1,
            &serde_json::json!({"bounds":[0.,1.,0.,1.],"limits":[0.,1.],"time":false,"xlabel":"Frequency offset / MHz","ylabel":"Slant range / km"}),
        )?;
    }
    profile(v);
    points(v);
    Ok(())
}
fn profile(v: &View) {
    let d = v.fit.as_ref().or(v.detail.as_ref());
    let Some(d) = d else { return };
    let values: Vec<Option<f64>> = if let Some(f) = &d.fit {
        f.parameters.iter().map(|p| p[0]).collect()
    } else {
        d.raw_density.clone()
    };
    if values.is_empty() {
        text("observation-profile", "No ion-line counterpart available");
        text(
            "observation-fit-context",
            "Select a fit integration from the same antenna.",
        );
        return;
    }
    let gain = input("observation-multiplier");
    if !gain.is_finite() || gain <= 0. {
        text(
            "observation-fit-context",
            "Trial multiplier must be finite and positive.",
        );
        return;
    }
    let alts = if d.altitude_m.len() == values.len() {
        &d.altitude_m
    } else {
        &d.range_m
    };
    if alts.len() != values.len() {
        return;
    }
    let mut logs: Vec<_> = values
        .iter()
        .flatten()
        .filter(|x| **x > 0.)
        .map(|x| (x * gain).log10())
        .collect();
    logs.sort_by(f64::total_cmp);
    let lo = logs.first().copied().unwrap_or(0.).floor();
    let hi = logs.last().copied().unwrap_or(1.).ceil().max(lo + 1.);
    let bottom = alts[0] / 1000.;
    let top = alts.last().unwrap() / 1000.;
    let points: Vec<_> = values
        .iter()
        .zip(alts)
        .flat_map(|(n, h)| {
            [
                n.filter(|n| *n > 0.).map(|n| (n * gain).log10()),
                Some(h / 1000.),
            ]
        })
        .collect();
    plotting::shared_curve("observation-profile",&serde_json::json!({"title":"Ion-line density profile","height":380,"bounds":[lo,hi,bottom,top],"xlabel":if d.density_semantics.as_deref()==Some("uncalibrated_ne_proxy"){"Density proxy"}else{"Ne / m⁻³"},"xformat":"log10","ylabel":if d.altitude_m.len()==values.len(){"Altitude / km"}else{"Slant range / km"},"series":[{"xy":points,"radius":2.}]}).to_string());
    text(
        "observation-fit-context",
        &format!(
            "{} · {} · {}–{} UTC. {}",
            d.antenna,
            d.mode,
            utc(d.start_unix_s),
            utc(d.end_unix_s),
            if let Some(f) = &d.fit {
                format!("Ptx {:.0} kW; Tsys {:.0} K. Source uncertainty retained; temperature covariance unavailable.",f.transmitter_w.unwrap_or(f64::NAN)/1000.,f.system_temperature_k.unwrap_or(f64::NAN))
            } else {
                "Signed raw measurements; nonpositive values are outside the logarithmic display."
                    .into()
            }
        ),
    );
}
fn candidate(v: &View) -> Result<Point, JsValue> {
    let d = v
        .detail
        .as_ref()
        .ok_or_else(|| error("Select an observation"))?;
    let fit = v
        .fit
        .as_ref()
        .ok_or_else(|| error("Select a matching ion-line fit"))?;
    let f = fit
        .fit
        .as_ref()
        .ok_or_else(|| error("Ion fit unavailable"))?;
    if fit.station != d.station || fit.antenna != d.antenna {
        return Err(error("Antenna mismatch"));
    }
    let dt = ((fit.start_unix_s + fit.end_unix_s - d.start_unix_s - d.end_unix_s) / 2.).abs();
    if dt > 300. {
        return Err(error(
            "Ion-line and plasma-line integration centres differ by more than 300 s",
        ));
    }
    let frequency = input("observation-frequency") * 1e6;
    let range = input("observation-range") * 1000.;
    let sigma = input("observation-frequency-sigma") * 1000.;
    if ![frequency, range, sigma].iter().all(|v| v.is_finite()) || sigma <= 0. {
        return Err(error(
            "Provide frequency, range, and positive frequency uncertainty",
        ));
    }
    let s = d
        .spectra
        .get(v.band)
        .ok_or_else(|| error("No plasma spectrum selected"))?;
    if frequency < s.frequency_hz[0]
        || frequency > *s.frequency_hz.last().unwrap()
        || range < s.range_m[0]
        || range > *s.range_m.last().unwrap()
    {
        return Err(error("Selection is outside measured plasma coordinates"));
    }
    let (gate, _) = fit
        .range_m
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| (*a - range).abs().total_cmp(&(*b - range).abs()))
        .ok_or_else(|| error("Fit range missing"))?;
    if range < fit.range_m[0] || range > *fit.range_m.last().unwrap() {
        return Err(error("Selection is outside measured ion-line range"));
    }
    let p = f.parameters[gate];
    let raw = p[0].ok_or_else(|| error("Density proxy unavailable"))?;
    let ti = p[1].ok_or_else(|| error("Ti unavailable"))?;
    let te = p[2].ok_or_else(|| error("Te unavailable"))?;
    let estimate =
        guisdap_core::calibration::estimate(frequency, d.radar_frequency_hz, te, ti, 1e11)
            .map_err(error)?;
    if fit.density_semantics.as_deref() != Some("uncalibrated_ne_proxy") || raw <= 0. {
        return Err(error(
            "This calibration requires a positive Millstone uncalibrated density proxy",
        ));
    }
    let multiplier =
        estimate.plasma_density_m3 * (1. + te / ti) / estimate.ion_power_denominator / raw;
    let mut flags = Vec::new();
    if estimate.k2_debye2 >= 0.1 {
        flags.push("Bohm–Gross approximation: k²λD² ≥ 0.1".into())
    }
    if dt > 60. {
        flags.push(format!("Integration centres differ by {dt:.0} s"))
    }
    if f.transmitter_w.is_none_or(|x| x < 400e3) {
        flags.push("Low or missing transmit power".into())
    }
    if f.system_temperature_k.is_none_or(|x| x > 2000.) {
        flags.push("High or missing system temperature".into())
    }
    if f.space_object_count.get(gate).is_some_and(|x| *x > 0.) {
        flags.push("Space object in source gate".into())
    }
    if f.standard_error[gate][0].is_none_or(|s| s / raw > 0.5) {
        flags.push("Density proxy uncertainty missing or >50%".into())
    }
    Ok(Point {
        observation_id: d.id.clone(),
        fit_id: fit.id.clone(),
        station: d.station.clone(),
        antenna: d.antenna.clone(),
        mode: fit.mode.clone(),
        frequency_hz: frequency,
        frequency_standard_error_hz: sigma,
        range_m: fit.range_m[gate],
        gate,
        plasma_density_m3: estimate.plasma_density_m3,
        k2_debye2: estimate.k2_debye2,
        multiplier,
        included: false,
        flags,
    })
}
fn points(v: &View) {
    let mut html=String::from("<tr><th>Use</th><th>Antenna / mode</th><th>Frequency / MHz</th><th>Range / km</th><th>Multiplier</th><th>Review</th><th></th></tr>");
    for (i, p) in v.points.iter().enumerate() {
        html.push_str(&format!("<tr><td><input type='checkbox' data-point='{i}' {} aria-label='Include calibration point {}'></td><td>{} / {}</td><td>{:+.4}</td><td>{:.1}</td><td>{:.5e}</td><td>{}</td><td><button data-delete='{i}'>Delete</button></td></tr>",if p.included{"checked"}else{""},i+1,escape(&p.antenna),escape(&p.mode),p.frequency_hz/1e6,p.range_m/1000.,p.multiplier,escape(&p.flags.join("; "))));
    }
    el("observation-points").set_inner_html(&html);
    let current = v.fit.as_ref();
    let mut selected: Vec<_> = v
        .points
        .iter()
        .filter(|p| {
            p.included
                && current.is_some_and(|d| {
                    d.station == p.station && d.antenna == p.antenna && d.mode == p.mode
                })
        })
        .map(|p| p.multiplier)
        .collect();
    selected.sort_by(f64::total_cmp);
    text(
        "observation-gain",
        &if selected.is_empty() {
            "No included points for the selected antenna and fit mode.".into()
        } else {
            format!("Median multiplier: {:.6e} · {} included points for this antenna and mode; conditional uncertainty not yet estimated.",selected[selected.len()/2],selected.len())
        },
    );
    let data = serde_json::json!({"schema":"ruisdap-calibration-review-v1","points":v.points});
    let encoded = js_sys::encode_uri_component(&data.to_string())
        .as_string()
        .unwrap();
    let _ = el("observation-export").set_attribute(
        "href",
        &format!("data:application/json;charset=utf-8,{encoded}"),
    );
    if let Ok(Some(storage)) = web_sys::window().unwrap().local_storage() {
        let _ = storage.set_item("ruisdap-isr-calibration-v1", &data["points"].to_string());
    }
}
#[wasm_bindgen]
pub async fn start_observations() -> Result<(), JsValue> {
    crate::navigation::setup()?;
    let cfg: Config = json("observations-config.json").await?;
    let mut entries = Vec::new();
    for url in cfg.indexes {
        let root = url
            .rsplit_once('/')
            .map(|(p, _)| format!("{p}/"))
            .unwrap_or_default();
        let idx: Index = json(&url).await?;
        if idx.schema != "ruisdap-isr-v1" {
            return Err(error("Unsupported observation schema"));
        }
        for mut e in idx.records {
            e.root = root.clone();
            entries.push(e)
        }
    }
    entries.sort_by(|a, b| a.start_unix_s.total_cmp(&b.start_unix_s));
    if entries.is_empty() {
        return Err(error("No observations exported"));
    }
    let requested = web_sys::window()
        .unwrap()
        .location()
        .search()?
        .trim_start_matches('?')
        .split('&')
        .find_map(|s| s.strip_prefix("record=").map(str::to_owned));
    let initial = requested
        .as_ref()
        .and_then(|id| entries.iter().position(|e| &e.id == id))
        .or_else(|| {
            entries
                .iter()
                .position(|e| e.has_spectra && e.station == "Millstone Hill")
        })
        .unwrap_or(0);
    options(
        "observation-record",
        entries
            .iter()
            .enumerate()
            .map(|(i, e)| (i, format!("{} · {}", label(e), e.kind))),
    );
    el("observation-record")
        .dyn_into::<HtmlSelectElement>()?
        .set_value(&initial.to_string());
    let saved = web_sys::window()
        .unwrap()
        .local_storage()?
        .and_then(|s| s.get_item("ruisdap-isr-calibration-v1").ok().flatten())
        .and_then(|s| serde_json::from_str(&s).ok())
        .unwrap_or_default();
    let state = Rc::new(RefCell::new(View {
        entries,
        selected: initial,
        band: 0,
        detail: None,
        fit: None,
        values: Vec::new(),
        points: saved,
        request: 0,
    }));
    for id in ["observation-record", "observation-band", "observation-fit"] {
        let state = state.clone();
        let callback = Closure::<dyn FnMut()>::new(move || {
            let state = state.clone();
            spawn_local(async move {
                let result = async {
                    match id {
                        "observation-record" => {
                            load(state.clone(), selector(id).parse().unwrap_or(0)).await
                        }
                        "observation-band" => {
                            state.borrow_mut().band = selector(id).parse().unwrap_or(0);
                            load_band(state.clone()).await
                        }
                        _ => {
                            let i = selector(id).parse::<usize>().unwrap_or(usize::MAX);
                            let e = state.borrow().entries.get(i).cloned();
                            let fit = if let Some(e) = e {
                                Some(json::<Detail>(&format!("{}{}", e.root, e.url)).await?)
                            } else {
                                None
                            };
                            state.borrow_mut().fit = fit;
                            render(&state.borrow())
                        }
                    }
                }
                .await;
                if let Err(e) = result {
                    text("status", &e.as_string().unwrap_or_default());
                }
            });
        });
        el(id).add_event_listener_with_callback("change", callback.as_ref().unchecked_ref())?;
        callback.forget();
    }
    for id in ["observation-baseline", "observation-multiplier"] {
        let state = state.clone();
        let callback = Closure::<dyn FnMut()>::new(move || {
            if let Err(e) = render(&state.borrow()) {
                crate::report(e)
            }
        });
        el(id).add_event_listener_with_callback("change", callback.as_ref().unchecked_ref())?;
        callback.forget();
    }
    let view = state.clone();
    let click = Closure::<dyn FnMut(web_sys::MouseEvent)>::new(move |e: web_sys::MouseEvent| {
        let v = view.borrow();
        let Some(s) = v.detail.as_ref().and_then(|d| d.spectra.get(v.band)) else {
            return;
        };
        let rect = el("observation-spectrum").get_bounding_client_rect();
        let fx = (e.client_x() as f64 - rect.x()) / rect.width();
        let fy = 1. - (e.client_y() as f64 - rect.y()) / rect.height();
        let f = s.frequency_hz[0] + fx * (s.frequency_hz.last().unwrap() - s.frequency_hz[0]);
        let r = s.range_m[0] + fy * (s.range_m.last().unwrap() - s.range_m[0]);
        for (id, value) in [
            ("observation-frequency", f / 1e6),
            ("observation-range", r / 1000.),
        ] {
            el(id)
                .dyn_into::<HtmlInputElement>()
                .unwrap()
                .set_value_as_number(value)
        }
        text(
            "observation-pick",
            &match candidate(&v) {
                Ok(p) => format!(
                    "Candidate Ne {:.3e} m⁻³, k²λD² {:.3}; multiplier {:.6e}. {}",
                    p.plasma_density_m3,
                    p.k2_debye2,
                    p.multiplier,
                    p.flags.join("; ")
                ),
                Err(e) => e.as_string().unwrap_or_default(),
            },
        );
    });
    el("observation-spectrum")
        .add_event_listener_with_callback("click", click.as_ref().unchecked_ref())?;
    click.forget();
    let view = state.clone();
    let save = Closure::<dyn FnMut()>::new(move || {
        let p = candidate(&view.borrow());
        match p {
            Ok(p) => {
                view.borrow_mut().points.push(p);
                points(&view.borrow());
                text(
                    "observation-pick",
                    "Point saved excluded. Review its source gate and quality before inclusion.",
                )
            }
            Err(e) => text("observation-pick", &e.as_string().unwrap_or_default()),
        }
    });
    el("observation-save")
        .add_event_listener_with_callback("click", save.as_ref().unchecked_ref())?;
    save.forget();
    let view = state.clone();
    let toggle = Closure::<dyn FnMut(web_sys::Event)>::new(move |e: web_sys::Event| {
        let Some(target) = e
            .target()
            .and_then(|x| x.dyn_into::<web_sys::Element>().ok())
        else {
            return;
        };
        if let Some(i) = target
            .get_attribute("data-point")
            .and_then(|v| v.parse::<usize>().ok())
        {
            if let Some(p) = view.borrow_mut().points.get_mut(i) {
                p.included = target.dyn_into::<HtmlInputElement>().unwrap().checked();
            }
            points(&view.borrow())
        }
    });
    el("observation-points")
        .add_event_listener_with_callback("change", toggle.as_ref().unchecked_ref())?;
    toggle.forget();
    let view = state.clone();
    let delete = Closure::<dyn FnMut(web_sys::Event)>::new(move |e: web_sys::Event| {
        let Some(target) = e
            .target()
            .and_then(|x| x.dyn_into::<web_sys::Element>().ok())
        else {
            return;
        };
        if let Some(i) = target
            .get_attribute("data-delete")
            .and_then(|v| v.parse::<usize>().ok())
        {
            if i < view.borrow().points.len() {
                view.borrow_mut().points.remove(i);
            }
            points(&view.borrow())
        }
    });
    el("observation-points")
        .add_event_listener_with_callback("click", delete.as_ref().unchecked_ref())?;
    delete.forget();
    load(state, initial).await
}
