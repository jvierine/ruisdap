//! Static-hosted browser application. Data selection, scaling, UI and WebGL live in Rust.
#![cfg(target_arch = "wasm32")]
mod calibration;
mod catalog;
mod ion;
mod live;
mod navigation;
mod observation;
mod plotting;
use guisdap_core::{merge_bands, native_spectra};
use serde::Deserialize;
use std::{cell::RefCell, rc::Rc};
use wasm_bindgen::{prelude::*, JsCast};
use wasm_bindgen_futures::{spawn_local, JsFuture};
use web_sys::{HtmlInputElement, HtmlSelectElement};
#[derive(Deserialize)]
struct Config {
    data_root: String,
    #[serde(default)]
    datasets: Vec<PlasmaDataset>,
}
#[derive(Deserialize)]
struct PlasmaDataset {
    label: String,
    data_root: String,
}
#[derive(Deserialize)]
struct Frame {
    #[serde(default)]
    source_count: usize,
    id: String,
    start_unix: f64,
    end_unix: f64,
    elevation_deg: f64,
    azimuth_deg: f64,
    transmitter_w: f64,
    centers_hz: [f64; 2],
}
#[derive(Deserialize)]
struct Manifest {
    schema: String,
    frames: Vec<Frame>,
    ion_range_km: Vec<f64>,
    decoded_range_km: Vec<f64>,
    #[serde(default)]
    overview_root: String,
}
struct App {
    calibration: Option<calibration::Cal>,
    profile_altitude: [f64; 2],
    profile_power: [f64; 2],
    broad_frequency: [f64; 2],
    root: String,
    m: Manifest,
    broad: Vec<f32>,
    raw: Vec<f32>,
    power: Vec<f32>,
    power_sigma: Vec<f32>,
    decoded: Vec<f32>,
    selected: usize,
    view_start: f64,
    view_end: f64,
    decoded_selected: Option<usize>,
    raw_view: bool,
    percentile: f64,
    busy: bool,
    request: u64,
    smoothing: usize,
    magnitude: bool,
    weighted: bool,
    limit: Option<f32>,
    native: Vec<f32>,
    cache: std::collections::VecDeque<(usize, Vec<f32>)>,
}
type State = Rc<RefCell<App>>;
fn experiment_config(default: &str) -> String {
    let query = web_sys::window()
        .unwrap()
        .location()
        .search()
        .unwrap_or_default();
    let id = query
        .trim_start_matches('?')
        .split('&')
        .find_map(|p| p.strip_prefix("experiment="));
    match id {
        Some(id) if id.len() == 24 && id.bytes().all(|b| b.is_ascii_hexdigit()) => {
            format!("catalog/{id}/{default}")
        }
        _ => default.into(),
    }
}
fn doc() -> web_sys::Document {
    web_sys::window().unwrap().document().unwrap()
}
fn el(id: &str) -> web_sys::Element {
    doc().get_element_by_id(id).unwrap()
}
fn text(id: &str, s: &str) {
    el(id).set_text_content(Some(s));
}
fn error(s: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&s.to_string())
}
fn report(e: JsValue) {
    text(
        "status",
        &format!(
            "Load failed: {}",
            e.as_string().unwrap_or_else(|| format!("{e:?}"))
        ),
    );
    web_sys::console::error_1(&e);
}
async fn response(url: &str) -> Result<web_sys::Response, JsValue> {
    let request = web_sys::RequestInit::new();
    if url.ends_with("config.json") || url.ends_with("index.json") || url.ends_with("status.json") {
        request.set_cache(web_sys::RequestCache::NoCache);
    }
    let r: web_sys::Response = JsFuture::from(
        web_sys::window()
            .unwrap()
            .fetch_with_str_and_init(url, &request),
    )
    .await?
    .dyn_into()?;
    if !r.ok() {
        return Err(error(format!("HTTP {}: {url}", r.status())));
    }
    Ok(r)
}
async fn json<T: serde::de::DeserializeOwned>(url: &str) -> Result<T, JsValue> {
    let r = response(url).await?;
    let s = JsFuture::from(r.text()?)
        .await?
        .as_string()
        .ok_or_else(|| error("expected text"))?;
    serde_json::from_str(&s).map_err(error)
}
async fn floats(url: &str) -> Result<Vec<f32>, JsValue> {
    let r = response(url).await?;
    let b = JsFuture::from(r.array_buffer()?).await?;
    let bytes = js_sys::Uint8Array::new(&b).to_vec();
    if !bytes.len().is_multiple_of(4) {
        return Err(error("truncated float product"));
    }
    let v: Vec<f32> = bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|b| f32::from_le_bytes(*b))
        .collect();
    if !v.iter().all(|v| v.is_finite()) {
        return Err(error("non-finite exported product"));
    }
    Ok(v)
}
fn utc(t: f64) -> String {
    let s: String = js_sys::Date::new(&JsValue::from_f64(t * 1000.))
        .to_iso_string()
        .into();
    s[11..19].into()
}
fn altitude(r: f64, e: f64) -> f64 {
    let earth = 6370.03;
    (earth * earth + r * r + 2. * earth * r * e.to_radians().sin()).sqrt() - 6370.
}
fn percentile(v: &[f32], p: f64) -> f32 {
    let mut a: Vec<f32> = v
        .iter()
        .map(|x| x.abs())
        .filter(|x| x.is_finite())
        .collect();
    if a.is_empty() {
        return 1.;
    }
    let index = ((a.len() - 1) as f64 * p) as usize;
    (*a.select_nth_unstable_by(index, f32::total_cmp).1).max(1e-12)
}
fn draw(
    id: &str,
    v: &[f32],
    width: usize,
    height: usize,
    p: f64,
    cursor: Option<f32>,
    scale: (Option<f32>, bool),
    bounds: [f64; 4],
) -> Result<(), JsValue> {
    if v.len() != width * height || width == 0 || height == 0 {
        return Err(error("invalid field dimensions"));
    }
    let limit = scale.0.unwrap_or_else(|| percentile(v, p));
    let positive = scale.1;
    let limits = if id == "power" {
        [1e8, 1e13]
    } else if positive {
        [0., limit as f64]
    } else {
        [-(limit as f64), limit as f64]
    };
    plotting::shared_field(
        id,
        v,
        width,
        height,
        &serde_json::json!({"bounds":bounds,"limits":limits,"mode":if id=="power"{1}else{2},"normWidth":limit as f64/8.,"palette":1,"colormap":"plasma","title":if id=="power"{"Raw Ne / m⁻³"}else{"Plasma-line power / source units"},"time":id!="decoded","external_navigation":true,"xlabel":"Frequency offset / MHz","ylabel":if id.starts_with("broad"){"Frequency offset / MHz"}else{"Altitude / km"}}),
    )?;
    let canvas = el(id);
    let cursor_id = format!("{id}-time-cursor");
    if doc().get_element_by_id(&cursor_id).is_none() {
        let marker = doc().create_element("div")?;
        marker.set_id(&cursor_id);
        canvas.parent_element().unwrap().append_child(&marker)?;
    }
    el(&cursor_id).set_attribute("style",&format!("position:absolute;z-index:3;top:22px;bottom:58px;left:calc(76px + (100% - 208px)*{});border-left:1px solid #986327;pointer-events:none;display:{}",cursor.unwrap_or(0.),if cursor.is_some_and(|v|(0. ..=1.).contains(&v)){"block"}else{"none"}))?;
    Ok(())
}

fn overview(a: &App) -> Result<(), JsValue> {
    let n = a.m.frames.len();
    let frame = &a.m.frames[a.selected];
    let cursor = Some(
        (((frame.start_unix + frame.end_unix) / 2. - a.view_start) / (a.view_end - a.view_start))
            as f32,
    );
    let dt = a.m.frames[0].end_unix - a.m.frames[0].start_unix;
    let width = ((a.view_end - a.view_start) / dt).ceil().max(1.) as usize;
    let columns: Vec<Option<usize>> = (0..width)
        .map(|i| {
            let t = a.view_start + (a.view_end - a.view_start) * (i as f64 + 0.5) / width as f64;
            let j = a.m.frames.partition_point(|f| f.end_unix <= t);
            a.m.frames.get(j).filter(|f| f.start_unix <= t).map(|_| j)
        })
        .collect();
    let axes = doc().get_elements_by_class_name("time-axis");
    for i in 0..axes.length() {
        axes.item(i).unwrap().set_inner_html(&format!(
            "<span>{}</span><span>UTC · gaps left empty</span><span>{}</span>",
            utc(a.view_start),
            utc(a.view_end)
        ));
    }
    text(
        "plasma-window",
        &format!("{}–{} UTC", utc(a.view_start), utc(a.view_end)),
    );
    let time: HtmlInputElement = el("time").dyn_into()?;
    time.set_min(
        &a.m.frames
            .partition_point(|f| f.end_unix <= a.view_start)
            .min(n - 1)
            .to_string(),
    );
    time.set_max(
        &a.m.frames
            .partition_point(|f| f.start_unix < a.view_end)
            .saturating_sub(1)
            .to_string(),
    );
    let data = if a.raw_view { &a.raw } else { &a.broad };
    // Convert band/frequency/time to band/time/frequency before frequency merging.
    let mut all_history = Vec::new();
    for gate in 0..4 {
        let mut bands = vec![0.; 2 * n * 100];
        for band in 0..2 {
            for t in 0..n {
                for f in 0..100 {
                    bands[(band * n + t) * 100 + f] = data[((band * 4 + gate) * 100 + f) * n + t];
                }
            }
        }
        let (merged, bins, low) =
            merge_bands(&bands, n, 100, a.m.frames[0].centers_hz, 25000.).map_err(error)?;
        let mut history = vec![0.; bins * n];
        for t in 0..n {
            for f in 0..bins {
                history[f * n + t] = merged[t * bins + f];
            }
        }
        all_history.push((history, bins, low));
    }
    let combined: Vec<f32> = all_history
        .iter()
        .flat_map(|(h, _, _)| h.iter().copied())
        .collect();
    let shared_limit = percentile(&combined, a.percentile);
    let bounds = [[47., 221.], [128., 302.], [209., 383.], [290., 464.]];
    let elevation = a.m.frames[a.selected].elevation_deg;
    for (gate, (history, bins, low)) in all_history.iter().enumerate() {
        let first =
            (((a.broad_frequency[0] - low) / 25000.).ceil().max(0.) as usize).min(*bins - 1);
        let last = (((a.broad_frequency[1] - low) / 25000.)
            .floor()
            .max(first as f64) as usize)
            .min(*bins - 1);
        text(
            &format!("broad-high-{gate}"),
            &format!("{:.3}", (low + last as f64 * 25000.) / 1e6),
        );
        text(
            &format!("broad-low-{gate}"),
            &format!("{:.3}", (low + first as f64 * 25000.) / 1e6),
        );
        let visible: Vec<f32> = (first..=last)
            .flat_map(|f| {
                columns
                    .iter()
                    .map(move |t| t.map_or(f32::NAN, |t| history[f * n + t]))
            })
            .collect();
        draw(
            &format!("broad{gate}"),
            &visible,
            width,
            last - first + 1,
            a.percentile,
            cursor,
            (Some(shared_limit), a.raw_view),
            [
                a.view_start,
                a.view_end,
                (low + first as f64 * 25000.) / 1e6,
                (low + last as f64 * 25000.) / 1e6,
            ],
        )?;
        text(
            &format!("gate-altitude{gate}"),
            &format!(
                "{:.0}–{:.0} km altitude · broad gate {}",
                altitude(bounds[gate][0], elevation),
                altitude(bounds[gate][1], elevation),
                gate + 1
            ),
        );
    }
    let elevation = a.m.frames[0].elevation_deg;
    let shift = 0.;
    let minimum = a.profile_altitude[0];
    let indices: Vec<usize> =
        a.m.ion_range_km
            .iter()
            .enumerate()
            .filter_map(|(i, &r)| {
                let h = altitude(r - shift, elevation);
                if (minimum..=a.profile_altitude[1]).contains(&h) {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();
    let mut power = vec![];
    for &r in &indices {
        power.extend(
            columns
                .iter()
                .map(|t| t.map_or(f32::NAN, |t| a.power[r * n + t])),
        );
    }
    draw(
        "power",
        &power,
        width,
        indices.len(),
        a.percentile,
        cursor,
        (Some(1e13), true),
        [
            a.view_start,
            a.view_end,
            altitude(a.m.ion_range_km[indices[0]], elevation),
            altitude(a.m.ion_range_km[*indices.last().unwrap()], elevation),
        ],
    )?;
    let low = altitude(a.m.ion_range_km[indices[0]] - shift, elevation);
    let high = altitude(
        a.m.ion_range_km[*indices.last().unwrap()] - shift,
        elevation,
    );
    text("power-axis-label", "Response-centroid altitude / km");
    text("power-axis-caption", "10 µs lag raw Ne / m⁻³ versus UTC");
    let r0 = a.m.ion_range_km[0];
    let half = 2.1724186 / 2.;
    text(
        "power-support",
        &format!(
            "First recorded gate response: {:.1}–{:.1} km altitude.",
            altitude(r0 - half, elevation),
            altitude(r0 + half, elevation)
        ),
    );
    let mut ticks = String::new();
    for h in (0..=1000).step_by(10) {
        let h = h as f64;
        if h < low || h > high {
            continue;
        }
        let top = 100. * (high - h) / (high - low);
        if (h as usize).is_multiple_of(50) {
            ticks += &format!("<span class='major' style='top:{top:.3}%'>{h:.0}</span>");
        } else {
            ticks += &format!("<span class='minor' style='top:{top:.3}%'></span>");
        }
    }
    el("power-ticks").set_inner_html(&ticks);
    overview_selection(a)
}
fn overview_selection(a: &App) -> Result<(), JsValue> {
    let frame = &a.m.frames[a.selected];
    let cursor =
        ((frame.start_unix + frame.end_unix) / 2. - a.view_start) / (a.view_end - a.view_start);
    for id in ["broad0", "broad1", "broad2", "broad3", "power"] {
        if let Some(marker) = doc().get_element_by_id(&format!("{id}-time-cursor")) {
            marker.set_attribute("style", &format!("position:absolute;top:0;bottom:0;left:{}%;border-left:1px solid white;pointer-events:none;display:{}", cursor*100., if (0. ..=1.).contains(&cursor) {"block"} else {"none"}))?;
        }
    }
    let bounds = [[47., 221.], [128., 302.], [209., 383.], [290., 464.]];
    for (gate, range) in bounds.iter().enumerate() {
        text(
            &format!("gate-altitude{gate}"),
            &format!(
                "{:.0}–{:.0} km altitude · broad gate {}",
                altitude(range[0], frame.elevation_deg),
                altitude(range[1], frame.elevation_deg),
                gate + 1
            ),
        );
    }
    let indices: Vec<_> =
        a.m.ion_range_km
            .iter()
            .enumerate()
            .filter_map(|(i, &r)| {
                (a.profile_altitude[0]..=a.profile_altitude[1])
                    .contains(&altitude(r, frame.elevation_deg))
                    .then_some(i)
            })
            .collect();
    profile(a, &indices);
    Ok(())
}
fn profile(a: &App, indices: &[usize]) {
    let n = a.m.frames.len();
    let f = &a.m.frames[a.selected];
    let [bottom, top] = a.profile_altitude;
    let [xmin, xmax] = a.profile_power;
    let mut points = Vec::new();
    let mut errors = Vec::new();
    let mut lower = Vec::new();
    let mut upper = Vec::new();
    let mut nonpositive = 0;
    for &i in indices {
        let h = altitude(a.m.ion_range_km[i], f.elevation_deg);
        if h < bottom || h > top {
            continue;
        }
        let v = a.power[i * n + a.selected] as f64;
        let se = a.power_sigma[i * n + a.selected] as f64;
        if v <= 0. {
            nonpositive += 1;
        }
        if v - se < xmax && v + se > xmin {
            errors.extend([(v - se).max(xmin).log10(), h, (v + se).min(xmax).log10(), h]);
        }
        if v < xmin {
            lower.extend([xmin.log10(), h]);
            points.extend([f64::NAN, f64::NAN]);
        } else if v > xmax {
            upper.extend([xmax.log10(), h]);
            points.extend([f64::NAN, f64::NAN]);
        } else {
            points.extend([v.log10(), h]);
        }
    }
    plotting::shared_curve("profile",&serde_json::json!({"title":"Raw Ne profile","height":380,"bounds":[xmin.log10(),xmax.log10(),bottom,top],"xlabel":"Raw Ne / m⁻³","xformat":"log10","ylabel":"Altitude / km","series":[{"xy":errors,"segments":true,"color":"#7b929e"},{"xy":points,"radius":1.5},{"xy":lower,"style":2,"radius":3},{"xy":upper,"style":3,"radius":3}]}).to_string());
    if let Some(c) = &a.calibration {
        if let Ok(c) = c.try_borrow() {
            calibration::markers(a, &c);
        }
    }
    text("power-uncertainty",&format!("Horizontal bars: approximate ±1σ from 14 short integrations plus calibration. Clutter-induced covariance and systematic errors are not included. {nonpositive} nonpositive estimates in the displayed altitude interval; limit markers show values outside the log axis. Signed estimates remain in HDF5."));
}

fn decoded(a: &mut App) -> Result<(), JsValue> {
    if a.decoded_selected != Some(a.selected) {
        return Ok(());
    }
    let bins = guisdap_core::NFFT;
    let mut data = a.decoded.clone();
    if a.magnitude {
        for x in &mut data {
            *x = x.abs();
        }
    }
    if a.smoothing > 1 {
        let source = data.clone();
        for (dest, row) in data.chunks_mut(bins).zip(source.chunks(bins)) {
            for (i, v) in dest.iter_mut().enumerate() {
                let left = i.saturating_sub(a.smoothing / 2);
                let right = (i + (a.smoothing - 1) / 2 + 1).min(bins);
                *v = row[left..right].iter().sum::<f32>() / (right - left) as f32;
            }
        }
    }
    if a.magnitude {
        // Detection contrast, as in ipypl: remove frequency-dependent noise amplitude.
        // This is explicitly not a calibrated spectrum or a power estimate.
        for b in 0..2 {
            for f in 0..bins {
                let mut column: Vec<f32> = (0..90).map(|r| data[(b * 90 + r) * bins + f]).collect();
                column.sort_by(f32::total_cmp);
                let median = ((column[44] + column[45]) * 0.5).max(1e-20);
                for r in 0..90 {
                    data[(b * 90 + r) * bins + f] /= median;
                }
            }
        }
    }
    let frame = &a.m.frames[a.selected];
    let (merged, width, low) =
        merge_bands(&data, 90, bins, frame.centers_hz, 2.5e6 / bins as f64).map_err(error)?;
    let limit = *a.limit.get_or_insert_with(|| {
        if a.magnitude {
            6.
        } else {
            percentile(&merged, a.percentile)
        }
    });
    draw(
        "decoded",
        &merged,
        width,
        90,
        a.percentile,
        None,
        (Some(limit), a.magnitude),
        [
            low / 1e6,
            (low + (width - 1) as f64 * 2.5e6 / bins as f64) / 1e6,
            altitude(a.m.decoded_range_km[0], frame.elevation_deg),
            altitude(a.m.decoded_range_km[89], frame.elevation_deg),
        ],
    )?;
    text(
        "decoded-range",
        &format!("{} — {} UTC", utc(frame.start_unix), utc(frame.end_unix)),
    );
    text("freq-min", &format!("{:.3}", low / 1e6));
    text(
        "freq-max",
        &format!(
            "{:.3}",
            (low + (width - 1) as f64 * 2.5e6 / bins as f64) / 1e6
        ),
    );
    for (class, i) in [("alt-min", 0), ("alt-max", 89)] {
        let collection = doc().get_elements_by_class_name(class);
        for k in 0..collection.length() {
            collection.item(k).unwrap().set_text_content(Some(&format!(
                "{:.0}",
                altitude(a.m.decoded_range_km[i], frame.elevation_deg)
            )));
        }
    }
    Ok(())
}
async fn select(state: State, index: usize) -> Result<(), JsValue> {
    let (url, request, cached) = {
        let mut a = state.borrow_mut();
        a.request += 1;
        a.busy = true;
        text("loading", "Loading next integration…");
        (
            format!("{}frames/{}/native-long.f32", a.root, a.m.frames[index].id),
            a.request,
            a.cache
                .iter()
                .find(|(i, _)| *i == index)
                .map(|(_, v)| v.clone()),
        )
    };
    let result = match cached {
        Some(v) => Ok(v),
        None => floats(&url).await,
    };
    let mut a = state.borrow_mut();
    if a.request != request {
        return Ok(());
    }
    a.busy = false;
    text("loading", "");
    let data = match result {
        Ok(d) => d,
        Err(e) => {
            return Err(e);
        }
    };
    if data.len() != 2 * 800 * 90 * 2 {
        return Err(error("invalid native lag dimensions"));
    }
    let render_start = js_sys::Date::now();
    a.decoded = native_spectra(&data, a.weighted).map_err(error)?;
    if !a.cache.iter().any(|(i, _)| *i == index) {
        a.cache.push_back((index, data.clone()));
    }
    while a.cache.len() > 6 {
        a.cache.pop_front();
    }
    let first_draw = a.native.is_empty();
    a.native = data;
    a.selected = index;
    a.decoded_selected = Some(index);
    let f = &a.m.frames[index];
    text(
        "selected",
        &format!("{} — {} UTC", utc(f.start_unix), utc(f.end_unix)),
    );
    text(
        "frameinfo",
        &format!(
            "Integration {} / {} · {} / {:.0} native records · {:.2} MW",
            index + 1,
            a.m.frames.len(),
            f.source_count.max(1),
            ((f.end_unix - f.start_unix) / 5.).max(1.),
            f.transmitter_w / 1e6
        ),
    );
    text(
        "pointing",
        &format!("{:.1}° / {:.1}°", f.azimuth_deg, f.elevation_deg),
    );
    text(
        "profile-time",
        &format!("{} — {} UTC", utc(f.start_unix), utc(f.end_unix)),
    );
    el("time")
        .dyn_into::<HtmlInputElement>()?
        .set_value(&index.to_string());
    // Commit all panels together after fetch and FFT; never display mixed timestamps.
    if first_draw {
        overview(&a)?;
    } else {
        overview_selection(&a)?;
    }
    decoded(&mut a)?;
    el("decoded").set_attribute(
        "data-update-ms",
        &(js_sys::Date::now() - render_start).to_string(),
    )?;
    Ok(())
}
fn on(
    id: &str,
    event: &str,
    callback: impl FnMut(web_sys::Event) + 'static,
) -> Result<(), JsValue> {
    let c = Closure::<dyn FnMut(web_sys::Event)>::new(callback);
    el(id).add_event_listener_with_callback(event, c.as_ref().unchecked_ref())?;
    c.forget();
    Ok(())
}
#[wasm_bindgen]
pub async fn start() -> Result<(), JsValue> {
    navigation::setup()?;
    let cfg: Config = json(&experiment_config("config.json")).await?;
    let query = web_sys::window().unwrap().location().search()?;
    let param = |key: &str| {
        query.trim_start_matches('?').split('&').find_map(|part| {
            part.strip_prefix(&format!("{key}="))
                .and_then(|v| v.parse::<f64>().ok())
        })
    };
    let chosen = param("resolution").unwrap_or(0.) as usize;
    if let Some(selector) = doc().get_element_by_id("plasma-resolution") {
        selector.set_inner_html(
            &cfg.datasets
                .iter()
                .enumerate()
                .map(|(i, d)| format!("<option value='{i}'>{}</option>", d.label))
                .collect::<String>(),
        );
        selector
            .clone()
            .dyn_into::<HtmlSelectElement>()?
            .set_value(&chosen.to_string());
    }
    let data_root = cfg
        .datasets
        .get(chosen)
        .map_or(cfg.data_root.as_str(), |d| d.data_root.as_str())
        .to_string();
    let root = if data_root.ends_with('/') {
        data_root
    } else {
        data_root + "/"
    };
    let m: Manifest = json(&(root.clone() + "index.json")).await?;
    if m.schema != "beata22-quicklook-v5"
        || m.frames.is_empty()
        || m.ion_range_km.len() != 431
        || m.decoded_range_km.len() != 90
    {
        return Err(error("unsupported or empty dataset"));
    }
    let overview_root = root.clone() + &m.overview_root;
    let urls = [
        "broad-excess.f32",
        "broad-raw.f32",
        "ion-short-raw-ne.f32",
        "ion-short-raw-ne-sigma.f32",
    ]
    .map(|name| overview_root.clone() + name);
    let (broad, raw, power, power_sigma) = futures_util::try_join!(
        floats(&urls[0]),
        floats(&urls[1]),
        floats(&urls[2]),
        floats(&urls[3])
    )?;
    if let Some(button) = doc().get_element_by_id("refresh-data") {
        let cb = Closure::<dyn FnMut(web_sys::Event)>::new(move |_| {
            let _ = web_sys::window().unwrap().location().reload();
        });
        button.add_event_listener_with_callback("click", cb.as_ref().unchecked_ref())?;
        cb.forget();
    }
    let n = m.frames.len();
    if broad.len() != 800 * n
        || raw.len() != 800 * n
        || power.len() != 431 * n
        || power_sigma.len() != 431 * n
    {
        return Err(error("invalid overview dimensions"));
    }
    let first = &m.frames[0];
    let last = m.frames.last().unwrap();
    text(
        "window",
        &format!("{} – {}", utc(first.start_unix), utc(last.end_unix)),
    );
    let dt = first.end_unix - first.start_unix;
    text("cadence", &format!("{dt:.1} s · {n} records"));
    let gaps = m
        .frames
        .windows(2)
        .filter(|w| (w[1].start_unix - w[0].end_unix).abs() > 0.1)
        .count();
    text(
        "status",
        &format!("{n} observed integrations · {gaps} time gaps"),
    );
    let axes = doc().get_elements_by_class_name("time-axis");
    for k in 0..axes.length() {
        axes.item(k).unwrap().set_inner_html(&format!(
            "<span>{}</span><span>UTC · one column per integration</span><span>{}</span>",
            utc(first.start_unix),
            utc(last.end_unix)
        ));
    }
    el("time")
        .dyn_into::<HtmlInputElement>()?
        .set_max(&(n - 1).to_string());
    let target = param("time").unwrap_or(1667408700.1);
    let initial = m
        .frames
        .partition_point(|f| f.end_unix <= target)
        .min(n - 1);
    let view_start = param("from")
        .unwrap_or(first.start_unix)
        .clamp(first.start_unix, last.end_unix - dt);
    let view_end = param("to")
        .unwrap_or(last.end_unix)
        .clamp(view_start + dt, last.end_unix);
    let state = Rc::new(RefCell::new(App {
        calibration: None,
        profile_altitude: [40., 400.],
        profile_power: [1e8, 1e13],
        broad_frequency: [-6.15e6, -1.275e6],
        root,
        m,
        broad,
        raw,
        power,
        power_sigma,
        decoded: vec![],
        selected: initial,
        view_start,
        view_end,
        decoded_selected: None,
        raw_view: false,
        percentile: 0.995,
        busy: false,
        request: 0,
        smoothing: 4,
        magnitude: true,
        weighted: true,
        limit: None,
        native: vec![],
        cache: std::collections::VecDeque::new(),
    }));
    if doc().get_element_by_id("plasma-resolution").is_some() {
        let s = state.clone();
        on("plasma-resolution", "change", move |_| {
            let a = s.borrow();
            let resolution = el("plasma-resolution")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value();
            let location = web_sys::window().unwrap().location();
            let query = location.search().unwrap_or_default();
            let experiment = query
                .trim_start_matches('?')
                .split('&')
                .find(|p| p.starts_with("experiment="))
                .unwrap_or("");
            let _ = location.set_search(&format!(
                "?resolution={resolution}&from={}&to={}&time={}&{experiment}",
                a.view_start,
                a.view_end,
                a.m.frames[a.selected].start_unix + 0.1
            ));
        })?;
    }
    for (id, factor) in [
        ("plasma-zoom-in", 0.5),
        ("plasma-zoom-out", 2.),
        ("plasma-view-all", 0.),
    ] {
        let s = state.clone();
        on(id, "click", move |_| {
            let mut a = s.borrow_mut();
            let first = a.m.frames[0].start_unix;
            let last = a.m.frames.last().unwrap().end_unix;
            let dt = a.m.frames[0].end_unix - first;
            let center = (a.m.frames[a.selected].start_unix + a.m.frames[a.selected].end_unix) / 2.;
            let width = if factor == 0. {
                last - first
            } else {
                ((a.view_end - a.view_start) * factor).clamp(dt, last - first)
            };
            a.view_start = (center - width / 2.).clamp(first, last - width);
            a.view_end = a.view_start + width;
            if let Err(e) = overview(&a) {
                report(e);
            }
        })?;
    }
    let s = state.clone();
    on("time", "input", move |_| {
        let s = s.clone();
        let i = el("time")
            .dyn_into::<HtmlInputElement>()
            .unwrap()
            .value()
            .parse::<usize>()
            .unwrap_or(0)
            .min(s.borrow().m.frames.len() - 1);
        spawn_local(async move {
            if let Err(e) = select(s, i).await {
                report(e);
            }
        });
    })?;
    for id in ["baseline", "stretch", "mode", "smooth", "weight"] {
        let s = state.clone();
        on(id, "change", move |_| {
            let mut a = s.borrow_mut();
            a.raw_view = el("baseline")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value()
                == "raw";
            a.percentile = match el("stretch")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value()
                .as_str()
            {
                "99" => 0.99,
                "999" => 0.999,
                "100" => 1.,
                _ => 0.995,
            };
            a.magnitude =
                el("mode").dyn_into::<HtmlSelectElement>().unwrap().value() == "magnitude";
            a.smoothing = el("smooth")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value()
                .parse()
                .unwrap_or(4);
            a.weighted = el("weight")
                .dyn_into::<HtmlSelectElement>()
                .unwrap()
                .value()
                == "sqrt";
            if matches!(id, "mode" | "smooth" | "weight" | "stretch") {
                a.limit = None;
            }
            if id == "weight" && !a.native.is_empty() {
                match native_spectra(&a.native, a.weighted) {
                    Ok(d) => a.decoded = d,
                    Err(e) => {
                        report(error(e));
                        return;
                    }
                }
            }
            if let Err(e) = overview(&a).and_then(|_| decoded(&mut a)) {
                report(e);
            }
        })?;
    }
    select(state.clone(), initial).await?;
    calibration::start(state.clone())?;
    Ok(())
}

/// Native measured-lag DFT shared by the EISCAT and Millstone inspectors.
#[wasm_bindgen]
pub fn acf_diagnostic(input: &str) -> Result<String, JsValue> {
    #[derive(serde::Deserialize)]
    struct Input {
        lag_s: Vec<f64>,
        observed: Vec<[Option<f64>; 2]>,
        sigma: Vec<Option<f64>>,
    }
    let v: Input = serde_json::from_str(input).map_err(error)?;
    let rows: Vec<_> = v
        .lag_s
        .iter()
        .zip(&v.observed)
        .zip(&v.sigma)
        .filter_map(|((t, a), s)| {
            let (Some(re), Some(im), Some(sigma)) = (a[0], a[1], *s) else {
                return None;
            };
            [*t, re, im, sigma]
                .iter()
                .all(|v| v.is_finite())
                .then_some((*t, guisdap_ion::Complex::new(re, im), [sigma, sigma]))
        })
        .collect();
    let lags: Vec<_> = rows.iter().map(|r| r.0).collect();
    let values: Vec<_> = rows.iter().map(|r| [r.1.re, r.1.im]).collect();
    let errors: Vec<_> = rows.iter().map(|r| r.2).collect();
    let freq: Vec<_> = (0..=800).map(|i| -40000. + i as f64 * 100.).collect();
    let op = guisdap_ion::dft::AcfDft::new(&lags, &freq, 1e-5).map_err(error)?;
    let power = op.transform(&values).map_err(error)?;
    let se = op.diagonal_standard_error(&errors).map_err(error)?;
    let acf = |imag: bool| {
        let f = |z: guisdap_ion::Complex| if imag { z.im } else { z.re };
        let high = rows
            .iter()
            .map(|r| f(r.1).abs() + r.2[0])
            .fold(1., f64::max)
            * 1.1;
        serde_json::json!({"title":if imag{"Imaginary ACF"}else{"Real ACF"},"bounds":[0.,lags.last().copied().unwrap_or(1e-5)*1e6,-high,high],"xlabel":"Measured lag / µs","ylabel":"ACF / K","series":[{"xy":rows.iter().flat_map(|r|[r.0*1e6,f(r.1)-r.2[0],r.0*1e6,f(r.1)+r.2[0]]).collect::<Vec<_>>(),"segments":true,"color":"#7b929e"},{"xy":rows.iter().flat_map(|r|[r.0*1e6,f(r.1)]).collect::<Vec<_>>(),"radius":2,"line":false}],"caption":"Native measured lags; source variance shown, lag covariance unavailable."})
    };
    let lo = power.iter().zip(&se).map(|(p, e)| p - e).fold(0., f64::min);
    let hi = power.iter().zip(&se).map(|(p, e)| p + e).fold(1., f64::max);
    let pad = (hi - lo) * 0.05;
    Ok(serde_json::json!({"real":acf(false),"imaginary":acf(true),"spectrum":{"title":"Measured ACF spectrum","bounds":[-40.,40.,lo-pad,hi+pad],"xlabel":"Frequency offset / kHz","ylabel":"ACF DFT / K Hz⁻¹","series":[{"xy":freq.iter().zip(&power).flat_map(|(f,p)|[*f/1000.,*p]).collect::<Vec<_>>()}],"bands":[{"values":(0..freq.len()).flat_map(|i|[freq[i]/1000.,power[i]-se[i],power[i]+se[i]]).collect::<Vec<_>>(),"color":"#236b8e4d"}],"caption":format!("{} measured lags · {} missing; no interpolation. Missing zero lag leaves the baseline undetermined. Diagonal source variance only.",op.unique_lags,op.missing_lags)}}).to_string())
}
