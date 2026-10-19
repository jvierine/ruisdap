//! Bounded, read-only access to native Millstone products. Source paths come
//! exclusively from the server index; URL arguments never become file paths.
use axum::{
    extract::{Path, Query, State},
    http::{header, StatusCode},
    response::{IntoResponse, Response},
    routing::get,
    Json, Router,
};
use serde::Deserialize;
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, VecDeque},
    path::PathBuf,
    sync::{Arc, Mutex},
};
#[derive(Deserialize, Clone)]
struct PlasmaSource {
    path: PathBuf,
    key: String,
    antenna: String,
    date: String,
    start: f64,
    end: f64,
}
#[derive(Deserialize, Clone)]
struct RawSource {
    time: f64,
    path: PathBuf,
}
#[derive(Deserialize, Clone)]
struct Dataset {
    raw: Vec<RawSource>,
    manifest: PathBuf,
}
#[derive(Deserialize, Default)]
struct Index {
    plasma: BTreeMap<String, PlasmaSource>,
    datasets: BTreeMap<String, Dataset>,
}
struct Plasma {
    meta: Value,
    bytes: Vec<u8>,
}
#[derive(Clone)]
struct Store {
    index: Arc<Index>,
    cache: Arc<Mutex<VecDeque<(String, Arc<Plasma>)>>>,
    public: PathBuf,
    calibration_root: PathBuf,
    token: Arc<String>,
}
type Error = (StatusCode, String);
fn fail(e: impl std::fmt::Display) -> Error {
    (StatusCode::BAD_REQUEST, e.to_string())
}
#[repr(C)]
#[derive(Clone, Copy)]
struct C32 {
    r: f32,
    i: f32,
}
unsafe impl hdf5::H5Type for C32 {
    fn type_descriptor() -> hdf5::types::TypeDescriptor {
        use hdf5::types::{CompoundField, CompoundType, TypeDescriptor};
        TypeDescriptor::Compound(CompoundType {
            fields: vec![
                CompoundField::typed::<f32>("r", 0, 0),
                CompoundField::typed::<f32>("i", 4, 1),
            ],
            size: 8,
        })
    }
}
fn plasma(s: &Store, id: &str) -> anyhow::Result<Arc<Plasma>> {
    if let Some((_, v)) = s.cache.lock().unwrap().iter().find(|(k, _)| k == id) {
        return Ok(v.clone());
    }
    let source = s
        .index
        .plasma
        .get(id)
        .ok_or_else(|| anyhow::anyhow!("Unknown plasma record"))?;
    let h = hdf5::File::open(&source.path)?;
    let g = h.group(&source.key)?;
    let d = g.dataset("spec")?;
    let shape = d.shape();
    anyhow::ensure!(shape.len() == 2, "Invalid plasma spectrum");
    let nf = shape[0];
    let nr = shape[1];
    anyhow::ensure!(nf <= 65536 && nr <= 2048, "Spectrum exceeds viewer bounds");
    let mut frequency = g.dataset("freqs")?.read_raw::<f64>()?;
    let ranges = g.dataset("ranges")?.read_raw::<f64>()?;
    let sr = g.dataset("sr")?.read_scalar::<f64>()?;
    let span = frequency.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - frequency.iter().copied().fold(f64::INFINITY, f64::min);
    if ((span * 1e6 - sr) / sr).abs() < 0.002 {
        for f in &mut frequency {
            *f *= 1e6;
        }
    } else {
        anyhow::ensure!(
            ((span - sr) / sr).abs() < 0.002,
            "Unverified frequency units"
        );
    }
    let raw = d.read_raw::<f32>()?;
    anyhow::ensure!(
        frequency.len() == nf && ranges.len() == nr,
        "Spectrum coordinates disagree"
    );
    let mut order: Vec<_> = (0..nf).collect();
    order.sort_by(|a, b| frequency[*a].total_cmp(&frequency[*b]));
    let avg = 5;
    let bins = nf / avg;
    let mut values = vec![f32::NAN; nr * bins];
    let mut freq = Vec::with_capacity(bins);
    for (fi, group) in order.chunks_exact(avg).enumerate() {
        freq.push(group.iter().map(|i| frequency[*i]).sum::<f64>() / avg as f64);
        for r in 0..nr {
            let mut sum = 0.;
            let mut count = 0;
            for i in group {
                let v = raw[*i * nr + r];
                if v.is_finite() {
                    sum += v;
                    count += 1;
                }
            }
            if count > 0 {
                values[r * bins + fi] = sum / count as f32;
            }
        }
    }
    let bytes = values.iter().flat_map(|v| v.to_le_bytes()).collect();
    let result = Arc::new(Plasma {
        meta: json!({"id":id,"antenna":source.antenna,"date":source.date,"start":source.start,"end":source.end,"frequency_hz":freq,"range_m":ranges.iter().map(|r|r*1000.).collect::<Vec<_>>(),"frequency_bins_averaged":avg,"units":"source spectral power","values_url":format!("api/catalog/millstone/plasma/{id}/values.f32")}),
        bytes,
    });
    let mut cache = s.cache.lock().unwrap();
    cache.push_back((id.into(), result.clone()));
    while cache.len() > 6 {
        cache.pop_front();
    }
    Ok(result)
}
async fn metadata(State(s): State<Store>, Path(id): Path<String>) -> Result<Json<Value>, Error> {
    let p = tokio::task::spawn_blocking(move || plasma(&s, &id))
        .await
        .map_err(fail)?
        .map_err(fail)?;
    Ok(Json(p.meta.clone()))
}
async fn values(State(s): State<Store>, Path(id): Path<String>) -> Result<Response, Error> {
    let p = tokio::task::spawn_blocking(move || plasma(&s, &id))
        .await
        .map_err(fail)?
        .map_err(fail)?;
    Ok((
        [
            (header::CONTENT_TYPE, "application/octet-stream"),
            (header::CACHE_CONTROL, "private, max-age=3600"),
        ],
        p.bytes.clone(),
    )
        .into_response())
}
#[derive(Deserialize)]
struct Selection {
    date: String,
    antenna: String,
    time: f64,
    range_m: f64,
}
async fn acf(State(s): State<Store>, Query(q): Query<Selection>) -> Result<Json<Value>, Error> {
    if !q.time.is_finite() || !q.range_m.is_finite() {
        return Err(fail("Invalid selection"));
    }
    let v=tokio::task::spawn_blocking(move||->anyhow::Result<Value>{
  let d=s.index.datasets.get(&format!("{}/{}",q.date,q.antenna)).ok_or_else(||anyhow::anyhow!("Unknown Millstone dataset"))?;
  let source=d.raw.iter().min_by(|a,b|(a.time-q.time).abs().total_cmp(&(b.time-q.time).abs())).ok_or_else(||anyhow::anyhow!("No native ACFs"))?;
  anyhow::ensure!((source.time-q.time).abs()<=300.,"No native ACF within 300 seconds");
  let h=hdf5::File::open(&source.path)?;let actual_time=h.dataset("i0")?.read_scalar::<f64>()?;anyhow::ensure!((actual_time-source.time).abs()<0.01,"Native ACF time differs from index");let ranges=h.dataset("rgs_km")?.read_raw::<f64>()?;anyhow::ensure!(!ranges.is_empty() && ranges.iter().all(|r|r.is_finite()),"Invalid native range coordinates");
  let ri=ranges.iter().enumerate().min_by(|(_,a),(_,b)|(**a*1000.-q.range_m).abs().total_cmp(&(**b*1000.-q.range_m).abs())).map(|(i,_)|i).unwrap();
  let lag=h.dataset("lags")?.read_raw::<f64>()?;let n=lag.len();let data=h.dataset("acfs_e")?.read_raw::<C32>()?;let variance=h.dataset("acfs_var")?.read_raw::<f32>()?;let alpha=h.dataset("alpha")?.read_scalar::<f64>()?;
  anyhow::ensure!(alpha.is_finite()&&alpha>0.&&n>0&&data.len()==ranges.len()*n&&variance.len()==data.len(),"Invalid native ACF dimensions or receiver gain");
  let observed:Vec<_>=data[ri*n..(ri+1)*n].iter().map(|v|[v.r as f64/alpha,v.i as f64/alpha]).collect();
  let sigma:Vec<_>=variance[ri*n..(ri+1)*n].iter().map(|v|f64::from(*v).max(0.).sqrt()/alpha).collect();
  Ok(json!({"time":source.time,"range_m":ranges[ri]*1000.,"lag_s":lag,"observed":observed,"sigma":sigma,"units":"K","receiver_gain":alpha,"source":"native un-subtracted ACF; receiver gain removed","fit_averaged":false}))
 }).await.map_err(fail)?.map_err(fail)?;
    Ok(Json(v))
}
pub fn router(state: PathBuf, public: PathBuf, calibration_root: PathBuf) -> Router {
    let token = Arc::new(
        std::fs::read_to_string(state.join("processing-key"))
            .unwrap_or_default()
            .trim()
            .to_string(),
    );
    let index = std::fs::read(state.join("millstone-sources.json"))
        .ok()
        .and_then(|b| serde_json::from_slice(&b).ok())
        .unwrap_or_default();
    Router::new()
        .route("/guisdap/api/catalog/millstone/plasma/{id}", get(metadata))
        .route(
            "/guisdap/api/catalog/millstone/plasma/{id}/values.f32",
            get(values),
        )
        .route("/guisdap/api/catalog/millstone/acf", get(acf))
        .route(
            "/guisdap/api/catalog/millstone/calibrate",
            axum::routing::post(recalibrate),
        )
        .with_state(Store {
            index: Arc::new(index),
            cache: Arc::new(Mutex::new(VecDeque::new())),
            public,
            calibration_root,
            token,
        })
}

include!("millstone_calibration.rs");

include!("millstone_calibration_tests.rs");
