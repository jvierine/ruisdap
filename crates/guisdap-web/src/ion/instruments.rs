//! Instrument diagnostics retain native time coverage and share the fit time window.
use super::*;
#[derive(Deserialize)]
pub(super) struct Data {
    records: Vec<Record>,
}
#[derive(Deserialize)]
struct Record {
    start_unix: f64,
    end_unix: f64,
    system_temperature_k: Option<f64>,
    transmitter_w: f64,
    azimuth_deg: f64,
    elevation_deg: f64,
}
pub(super) fn render(v: &Viewer) {
    let Some(data)=&v.instruments else{text("fit-instrument-summary","Instrument diagnostics unavailable.");return;};
    let records:Vec<_>=data.records.iter().map(|r|serde_json::json!({"start":r.start_unix,"end":r.end_unix,"azimuth_deg":r.azimuth_deg,"elevation_deg":r.elevation_deg,"system_temperature_k":r.system_temperature_k,"transmitter_w":r.transmitter_w})).collect();
    let p=&v.data.profiles[v.time];
    crate::plotting::shared_curve("fit-instrument-0",&serde_json::json!({"instrument_records":records,"bounds":[v.window_start,v.window_end,0.,1.],"selected":(p.start_unix+p.end_unix)/2.}).to_string());
    text("fit-instrument-summary","");
}
