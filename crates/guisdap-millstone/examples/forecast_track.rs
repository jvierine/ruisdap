//! Diagnostic metadata-only line protocol for the native track forecaster.
use anyhow::Result;
use guisdap_millstone::track::{self, Event, Forecaster, Settings};
use serde::Deserialize;
use std::io::{BufRead, Write};
#[derive(Deserialize)]
struct Request {
    events: Vec<Event>,
    held_sample: u64,
    tx_centre: f64,
    #[serde(default)]
    settings: Settings,
}
fn main() -> Result<()> {
    let mut forecaster = Forecaster::default();
    let mut out = std::io::stdout().lock();
    for line in std::io::stdin().lock().lines() {
        let q: Request = serde_json::from_str(&line?)?;
        let groups = track::groups(&q.events, q.held_sample, &q.settings)?;
        let predictions: Vec<_> = groups
            .iter()
            .map(
                |g| match forecaster.forecast(g, q.held_sample, q.tx_centre, &q.settings) {
                    Ok(f) => serde_json::to_value(f).unwrap(),
                    Err(e) => serde_json::json!({"supported":false,"reason":e.to_string()}),
                },
            )
            .collect();
        writeln!(out, "{}", serde_json::to_string(&predictions)?)?;
        out.flush()?;
    }
    Ok(())
}
