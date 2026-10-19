use crate::{doc, el, json, plotting, text};
use serde::Deserialize;
use std::{cell::Cell, rc::Rc};
use wasm_bindgen::{prelude::*, JsCast};
use wasm_bindgen_futures::spawn_local;
#[derive(Deserialize)]
struct Rank {
    #[serde(default)]
    host: String,
    rank: usize,
    ranks: usize,
    phase: String,
    completed: usize,
    total: usize,
    elapsed_s: f64,
    updated_unix_s: f64,
}
#[derive(Deserialize)]
struct Preview {
    rank: usize,
    updated_unix_s: f64,
    station: String,
    antenna: String,
    mode: String,
    first_lag_real: Vec<Option<f64>>,
    retained_fraction: Vec<f64>,
    cancelled: usize,
    failed_removals: usize,
}
#[derive(Deserialize)]
struct Snapshot {
    server_unix_s: f64,
    ranks: Vec<Rank>,
    previews: Vec<Preview>,
    legacy: Vec<serde_json::Value>,
}
fn escaped(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}
async fn poll() -> Result<(), JsValue> {
    let s: Snapshot = json("api/catalog/live").await?;
    let mut html = String::new();
    let total: usize = s.ranks.iter().map(|r| r.total).sum();
    let done: usize = s.ranks.iter().map(|r| r.completed).sum();
    text(
        "live-status",
        &format!(
            "{done} / {total} integration periods · {} reporting ranks",
            s.ranks.len()
        ),
    );
    for r in s.ranks {
        html.push_str(&format!("<article class='plot'><b>{} · rank {} / {}</b><p>{}</p><progress value='{}' max='{}'></progress><p>{} / {} periods · {:.0} s elapsed</p><small>Updated {:.0} s ago{}</small></article>",escaped(&r.host),r.rank,r.ranks,escaped(&r.phase),r.completed,r.total.max(1),r.completed,r.total,r.elapsed_s,(s.server_unix_s-r.updated_unix_s).max(0.),if s.server_unix_s-r.updated_unix_s>30.{" · stale"}else{""}));
    }
    el("live-ranks").set_inner_html(&html);
    if let Some(p) = s
        .previews
        .iter()
        .max_by(|a, b| a.updated_unix_s.total_cmp(&b.updated_unix_s))
    {
        let max = p
            .first_lag_real
            .iter()
            .flatten()
            .map(|x| x.abs())
            .fold(1., f64::max);
        let n = p.first_lag_real.len();
        let points: Vec<_> = p
            .first_lag_real
            .iter()
            .enumerate()
            .flat_map(|(i, v)| [Some(i as f64), *v])
            .collect();
        plotting::shared_curve("live-preview",&serde_json::json!({"title":"Latest first-lag profile","bounds":[0.,n.max(1) as f64,-max*1.1,max*1.1],"xlabel":"Range gate index","ylabel":"First-lag scattering strength / K","series":[{"xy":points,"radius":1.5}]}).to_string());
        let retention =
            p.retained_fraction.iter().sum::<f64>() / p.retained_fraction.len().max(1) as f64;
        text("live-preview-context",&format!("{} · {} · {} · rank {} · age {:.0} s · {:.1}% retained · {} accepted cancellations, {} failed removals",p.station,p.antenna,p.mode,p.rank,(s.server_unix_s-p.updated_unix_s).max(0.),100.*retention,p.cancelled,p.failed_removals));
    } else {
        text("live-preview", "Waiting for the next integration preview.");
    }
    let mut legacy = String::new();
    for v in s.legacy {
        let phase = v["phase"]
            .as_str()
            .or_else(|| v["stage"].as_str())
            .unwrap_or("processing");
        let state = v["state"].as_str().unwrap_or("");
        let date = v["date"].as_str().unwrap_or("");
        let experiment = v["experiment"].as_str().unwrap_or("EISCAT processing");
        let cores = v["cores"].as_u64().unwrap_or(0);
        let done = v["completed"].as_u64().unwrap_or(0);
        let total = v["total"].as_u64().unwrap_or(0);
        legacy.push_str(&format!("<article class='plot'><b>Avaruus · {}</b><p>{} · {}</p><p>{} · {done} / {total} stages · {cores} threads</p><progress value='{done}' max='{}'></progress></article>",escaped(experiment),escaped(date),escaped(state),escaped(phase),total.max(1)));
    }
    el("live-legacy").set_inner_html(&legacy);
    Ok(())
}
#[wasm_bindgen]
pub async fn start_live() -> Result<(), JsValue> {
    crate::navigation::setup()?;
    let enabled = Rc::new(Cell::new(true));
    let busy = Rc::new(Cell::new(false));
    let active = enabled.clone();
    let button = Closure::<dyn FnMut()>::new(move || {
        active.set(!active.get());
        text(
            "live-toggle",
            if active.get() {
                "Pause monitoring"
            } else {
                "Resume monitoring"
            },
        );
        if !active.get() {
            text("live-status", "Monitoring paused; processing continues.");
        }
    });
    el("live-toggle").add_event_listener_with_callback("click", button.as_ref().unchecked_ref())?;
    button.forget();
    let callback = Closure::<dyn FnMut()>::new(move || {
        if !enabled.get() || doc().hidden() || busy.get() {
            return;
        }
        busy.set(true);
        let busy = busy.clone();
        spawn_local(async move {
            if let Err(e) = poll().await {
                text(
                    "live-status",
                    &format!(
                        "Progress unavailable: {}",
                        e.as_string().unwrap_or_default()
                    ),
                )
            }
            busy.set(false);
        });
    });
    web_sys::window()
        .unwrap()
        .set_interval_with_callback_and_timeout_and_arguments_0(
            callback.as_ref().unchecked_ref(),
            2000,
        )?;
    callback.forget();
    poll().await
}
