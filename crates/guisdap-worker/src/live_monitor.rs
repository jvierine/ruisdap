//! The lease is the only trigger for workers to export detailed live arrays.
use axum::{extract::State, http::header, routing::get, Json, Router};
use serde_json::{json, Value};
use std::{
    path::PathBuf,
    sync::Arc,
    time::{SystemTime, UNIX_EPOCH},
};
#[derive(Clone)]
struct Store {
    roots: Vec<PathBuf>,
    legacy: Vec<PathBuf>,
}
fn now() -> f64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_secs_f64()
}
fn paths(name: &str) -> Vec<PathBuf> {
    std::env::var(name)
        .unwrap_or_default()
        .split(':')
        .filter(|s| !s.is_empty())
        .map(PathBuf::from)
        .collect()
}
async fn probe(State(s): State<Arc<Store>>) -> impl axum::response::IntoResponse {
    let expires = now() + 12.;
    let bytes = serde_json::to_vec(&json!({"expires_unix_s":expires})).unwrap();
    let mut ranks = Vec::new();
    let mut previews = Vec::new();
    let mut unavailable = 0;
    for root in &s.roots {
        let lease = (|| -> std::io::Result<()> {
            let tmp = tempfile::NamedTempFile::new_in(root)?;
            std::fs::write(tmp.path(), &bytes)?;
            // The worker may use a different numeric UID on the shared mount.
            // Only the public expiry timestamp is exposed, never a secret.
            #[cfg(unix)]
            {
                use std::os::unix::fs::PermissionsExt;
                std::fs::set_permissions(tmp.path(), std::fs::Permissions::from_mode(0o644))?;
            }
            tmp.persist(root.join("lease.json"))?;
            Ok(())
        })();
        if lease.is_err() {
            unavailable += 1;
            continue;
        }
        let Ok(entries) = std::fs::read_dir(root) else {
            unavailable += 1;
            continue;
        };
        for entry in entries.flatten().take(4096) {
            let name = entry.file_name().to_string_lossy().into_owned();
            if !((name.starts_with("rank-") || name.starts_with("preview-"))
                && name.ends_with(".json"))
            {
                continue;
            }
            let Ok(meta) = entry.metadata() else {
                continue;
            };
            if !meta.is_file() || meta.len() > 2 * 1024 * 1024 {
                continue;
            }
            if let Ok(bytes) = std::fs::read(entry.path()) {
                if let Ok(value) = serde_json::from_slice::<Value>(&bytes) {
                    if name.starts_with("rank-") {
                        ranks.push(value);
                    } else {
                        previews.push(value);
                    }
                }
            }
        }
    }
    let mut legacy = Vec::new();
    for path in &s.legacy {
        let Ok(meta) = std::fs::metadata(path) else {
            continue;
        };
        if meta.len() > 1024 * 1024 {
            continue;
        }
        if let Ok(bytes) = std::fs::read(path) {
            if let Ok(value) = serde_json::from_slice::<Value>(&bytes) {
                if let Value::Array(items) = value {
                    legacy.extend(items);
                } else {
                    legacy.push(value);
                }
            }
        }
    }
    (
        [(header::CACHE_CONTROL, "no-store")],
        Json(
            json!({"schema":"ruisdap-live-v1","server_unix_s":now(),"lease_expires_unix_s":expires,"ranks":ranks,"previews":previews,"legacy":legacy,"unavailable_resources":unavailable}),
        ),
    )
}
pub fn router(route: &str) -> Router {
    let mut roots = paths("RUISDAP_LIVE_ROOTS");
    if roots.is_empty() {
        roots = paths("RUISDAP_LIVE_ROOT");
    }
    Router::new()
        .route(route, get(probe))
        .with_state(Arc::new(Store {
            roots,
            legacy: paths("RUISDAP_LEGACY_STATUS"),
        }))
}
