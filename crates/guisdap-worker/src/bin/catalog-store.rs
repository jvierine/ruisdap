//! Archive metadata catalogue and bounded, persistent processing queue.
#[path = "../live_monitor.rs"]
mod live_monitor;
#[path = "../millstone_api.rs"]
mod millstone_api;
use anyhow::{ensure, Context, Result};
use axum::{
    extract::{DefaultBodyLimit, Path as ApiPath, State},
    http::{HeaderMap, StatusCode},
    routing::{get, post},
    Json, Router,
};
use chrono::{NaiveDate, Utc};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Write,
    path::{Path, PathBuf},
    process::Command,
    sync::{Arc, Mutex},
};
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    roots: Vec<Root>,
    public: PathBuf,
    state: PathBuf,
    repository: PathBuf,
    worker_threads: usize,
    max_depth: usize,
    recipes: Vec<Recipe>,
    #[serde(default)]
    calibration_root: Option<PathBuf>,
}
#[derive(Clone, Deserialize, Serialize)]
struct Root {
    label: String,
    path: PathBuf,
}
#[derive(Clone, Deserialize)]
struct Recipe {
    date: String,
    experiment: String,
    #[serde(default)]
    calibration_binding: Option<PathBuf>,
    steps: Vec<Step>,
}
#[derive(Clone, Deserialize, Serialize)]
struct Step {
    label: String,
    commands: Vec<Vec<String>>,
    completion: PathBuf,
    link: String,
}
#[derive(Clone, Serialize, Deserialize, Default)]
struct Entry {
    id: String,
    date: String,
    experiment: String,
    radar: String,
    sources: BTreeSet<String>,
    hours: BTreeSet<String>,
    products: Vec<Product>,
    supported: bool,
    #[serde(default)]
    calibration: Option<Calibration>,
}
#[derive(Clone, Serialize, Deserialize)]
struct Product {
    label: String,
    link: String,
    complete: bool,
}
#[derive(Clone, Serialize, Deserialize, Default)]
struct Inventory {
    updated: String,
    scanning: bool,
    directories: usize,
    errors: usize,
    skipped_at_depth: usize,
    sources: Vec<String>,
    experiments: Vec<Entry>,
}
#[derive(Clone, Serialize, Deserialize)]
struct Job {
    #[serde(default)]
    calibration: Option<Calibration>,
    id: String,
    experiment_id: String,
    date: String,
    experiment: String,
    state: String,
    stage: String,
    completed: usize,
    total: usize,
    cores: usize,
    created: String,
    updated: String,
    message: String,
}
#[derive(Clone)]
struct App {
    config: Arc<Config>,
    jobs: Arc<Mutex<Vec<Job>>>,
    token: Arc<String>,
    scanning: Arc<std::sync::atomic::AtomicBool>,
}
fn write_json(path: &Path, v: &impl Serialize) -> Result<()> {
    fs::create_dir_all(path.parent().context("missing parent")?)?;
    if path.starts_with("/srv/guisdap") {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(path.parent().unwrap(), fs::Permissions::from_mode(0o755))?;
    }
    let mut tmp = tempfile::NamedTempFile::new_in(path.parent().unwrap())?;
    serde_json::to_writer(&mut tmp, v)?;
    tmp.flush()?;
    tmp.as_file().sync_all()?;
    tmp.persist(path)?;
    if path.starts_with("/srv/guisdap") {
        use std::os::unix::fs::PermissionsExt;
        fs::set_permissions(path, fs::Permissions::from_mode(0o644))?;
    }
    Ok(())
}
fn id(date: &str, name: &str) -> String {
    format!("{:x}", Sha256::digest(format!("{date}/{name}")))[..24].to_owned()
}
fn date_hour(name: &str) -> Option<(String, String)> {
    let bytes = name.as_bytes();
    if bytes.len() != 11
        || bytes[8] != b'_'
        || !bytes[..8].iter().all(u8::is_ascii_digit)
        || !bytes[9..].iter().all(u8::is_ascii_digit)
    {
        return None;
    }
    let d = NaiveDate::parse_from_str(&name[..8], "%Y%m%d").ok()?;
    let hour = name[9..].parse::<u8>().ok()?;
    if hour > 23 {
        return None;
    }
    Some((d.to_string(), format!("{hour:02}:00")))
}
#[derive(Clone, Serialize, Deserialize)]
struct Calibration {
    revision: u64,
    magic_constant: f64,
    standard_error: Option<f64>,
    included: usize,
    total: usize,
    sha256: String,
}
// Plasma products use the same experiment namespace as the ion products.
// Build their configuration from this recipe, never from another experiment's default.
fn publish_plasma_config(c: &Config, r: &Recipe) -> Result<()> {
    let mut datasets = BTreeMap::new();
    for step in &r.steps {
        let Some((page, query)) = step.link.split_once('?') else {
            continue;
        };
        if page != "plasma.html" && page != "./" {
            continue;
        }
        let resolution = query
            .split('&')
            .find_map(|p| p.strip_prefix("resolution="))
            .unwrap_or("0")
            .parse::<usize>()?;
        let directory = step
            .completion
            .parent()
            .and_then(Path::file_name)
            .context("plasma export directory missing")?
            .to_string_lossy();
        datasets.insert(
            resolution,
            serde_json::json!({
                "label": step.label,
                "data_root": format!("/guisdap/data/{directory}/")
            }),
        );
    }
    if datasets.is_empty()
        || !r.steps.iter().any(|s| {
            (s.link.starts_with("plasma.html?") || s.link.starts_with("./?")) && complete(s)
        })
    {
        return Ok(());
    }
    ensure!(
        datasets.keys().copied().eq(0..datasets.len()),
        "non-contiguous plasma resolutions"
    );
    let datasets: Vec<_> = datasets.into_values().collect();
    let public = c.public.join(id(&r.date, &r.experiment));
    fs::create_dir_all(&public)?;
    write_json(
        &public.join("config.json"),
        &serde_json::json!({
            "data_root": datasets[0]["data_root"], "datasets": datasets
        }),
    )
}
fn steps(c: &Config, r: &Recipe) -> Result<(Vec<Step>, Option<Calibration>)> {
    publish_plasma_config(c, r)?;
    let Some(binding) = &r.calibration_binding else {
        return Ok((r.steps.clone(), None));
    };
    let source = fs::read_to_string(binding)?.trim().to_owned();
    let bytes = fs::read(source)?;
    let notebook: serde_json::Value = serde_json::from_slice(&bytes)?;
    let points = notebook["points"]
        .as_array()
        .context("invalid calibration notebook")?;
    let values: Vec<_> = points
        .iter()
        .filter(|p| p["included"] == true)
        .map(|p| -> Result<_> {
            let raw = p["raw_ne_proxy_m3"].as_f64().context("missing raw power")?;
            let old = p["original_calibration_multiplier"]
                .as_f64()
                .context("missing old constant")?;
            let ratio = p["estimate"]["calibration_multiplier_ratio"]
                .as_f64()
                .context("missing correction")?;
            let time = p["ion_start_unix"].as_f64().context("missing time")?;
            Ok((raw / old, raw * ratio, time.to_bits()))
        })
        .collect::<Result<_>>()?;
    let fit = guisdap_core::calibration::fit_gain(&values)
        .context("No valid included calibration points")?;
    let cal = Calibration {
        revision: notebook["revision"].as_u64().context("missing revision")?,
        magic_constant: fit.gain,
        standard_error: fit.standard_error,
        included: fit.points,
        total: points.len(),
        sha256: format!("{:x}", Sha256::digest(&bytes)),
    };
    let dir = c.state.join("recipes").join(format!(
        "{}-{}",
        id(&r.date, &r.experiment),
        &cal.sha256[..12]
    ));
    fs::create_dir_all(&dir)?;
    if !dir.join("calibration.json").exists() {
        fs::write(dir.join("calibration.json"), &bytes)?;
    }
    let mut result = r.steps.clone();
    for (i, s) in result.iter_mut().enumerate() {
        for args in &mut s.commands {
            if args.first().is_some_and(|a| a == "fit-ion") {
                let pos = args
                    .iter()
                    .position(|a| a == "--config")
                    .context("fit config missing")?
                    + 1;
                let mut v: serde_json::Value =
                    serde_json::from_slice(&fs::read(c.repository.join(&args[pos]))?)?;
                v["calibration_multiplier"] = cal.magic_constant.into();
                v["threads"] = c.worker_threads.into();
                v["calibration_provenance"] = serde_json::to_value(&cal)?;
                for field in ["output", "export"] {
                    let base = v[field].as_str().context("missing fit destination")?;
                    v[field] = format!("{base}-cal-{}", &cal.sha256[..12]).into();
                }
                s.link = format!(
                    "ion.html?experiment={}&resolution={i}",
                    id(&r.date, &r.experiment)
                );
                s.label = format!(
                    "{} · C={:.4} · revision {}",
                    s.label, cal.magic_constant, cal.revision
                );
                s.completion = PathBuf::from(v["export"].as_str().unwrap()).join("index.json");
                let config = dir.join(format!("fit-{i}.json"));
                if !config.exists() {
                    write_json(&config, &v)?;
                }
                args[pos] = config.to_string_lossy().into_owned();
            }
        }
    }
    // Per-experiment browser configuration; generated independently of live global defaults.
    let public = c.public.join(id(&r.date, &r.experiment));
    fs::create_dir_all(&public)?;
    let datasets:Vec<_>=result.iter().filter(|s|s.label.starts_with("Ion-line")).map(|s|serde_json::json!({"label":s.label,"integration_s":if s.label.contains("15 s"){15}else{60},"data_url":format!("/guisdap/data/{}/index.json",s.completion.parent().unwrap().file_name().unwrap().to_string_lossy())})).collect();
    let mut cfg: serde_json::Value =
        serde_json::from_slice(&fs::read(c.repository.join("web/ion-config.json"))?)?;
    cfg["datasets"] = datasets.into();
    if let Some(url) = cfg["instruments_url"].as_str() {
        if !url.starts_with('/') {
            cfg["instruments_url"] = format!("/guisdap/{url}").into();
        }
    }
    if result
        .iter()
        .any(|s| s.label.starts_with("Ion-line") && s.completion.is_file())
    {
        write_json(&public.join("ion-config.json"), &cfg)?;
    }
    Ok((result, Some(cal)))
}
fn complete(step: &Step) -> bool {
    if !step.completion.is_file() {
        return false;
    }
    if let Ok(bytes) = fs::read(step.completion.with_file_name("status.json")) {
        if let Ok(v) = serde_json::from_slice::<serde_json::Value>(&bytes) {
            if v.get("phase").is_some() {
                return v["phase"] == "complete";
            }
        }
    }

    let Ok(b) = fs::read(&step.completion) else {
        return false;
    };
    let Ok(v) = serde_json::from_slice::<serde_json::Value>(&b) else {
        return false;
    };
    if let Some(p) = v.get("processing") {
        return p["phase"] == "complete";
    }
    v.get("frames")
        .and_then(|a| a.as_array())
        .is_some_and(|a| !a.is_empty())
}
fn recipe<'a>(c: &'a Config, e: &Entry) -> Option<&'a Recipe> {
    c.recipes
        .iter()
        .find(|r| r.date == e.date && r.experiment == e.experiment)
}
fn add_products(c: &Config, entries: &mut BTreeMap<String, Entry>) {
    for r in &c.recipes {
        let key = id(&r.date, &r.experiment);
        let e = entries.entry(key.clone()).or_insert_with(|| Entry {
            id: key,
            date: r.date.clone(),
            experiment: r.experiment.clone(),
            radar: r.experiment.rsplit('@').next().unwrap_or("").into(),
            ..Entry::default()
        });
        match steps(c, r) {
            Ok((steps, cal)) => {
                e.supported = true;
                e.calibration = cal;
                e.products = r
                    .steps
                    .iter()
                    .chain(
                        steps
                            .iter()
                            .filter(|s| !r.steps.iter().any(|old| old.completion == s.completion)),
                    )
                    .filter(|s| s.completion.exists())
                    .map(|s| Product {
                        label: s.label.clone(),
                        link: if s.link.starts_with("./?") {
                            s.link.replacen("./?", "plasma.html?", 1)
                        } else {
                            s.link.clone()
                        },
                        complete: complete(s),
                    })
                    .collect();
            }
            Err(err) => {
                e.supported = false;
                eprintln!("recipe {}: {err}", r.experiment);
            }
        }
    }
}
fn scan(app: &App) -> Result<()> {
    let c = &app.config;
    let mut inventory = Inventory {
        scanning: true,
        sources: c.roots.iter().map(|r| r.label.clone()).collect(),
        ..Inventory::default()
    };
    let mut entries = BTreeMap::new();
    let publish = |inv: &mut Inventory, entries: &BTreeMap<String, Entry>| -> Result<()> {
        inv.updated = Utc::now().to_rfc3339();
        inv.experiments = entries.values().cloned().collect();
        inv.experiments
            .sort_by(|a, b| (&a.date, &a.experiment).cmp(&(&b.date, &b.experiment)));
        write_json(&c.public.join("catalog.json"), inv)
    };
    add_products(c, &mut entries);
    publish(&mut inventory, &entries)?;
    let mut last_publish = std::time::Instant::now();
    for root in &c.roots {
        let mut pending = vec![(root.path.clone(), 0usize, None::<String>)];
        while let Some((dir, depth, mut experiment)) = pending.pop() {
            inventory.directories += 1;
            if last_publish.elapsed().as_secs() >= 3 {
                add_products(c, &mut entries);
                publish(&mut inventory, &entries)?;
                last_publish = std::time::Instant::now();
            }
            let name = dir.file_name().unwrap_or_default().to_string_lossy();
            if name.contains('@') {
                experiment = Some(name.to_string());
            }
            if let (Some(exp), Some((date, hour))) = (&experiment, date_hour(&name)) {
                let key = id(&date, exp);
                let e = entries.entry(key.clone()).or_insert_with(|| Entry {
                    id: key,
                    date,
                    experiment: exp.clone(),
                    radar: exp.rsplit('@').next().unwrap_or("").into(),
                    ..Entry::default()
                });
                e.sources.insert(root.label.clone());
                e.hours.insert(hour);
                continue; // Do not enumerate the thousands of recordings in an hour directory.
            }
            if depth >= c.max_depth {
                inventory.skipped_at_depth += 1;
                continue;
            }
            let children = match fs::read_dir(&dir) {
                Ok(v) => v,
                Err(_) => {
                    inventory.errors += 1;
                    continue;
                }
            };
            for child in children {
                let child = match child {
                    Ok(v) => v,
                    Err(_) => {
                        inventory.errors += 1;
                        continue;
                    }
                };
                let name = child.file_name().to_string_lossy().into_owned();
                if name.starts_with('.')
                    || [
                        "lost+found",
                        "maarsy",
                        "syisr",
                        "digital_rf",
                        "raw_voltage",
                        "node_modules",
                        "target",
                        "venv",
                        "digital_rf_data",
                        "pansy_data",
                        "kaira_meteor",
                        "guisdap_processed_data",
                    ]
                    .contains(&name.as_str())
                {
                    continue;
                }
                match child.file_type() {
                    Ok(t) if t.is_dir() => {
                        pending.push((child.path(), depth + 1, experiment.clone()))
                    }
                    Err(_) => inventory.errors += 1,
                    _ => {}
                }
            }
            if inventory.directories.is_multiple_of(1000) {
                add_products(c, &mut entries);
                publish(&mut inventory, &entries)?;
            }
        }
    }
    add_products(c, &mut entries);
    inventory.scanning = false;
    publish(&mut inventory, &entries)
}
fn load_inventory(c: &Config) -> Result<Inventory> {
    Ok(serde_json::from_slice(&fs::read(
        c.public.join("catalog.json"),
    )?)?)
}
fn start_scan(app: App) {
    use std::sync::atomic::Ordering;
    if app.scanning.swap(true, Ordering::SeqCst) {
        return;
    }
    std::thread::spawn(move || {
        if let Err(e) = scan(&app) {
            eprintln!("inventory: {e:#}");
        }
        app.scanning.store(false, Ordering::SeqCst);
    });
}
fn save_jobs(app: &App, jobs: &Vec<Job>) -> Result<()> {
    write_json(&app.config.state.join("jobs.json"), jobs)?;
    write_json(&app.config.public.join("queue.json"), jobs)
}
fn authorized(h: &HeaderMap, app: &App) -> bool {
    h.get("authorization")
        .and_then(|v| v.to_str().ok())
        .is_some_and(|v| v == format!("Bearer {}", app.token))
}
#[derive(Deserialize)]
struct RecalibrationRequest {
    notebook_key: String,
    revision: u64,
}
async fn enqueue(
    State(app): State<App>,
    ApiPath(key): ApiPath<String>,
    headers: HeaderMap,
    request: Option<Json<RecalibrationRequest>>,
) -> Result<Json<Job>, (StatusCode, String)> {
    if !authorized(&headers, &app) {
        return Err((
            StatusCode::UNAUTHORIZED,
            "Enter the server processing key.".into(),
        ));
    }
    let inv = load_inventory(&app.config)
        .map_err(|e| (StatusCode::SERVICE_UNAVAILABLE, e.to_string()))?;
    let e = inv
        .experiments
        .iter()
        .find(|e| e.id == key)
        .ok_or((StatusCode::NOT_FOUND, "Experiment not in inventory".into()))?;
    let mut recipe = recipe(&app.config, e).cloned().ok_or((
        StatusCode::UNPROCESSABLE_ENTITY,
        "No validated processing recipe for this experiment/date.".into(),
    ))?;
    if let Some(Json(request)) = request {
        let prepare = (|| -> Result<PathBuf> {
            ensure!(
                request.notebook_key.len() == 64
                    && request.notebook_key.bytes().all(|b| b.is_ascii_hexdigit()),
                "Invalid calibration notebook key"
            );
            let root = app
                .config
                .calibration_root
                .clone()
                .unwrap_or_else(|| app.config.state.parent().unwrap().join("calibrations"));
            let bytes = fs::read(root.join(&request.notebook_key).join("latest.json"))?;
            let notebook: serde_json::Value = serde_json::from_slice(&bytes)?;
            ensure!(
                notebook["revision"].as_u64() == Some(request.revision),
                "Calibration changed; refresh and review the latest revision"
            );
            let points = notebook["points"]
                .as_array()
                .context("Missing calibration points")?;
            ensure!(
                points.iter().any(|p| p["included"] == true),
                "No included calibration points"
            );
            for p in points.iter().filter(|p| p["included"] == true) {
                let time = p["ion_start_unix"]
                    .as_f64()
                    .context("Missing calibration time")?;
                ensure!(
                    chrono::DateTime::from_timestamp(time as i64, 0)
                        .context("Invalid calibration time")?
                        .format("%Y-%m-%d")
                        .to_string()
                        == e.date,
                    "Calibration point belongs to another UTC date"
                );
            }
            let sha = format!("{:x}", Sha256::digest(&bytes));
            let directory = app.config.state.join("reviewed-calibrations");
            fs::create_dir_all(&directory)?;
            let snapshot = directory.join(format!("{sha}.json"));
            if !snapshot.exists() {
                fs::write(&snapshot, &bytes)?;
            }
            let binding = directory.join(format!("{sha}.binding"));
            fs::write(&binding, snapshot.to_string_lossy().as_bytes())?;
            Ok(binding)
        })()
        .map_err(|e| (StatusCode::UNPROCESSABLE_ENTITY, e.to_string()))?;
        recipe.calibration_binding = Some(prepare);
        // Recalibration reuses plasma products and queues the ion fitting steps.
        recipe.steps.retain(|s| {
            s.commands
                .iter()
                .any(|args| args.first().is_some_and(|a| a == "fit-ion"))
        });
    }
    let mut jobs = app.jobs.lock().unwrap();
    if let Some(j) = jobs
        .iter()
        .find(|j| j.experiment_id == key && ["queued", "running"].contains(&j.state.as_str()))
    {
        return Ok(Json(j.clone()));
    }
    let (steps, calibration) = steps(&app.config, &recipe)
        .map_err(|e| (StatusCode::UNPROCESSABLE_ENTITY, e.to_string()))?;
    if steps.iter().all(complete) {
        return Err((
            StatusCode::CONFLICT,
            "All configured products are already complete.".into(),
        ));
    }
    let now = Utc::now().to_rfc3339();
    let job = Job {
        calibration,
        id: format!("{}-{}", key, Utc::now().timestamp_millis()),
        experiment_id: key,
        date: e.date.clone(),
        experiment: e.experiment.clone(),
        state: "queued".into(),
        stage: "Waiting for compute slot".into(),
        completed: 0,
        total: steps.len(),
        cores: app.config.worker_threads,
        created: now.clone(),
        updated: now,
        message: String::new(),
    };
    write_json(
        &app.config.state.join(format!("{}.steps.json", job.id)),
        &steps,
    )
    .map_err(|e| (StatusCode::INTERNAL_SERVER_ERROR, e.to_string()))?;
    jobs.push(job.clone());
    save_jobs(&app, &jobs).map_err(|e| (StatusCode::INTERNAL_SERVER_ERROR, e.to_string()))?;
    Ok(Json(job))
}
async fn rescan(State(app): State<App>, headers: HeaderMap) -> Result<StatusCode, StatusCode> {
    if !authorized(&headers, &app) {
        return Err(StatusCode::UNAUTHORIZED);
    }
    start_scan(app);
    Ok(StatusCode::ACCEPTED)
}
fn update(app: &App, key: &str, modify: impl FnOnce(&mut Job)) -> Result<()> {
    let mut jobs = app.jobs.lock().unwrap();
    let j = jobs
        .iter_mut()
        .find(|j| j.id == key)
        .context("missing job")?;
    modify(j);
    j.updated = Utc::now().to_rfc3339();
    save_jobs(app, &jobs)
}
fn run_job(app: &App, job: &Job) -> Result<()> {
    let c = &app.config;
    let steps: Vec<Step> =
        serde_json::from_slice(&fs::read(c.state.join(format!("{}.steps.json", job.id)))?)?;
    let log = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(c.state.join(format!("{}.log", job.id)))?;
    for (i, step) in steps.iter().enumerate() {
        update(app, &job.id, |j| {
            j.stage = step.label.clone();
            j.completed = i;
        })?;
        if !complete(step) {
            for args in &step.commands {
                ensure!(!args.is_empty(), "empty command");
                let status = Command::new(c.repository.join("target/release/guisdap-worker"))
                    .args(args)
                    .current_dir(&c.repository)
                    .env("RAYON_NUM_THREADS", c.worker_threads.to_string())
                    .stdout(log.try_clone()?)
                    .stderr(log.try_clone()?)
                    .status()?;
                ensure!(
                    status.success(),
                    "{} failed ({status}); inspect the server job log",
                    step.label
                );
            }
            ensure!(
                complete(step),
                "{} exited without a complete product",
                step.label
            );
        }
        update(app, &job.id, |j| j.completed = i + 1)?;
    }
    Ok(())
}
fn worker(app: App) {
    std::thread::spawn(move || loop {
        let next = {
            let mut jobs = app.jobs.lock().unwrap();
            let j = jobs.iter_mut().find(|j| j.state == "queued");
            if let Some(j) = j {
                j.state = "running".into();
                j.updated = Utc::now().to_rfc3339();
                let next = Some(j.clone());
                if let Err(e) = save_jobs(&app, &jobs) {
                    eprintln!("queue write: {e}");
                    None
                } else {
                    next
                }
            } else {
                None
            }
        };
        if let Some(j) = next {
            let result = run_job(&app, &j);
            let _ = update(&app, &j.id, |j| match result {
                Ok(()) => {
                    j.state = "complete".into();
                    j.stage = "Complete".into();
                }
                Err(e) => {
                    j.state = "failed".into();
                    j.message = e.to_string();
                }
            });
            // Refresh availability without rescanning the archive.
            if let Ok(mut inv) = load_inventory(&app.config) {
                let mut entries = inv
                    .experiments
                    .into_iter()
                    .map(|e| (e.id.clone(), e))
                    .collect();
                add_products(&app.config, &mut entries);
                inv.experiments = entries.into_values().collect();
                let _ = write_json(&app.config.public.join("catalog.json"), &inv);
            }
        } else {
            std::thread::sleep(std::time::Duration::from_secs(2));
        }
    });
}
fn monitor(app: App) {
    std::thread::spawn(move || loop {
        let publish = || -> Result<()> {
            let mut jobs = app.jobs.lock().unwrap().clone();
            for job in &mut jobs {
                if job.state != "running" {
                    continue;
                }
                let steps: Vec<Step> = serde_json::from_slice(&fs::read(
                    app.config.state.join(format!("{}.steps.json", job.id)),
                )?)?;
                if let Some(step) = steps.get(job.completed) {
                    if let Ok(bytes) = fs::read(step.completion.with_file_name("status.json")) {
                        if let Ok(v) = serde_json::from_slice::<serde_json::Value>(&bytes) {
                            job.stage = format!(
                                "{} · {} / {} intervals",
                                step.label, v["completed"], v["total"]
                            );
                        }
                    }
                }
            }
            write_json(&app.config.public.join("queue.json"), &jobs)?;
            // Switch each default resolution only after its calibrated index exists.
            for r in &app.config.recipes {
                let (steps, _) = steps(&app.config, r)?;
                let mut cfg: serde_json::Value = serde_json::from_slice(&fs::read(
                    app.config.repository.join("web/ion-config.json"),
                )?)?;
                let mut changed = false;
                for (i, s) in steps
                    .iter()
                    .filter(|s| s.label.starts_with("Ion-line"))
                    .enumerate()
                {
                    if s.completion.exists() {
                        let dt = if s.label.contains("15 s") { 15 } else { 60 };
                        let native = s.label.contains("native");
                        if let Some(datasets) = cfg["datasets"].as_array_mut() {
                            for d in datasets {
                                if d["integration_s"] == dt
                                    && d["label"].as_str().unwrap_or("").contains("native")
                                        == native
                                {
                                    d["data_url"] = format!(
                                        "data/{}/index.json",
                                        s.completion
                                            .parent()
                                            .unwrap()
                                            .file_name()
                                            .unwrap()
                                            .to_string_lossy()
                                    )
                                    .into();
                                    changed = true;
                                }
                            }
                        }
                        let _ = i;
                    }
                }
                if changed {
                    write_json(&app.config.public.join("default-ion-config.json"), &cfg)?;
                    // Keep older browser builds on the same calibrated default.
                    write_json(
                        &app.config
                            .public
                            .parent()
                            .context("catalog public parent")?
                            .join("ion-config.json"),
                        &cfg,
                    )?;
                }
            }
            if !app.scanning.load(std::sync::atomic::Ordering::SeqCst) {
                let mut inv = load_inventory(&app.config)?;
                let mut entries = inv
                    .experiments
                    .into_iter()
                    .map(|e| (e.id.clone(), e))
                    .collect();
                add_products(&app.config, &mut entries);
                inv.experiments = entries.into_values().collect();
                inv.updated = Utc::now().to_rfc3339();
                write_json(&app.config.public.join("catalog.json"), &inv)?;
            }
            Ok(())
        };
        if let Err(e) = publish() {
            eprintln!("queue monitor: {e}");
        }
        std::thread::sleep(std::time::Duration::from_secs(5));
    });
}

#[tokio::main(worker_threads = 2)]
async fn main() -> Result<()> {
    let path = std::env::args()
        .nth(1)
        .context("usage: catalog-store CONFIG")?;
    let c: Config = serde_json::from_slice(&fs::read(path)?)?;
    ensure!(
        c.worker_threads > 0 && c.worker_threads <= 128,
        "worker_threads must be 1..128"
    );
    fs::create_dir_all(&c.state)?;
    fs::create_dir_all(&c.public)?;
    let token = fs::read_to_string(c.state.join("processing-key"))?
        .trim()
        .to_owned();
    ensure!(token.len() >= 32, "processing key too short");
    let mut jobs: Vec<Job> = match fs::read(c.state.join("jobs.json")) {
        Ok(b) => serde_json::from_slice(&b)?,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => vec![],
        Err(e) => return Err(e.into()),
    };
    for j in &mut jobs {
        if j.state == "running" {
            j.state = "interrupted".into();
            j.message = "Service restarted; inspect outputs before retrying".into();
        }
    }
    let app = App {
        config: Arc::new(c),
        jobs: Arc::new(Mutex::new(jobs)),
        token: Arc::new(token),
        scanning: Arc::new(std::sync::atomic::AtomicBool::new(false)),
    };
    save_jobs(&app, &app.jobs.lock().unwrap())?;
    start_scan(app.clone());
    worker(app.clone());
    monitor(app.clone());
    let millstone_routes = millstone_api::router(
        app.config.state.clone(),
        app.config.public.parent().unwrap().to_path_buf(),
        app.config
            .calibration_root
            .clone()
            .unwrap_or_else(|| app.config.state.join("calibration")),
    );
    let router = Router::new()
        .route("/guisdap/api/catalog/jobs/{id}", post(enqueue))
        .route("/guisdap/api/catalog/rescan", post(rescan))
        .route("/guisdap/api/catalog/health", get(|| async { "ok" }))
        .layer(DefaultBodyLimit::max(1024))
        .with_state(app)
        .merge(live_monitor::router("/guisdap/api/catalog/live"))
        .merge(millstone_routes);
    axum::serve(
        tokio::net::TcpListener::bind("127.0.0.1:8092").await?,
        router,
    )
    .await?;
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dates_are_calendar_dates_not_arbitrary_digits() {
        assert_eq!(
            date_hour("20221102_17"),
            Some(("2022-11-02".into(), "17:00".into()))
        );
        assert!(date_hour("20220230_12").is_none());
        assert!(date_hour("20221102_24").is_none());
        assert!(date_hour("20221102_17_extra").is_none());
        assert!(date_hour("é20221102_1").is_none());
    }
    #[test]
    fn duplicate_sources_share_identity_but_versions_do_not() {
        assert_eq!(id("2022-11-02", "beata@uhf"), id("2022-11-02", "beata@uhf"));
        assert_ne!(
            id("2022-11-02", "beata@uhf"),
            id("2022-11-02", "beata2@uhf")
        );
    }
    #[test]
    fn archive_scan_deduplicates_sources_and_stops_before_record_files() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path();
        for disk in ["a", "b"] {
            fs::create_dir_all(root.join(disk).join("beata@uhf/20221102_17/do-not-descend"))
                .unwrap();
        }
        let c = Config {
            roots: vec![
                Root {
                    label: "a".into(),
                    path: root.join("a"),
                },
                Root {
                    label: "b".into(),
                    path: root.join("b"),
                },
                Root {
                    label: "offline".into(),
                    path: root.join("absent"),
                },
            ],
            public: root.join("public"),
            state: root.join("state"),
            repository: root.into(),
            worker_threads: 128,
            max_depth: 8,
            recipes: vec![],
            calibration_root: None,
        };
        let app = App {
            config: Arc::new(c),
            jobs: Arc::new(Mutex::new(vec![])),
            token: Arc::new("secret".into()),
            scanning: Arc::new(std::sync::atomic::AtomicBool::new(false)),
        };
        scan(&app).unwrap();
        let inv = load_inventory(&app.config).unwrap();
        assert!(!inv.scanning);
        assert_eq!(inv.errors, 1);
        assert_eq!(inv.experiments.len(), 1);
        assert_eq!(inv.experiments[0].sources.len(), 2);
        assert_eq!(inv.experiments[0].hours.len(), 1);
        assert_eq!(inv.directories, 7);
        assert!(!authorized(&HeaderMap::new(), &app));
    }
    #[test]
    fn plasma_config_is_published_for_its_experiment_without_calibration() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path();
        let mut c = Config {
            roots: vec![],
            public: root.join("public"),
            state: root.join("state"),
            repository: root.into(),
            worker_threads: 1,
            max_depth: 8,
            recipes: vec![],
            calibration_root: None,
        };
        let mut r = Recipe {
            date: "2022-11-02".into(),
            experiment: "beata@uhf".into(),
            calibration_binding: None,
            steps: vec![],
        };
        for (resolution, cadence) in [(1, 60), (0, 15)] {
            let export = root.join(format!("plasma-{cadence}s"));
            fs::create_dir_all(&export).unwrap();
            r.steps.push(Step {
                label: format!("Plasma-line {cadence} s"),
                commands: vec![],
                completion: export.join("index.json"),
                link: format!("plasma.html?resolution={resolution}"),
            });
        }
        let config = c
            .public
            .join(id(&r.date, &r.experiment))
            .join("config.json");
        steps(&c, &r).unwrap();
        assert!(!config.exists()); // No published observation until data exists.
        write_json(
            &r.steps[0].completion,
            &serde_json::json!({"frames":[{"id":"one"}]}),
        )
        .unwrap();
        c.recipes.push(r.clone());
        let mut entries = BTreeMap::new();
        add_products(&c, &mut entries);
        let v: serde_json::Value = serde_json::from_slice(&fs::read(config).unwrap()).unwrap();
        assert_eq!(v["datasets"][0]["data_root"], "/guisdap/data/plasma-15s/");
        assert_eq!(v["datasets"][1]["data_root"], "/guisdap/data/plasma-60s/");
        assert_eq!(v["data_root"], v["datasets"][0]["data_root"]);
        assert!(!c.public.join("config.json").exists());
        assert!(entries[&id(&r.date, &r.experiment)]
            .products
            .iter()
            .any(|p| p.complete));
    }
    #[test]
    fn partial_products_do_not_count_as_complete() {
        let tmp = tempfile::tempdir().unwrap();
        let index = tmp.path().join("index.json");
        let s = Step {
            label: "test".into(),
            commands: vec![],
            completion: index.clone(),
            link: String::new(),
        };
        assert!(!complete(&s));
        write_json(
            &index,
            &serde_json::json!({"processing":{"phase":"running"}}),
        )
        .unwrap();
        assert!(!complete(&s));
        write_json(
            &index,
            &serde_json::json!({"processing":{"phase":"complete"}}),
        )
        .unwrap();
        assert!(complete(&s));
        write_json(&index, &serde_json::json!({"frames":[]})).unwrap();
        assert!(!complete(&s));
    }
}
