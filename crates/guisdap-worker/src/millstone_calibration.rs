// Density calibration is a separate immutable product. Temperatures and velocity
// retain their source fit; this operation does not claim a new ACF fit.
#[derive(Deserialize)]
struct CalibrationRequest {
    notebook_key: String,
    revision: u64,
    date: String,
    antenna: String,
}
fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    (v[(v.len() - 1) / 2] + v[v.len() / 2]) / 2.
}
fn ion_proxy(ne: f64, te: f64, ti: f64, k: f64) -> f64 {
    let a = k * k * 8.8541878188e-12 * 1.380649e-23 * te / (1.602176634e-19_f64).powi(2);
    let q = a / ne;
    ne * (1. + te / ti) / ((1. + q) * (1. + q + te / ti))
}
fn density(proxy: f64, te: f64, ti: f64, k: f64) -> f64 {
    if ![proxy, te, ti].iter().all(|v| v.is_finite() && *v > 0.) {
        return f64::NAN;
    }
    let (mut lo, mut hi) = (0., proxy.max(1.));
    while ion_proxy(hi, te, ti, k) < proxy {
        hi *= 2.;
    }
    for _ in 0..80 {
        let mid = (lo + hi) / 2.;
        if ion_proxy(mid, te, ti, k) < proxy {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo + hi) / 2.
}
fn atomic_json(path: &std::path::Path, v: &Value) -> anyhow::Result<()> {
    use std::io::Write;
    let mut f = tempfile::NamedTempFile::new_in(path.parent().unwrap())?;
    f.write_all(&serde_json::to_vec(v)?)?;
    f.as_file().sync_all()?;
    f.persist(path)?;
    Ok(())
}
fn calibrate_density(s: &Store, q: CalibrationRequest) -> anyhow::Result<Value> {
    use sha2::{Digest, Sha256};
    anyhow::ensure!(
        q.notebook_key.len() == 64 && q.notebook_key.bytes().all(|v| v.is_ascii_hexdigit()),
        "Invalid notebook key"
    );
    let notebook = std::fs::read(s.calibration_root.join(&q.notebook_key).join("latest.json"))?;
    let review: Value = serde_json::from_slice(&notebook)?;
    anyhow::ensure!(
        review["revision"].as_u64() == Some(q.revision),
        "Calibration revision changed; reload before recalculating"
    );
    let ds = s
        .index
        .datasets
        .get(&format!("{}/{}", q.date, q.antenna))
        .ok_or_else(|| anyhow::anyhow!("Unknown dataset"))?;
    let source = std::fs::read(&ds.manifest)?;
    let mut m: Value = serde_json::from_slice(&source)?;
    let k = 4. * std::f64::consts::PI * m["radar_frequency_hz"].as_f64().unwrap() / 299792458.;
    let fits = m["fits"].as_array().unwrap();
    let mut groups: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut checked = Vec::new();
    for p in review["points"]
        .as_array()
        .ok_or_else(|| anyhow::anyhow!("No points"))?
        .iter()
        .filter(|p| p["included"] == true)
    {
        anyhow::ensure!(
            p["antenna"] == q.antenna && p["date"] == q.date,
            "Calibration mixes datasets"
        );
        let source = s
            .index
            .plasma
            .get(p["observation_id"].as_str().unwrap_or(""))
            .ok_or_else(|| anyhow::anyhow!("Unknown plasma record"))?;
        anyhow::ensure!(
            source.antenna == q.antenna && source.date == q.date,
            "Calibration mixes instruments"
        );
        let fit = fits
            .iter()
            .find(|f| f["id"] == p["fit_id"])
            .ok_or_else(|| anyhow::anyhow!("Unknown ion fit"))?;
        let dt = ((fit["start_unix_s"].as_f64().unwrap() + fit["end_unix_s"].as_f64().unwrap()
            - source.start
            - source.end)
            / 2.)
            .abs();
        anyhow::ensure!(dt <= 300., "Ion/plasma pair is more than 300 seconds apart");
        let gate = p["gate"]
            .as_u64()
            .ok_or_else(|| anyhow::anyhow!("Missing gate"))? as usize;
        let ranges = m["range_m"].as_array().unwrap();
        anyhow::ensure!(gate < ranges.len(), "Gate outside product");
        anyhow::ensure!(
            (ranges[gate].as_f64().unwrap() - p["range_m"].as_f64().unwrap_or(f64::NAN)).abs()
                < 0.01,
            "Gate and range disagree"
        );
        let row = &fit["parameters"][gate];
        let raw = row[0].as_f64().unwrap_or(f64::NAN);
        let ti = row[1].as_f64().unwrap_or(f64::NAN);
        let te = row[2].as_f64().unwrap_or(f64::NAN);
        let frequency = p["frequency_hz"].as_f64().unwrap_or(f64::NAN);
        let sigma = p["frequency_standard_error_hz"]
            .as_f64()
            .unwrap_or(f64::NAN);
        anyhow::ensure!(
            [raw, te, ti, sigma]
                .iter()
                .all(|v| v.is_finite() && *v > 0.)
                && frequency.is_finite()
                && (1e6..=12.5e6).contains(&frequency.abs()),
            "Nonphysical calibration point"
        );
        let omega2 = (std::f64::consts::TAU * frequency).powi(2)
            - 3. * k * k * 1.380649e-23 * te / 9.1093837139e-31;
        anyhow::ensure!(omega2 > 0., "Thermal correction exceeds frequency squared");
        let ne = omega2 * 8.8541878188e-12 * 9.1093837139e-31 / (1.602176634e-19_f64).powi(2);
        let gain = ion_proxy(ne, te, ti, k) / raw;
        groups
            .entry(fit["id"].as_str().unwrap().into())
            .or_default()
            .push(gain);
        let mut checked_point = p.clone();
        checked_point["multiplier"] = gain.into();
        checked_point["plasma_density_m3"] = ne.into();
        checked.push(checked_point);
    }
    anyhow::ensure!(!groups.is_empty(), "Include at least one reviewed point");
    let independent: Vec<_> = groups.values().map(|v| median(v.clone())).collect();
    let gain = median(independent.clone());
    anyhow::ensure!(gain.is_finite() && gain > 0., "Invalid calibration gain");
    let scatter = median(independent.iter().map(|v| (v - gain).abs()).collect()) * 1.4826;
    let se = if independent.len() > 1 {
        Some(scatter / (independent.len() as f64).sqrt())
    } else {
        None
    };
    let hash = format!(
        "{:x}",
        Sha256::digest([source.as_slice(), notebook.as_slice()].concat())
    );
    let name = format!(
        "millstone-{}-{}-cal-{}",
        q.date,
        q.antenna.to_lowercase(),
        &hash[..12]
    );
    let dest = s.public.join("data").join(&name);
    if !dest.exists() {
        let temp = tempfile::tempdir_in(s.public.join("data"))?;
        let out = temp.path();
        let nt = m["fits"].as_array().unwrap().len();
        let nr = m["range_m"].as_array().unwrap().len();
        let mut parameters = vec![f32::NAN; 4 * nt * nr];
        let mut original_parameters: Vec<f64> = Vec::new();
        let mut calibrated = Vec::new();
        for (t, fit) in m["fits"].as_array_mut().unwrap().iter_mut().enumerate() {
            let mut proxy = Vec::new();
            for (r, row) in fit["parameters"]
                .as_array_mut()
                .unwrap()
                .iter_mut()
                .enumerate()
            {
                let values: Vec<_> = row
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap_or(f64::NAN))
                    .collect();
                original_parameters.extend(values.iter().copied());
                proxy.push(row[0].clone());
                let ne = density(values[0] * gain, values[2], values[1], k);
                row[0] = serde_json::to_value(ne)?;
                calibrated.push(ne);
                for p in 0..4 {
                    parameters[(p * nr + r) * nt + t] =
                        if p == 0 { ne as f32 } else { values[p] as f32 };
                }
            }
            fit["uncalibrated_density_proxy"] = proxy.into();
            if let Some(errors) = fit["standard_error"].as_array().cloned() {
                fit["uncalibrated_standard_error"] = errors.clone().into();
                let corrected: Vec<Value> = errors
                    .iter()
                    .enumerate()
                    .map(|(ri, errors)| {
                        let mut e = errors.clone();
                        let ne = fit["parameters"][ri][0].as_f64().unwrap_or(f64::NAN);
                        let te = fit["parameters"][ri][2].as_f64().unwrap_or(f64::NAN);
                        let ti = fit["parameters"][ri][1].as_f64().unwrap_or(f64::NAN);
                        let raw = fit["uncalibrated_density_proxy"][ri]
                            .as_f64()
                            .unwrap_or(f64::NAN);
                        let a = k * k * 8.8541878188e-12 * 1.380649e-23 * te
                            / (1.602176634e-19_f64).powi(2);
                        let ratio = 1. + te / ti;
                        let derivative = ion_proxy(ne, te, ti, k)
                            * (3. / ne - 1. / (ne + a) - ratio / (ratio * ne + a));
                        let sigma = errors[0].as_f64().map(|sigma| {
                            ((gain * sigma).powi(2) + (raw * se.unwrap_or(0.)).powi(2)).sqrt()
                                / derivative
                        });
                        e[0] = sigma
                            .and_then(|v| v.is_finite().then_some(v))
                            .map_or(Value::Null, Value::from);
                        e
                    })
                    .collect();
                fit["standard_error"] = corrected.into();
            }
        }
        let relative = ds
            .manifest
            .parent()
            .unwrap()
            .strip_prefix(s.public.join("data"))?
            .to_string_lossy();
        for field in ["raw_url", "plasma_url"] {
            let url = m[field].as_str().unwrap();
            m[field] = format!("/guisdap/data/{relative}/{url}").into();
        }
        m["parameters_url"] = "parameters.f32".into();
        m["parameter_names"][0] = "Ne / m⁻³".into();
        m["density_semantics"] = "plasma_line_calibrated_density".into();
        m["calibration"] = json!({"gain":gain,"revision":q.revision,"independent_fit_count":independent.len(),"gain_standard_error":se,"density_uncertainty":"Conditional source-amplitude and independent-fit gain error; Ti/Te covariance unavailable","method":"Bohm-Gross plasma density; finite-Debye integrated ion power. Density recomputed with source Ti/Te fixed; no new ACF fit.","sha256":hash});
        std::fs::write(
            out.join("parameters.f32"),
            parameters
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect::<Vec<_>>(),
        )?;
        let h = hdf5::File::create_excl(out.join("calibration.h5"))?;
        h.new_dataset_builder()
            .with_data(notebook.as_slice())
            .create("review_json_utf8")?;
        h.new_dataset_builder()
            .with_data(source.as_slice())
            .create("source_manifest_json_utf8")?;
        h.new_dataset_builder()
            .with_data(serde_json::to_vec(&checked)?.as_slice())
            .create("checked_points_json_utf8")?;
        h.new_dataset::<f64>()
            .shape([nt, nr, 4])
            .create("source_fit_parameters")?
            .write_raw(&original_parameters)?;
        h.new_dataset::<f64>()
            .shape([nt, nr])
            .create("calibrated_ne_m3")?
            .write_raw(&calibrated)?;
        h.new_dataset_builder()
            .with_data(&[gain])
            .create("calibration_gain")?; // Calibrated fits use the same interchange records as both source adapters.
        h.new_attr_builder()
            .with_data("ruisdap-isr-v1")
            .create("schema")?;
        h.new_attr_builder()
            .with_data("UTC Unix seconds; half-open intervals")
            .create("time_convention")?;
        let records = h.create_group("records")?;
        let sources = h.create_group("sources")?;
        let calibration = h.create_group("calibration")?;
        let sid = &hash[..20];
        let source_group = sources.create_group(sid)?;
        source_group
            .new_dataset_builder()
            .with_data(source.as_slice())
            .create("source_manifest_json_utf8")?;
        let calibration_json = serde_json::to_vec(&m["calibration"])?;
        calibration
            .new_dataset_builder()
            .with_data(calibration_json.as_slice())
            .create("metadata_json_utf8")?;
        for fit in m["fits"].as_array().unwrap() {
            let id = fit["id"].as_str().unwrap();
            let record = records.create_group(id)?;
            let meta = json!({"id":id,"source_id":sid,"source_key":id,"station":"Millstone Hill","antenna":q.antenna,"channel":format!("{}-l",q.antenna.to_lowercase()),"mode":"lpi-30","kind":"ion_fit","start_unix_s":fit["start_unix_s"],"end_unix_s":fit["end_unix_s"],"radar_frequency_hz":m["radar_frequency_hz"],"density_semantics":"plasma_line_calibrated_density","density_units":"m^-3","calibration":m["calibration"]});
            record
                .new_attr_builder()
                .with_data(meta.to_string().as_str())
                .create("metadata")?;
            let ranges: Vec<_> = m["range_m"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            record
                .new_dataset_builder()
                .with_data(&ranges)
                .create("range_m")?
                .new_attr_builder()
                .with_data("m")
                .create("units")?;
            let f = record.create_group("fit")?;
            let parameters: Vec<_> = fit["parameters"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|v| v.as_f64().unwrap_or(f64::NAN))
                })
                .collect();
            let errors: Vec<f64> = (0..ranges.len())
                .flat_map(|ri| {
                    (0..4).map(move |p| fit["standard_error"][ri][p].as_f64().unwrap_or(f64::NAN))
                })
                .collect();
            for (name, values) in [("parameters", parameters), ("standard_error", errors)] {
                f.new_dataset::<f64>()
                    .shape([ranges.len(), 4])
                    .create(name)?
                    .write_raw(&values)?;
                f.dataset(name)?
                    .new_attr_builder()
                    .with_data("m^-3,K,K,m s^-1")
                    .create("units")?;
            }
        }
        h.flush()?;
        drop(h);
        atomic_json(&out.join("manifest.json"), &m)?;
        std::fs::rename(out, &dest)?;
    }
    let configuration = s.public.join("catalog/millstone-config.json");
    let base = if configuration.is_file() {
        configuration.clone()
    } else {
        s.public.join("millstone-config.json")
    };
    let mut config: Value = serde_json::from_slice(&std::fs::read(base)?)?;
    for d in config["datasets"].as_array_mut().unwrap() {
        if d["date"] == q.date && d["antenna"] == q.antenna {
            d["data_root"] = format!("data/{name}/").into();
        }
    }
    atomic_json(&configuration, &config)?;
    Ok(
        json!({"gain":gain,"revision":q.revision,"independent_fit_count":independent.len(),"gain_standard_error":se,"density_uncertainty":"Conditional source-amplitude and independent-fit gain error; Ti/Te covariance unavailable","hdf5_url":format!("data/{name}/calibration.h5"),"manifest_url":format!("data/{name}/manifest.json"),"message":"Density recalculated. Source temperatures, velocity and ACF fits retained."}),
    )
}
async fn recalibrate(
    State(s): State<Store>,
    headers: axum::http::HeaderMap,
    Json(q): Json<CalibrationRequest>,
) -> Result<Json<Value>, Error> {
    let expected = format!("Bearer {}", s.token);
    if s.token.len() < 32
        || headers
            .get(header::AUTHORIZATION)
            .and_then(|h| h.to_str().ok())
            != Some(expected.as_str())
    {
        return Err((
            StatusCode::UNAUTHORIZED,
            "Enter the processing key under Processing access.".into(),
        ));
    }
    let v = tokio::task::spawn_blocking(move || calibrate_density(&s, q))
        .await
        .map_err(fail)?
        .map_err(fail)?;
    Ok(Json(v))
}
#[cfg(test)]
mod calibration_tests {
    use super::*;
    #[test]
    fn finite_debye_density_inverts_power() {
        let k = 4. * std::f64::consts::PI * 440.2e6 / 299792458.;
        for ne in [1e9, 1e10, 1e11, 1e12] {
            for te in [500., 2000., 5000.] {
                let power = ion_proxy(ne, te, 1000., k);
                assert!((density(power, te, 1000., k) / ne - 1.).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn nonphysical_density_is_missing() {
        assert!(density(-1., 2000., 1000., 1.).is_nan());
        assert!(density(1e11, 0., 1000., 1.).is_nan());
    }
}
