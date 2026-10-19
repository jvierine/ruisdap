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
impl Record {
    fn value(&self, column: usize) -> Option<f64> {
        match column {
            0 => self.system_temperature_k,
            1 => Some(self.transmitter_w / 1e6),
            2 => Some(self.azimuth_deg),
            _ => Some(self.elevation_deg),
        }
        .filter(|v| v.is_finite())
    }
}
pub(super) fn render(v: &Viewer) {
    let Some(data) = &v.instruments else {
        text(
            "fit-instrument-summary",
            "Instrument diagnostics unavailable.",
        );
        return;
    };
    let selected = &v.data.profiles[v.time];
    let mut summary = Vec::new();
    for (column, title) in [
        "System noise temperature / K",
        "Transmit power / MW",
        "Azimuth / degrees",
        "Elevation / degrees",
    ]
    .iter()
    .enumerate()
    {
        let mut values: Vec<_> = data
            .records
            .iter()
            .filter(|r| r.end_unix > v.window_start && r.start_unix < v.window_end)
            .filter_map(|r| r.value(column).filter(|a| column != 0 || *a > 0.))
            .collect();
        let transform = |a: f64| if column == 0 { a.log10() } else { a };
        let low = values
            .iter()
            .copied()
            .map(transform)
            .fold(f64::INFINITY, f64::min);
        let high = values
            .iter()
            .copied()
            .map(transform)
            .fold(f64::NEG_INFINITY, f64::max);
        let pad = ((high - low) * 0.1).max(high.abs() * 0.002).max(0.01);
        let (lo, hi) = if values.is_empty() {
            (0., 1.)
        } else if column == 0 {
            values.sort_by(f64::total_cmp);
            let n = values.len();
            (
                low - pad,
                (2. * (values[(n - 1) / 2] + values[n / 2]) / 2.).log10(),
            )
        } else {
            (low - pad, high + pad)
        };
        let mut segments = Vec::new();
        let mut points = Vec::new();
        let mut above = 0;
        for r in &data.records {
            if r.end_unix <= v.window_start || r.start_unix >= v.window_end {
                continue;
            }
            if let Some(value) = r.value(column).filter(|a| column != 0 || *a > 0.) {
                let y = transform(value);
                if y > hi {
                    above += 1;
                    continue;
                }
                segments.extend([
                    r.start_unix.max(v.window_start),
                    y,
                    r.end_unix.min(v.window_end),
                    y,
                ]);
                points.extend([(r.start_unix + r.end_unix) / 2., y]);
            }
        }
        let spec = serde_json::json!({"title":title,"height":220,"bounds":[v.window_start,v.window_end,lo,hi],"xlabel":"Time / UTC","xformat":"utc","ylabel":title,"yformat":if column==0{"log10"}else{""},"series":[{"xy":segments,"segments":true},{"xy":points,"line":false,"radius":1.2},{"xy":[(selected.start_unix+selected.end_unix)/2.,lo,(selected.start_unix+selected.end_unix)/2.,hi],"color":"#986327"}],"caption":if above>0{format!("{above} samples above the displayed limit (2 × median). Source values retained.")}else{String::new()}});
        crate::plotting::shared_curve(&format!("fit-instrument-{column}"), &spec.to_string());
        let mut sum = 0.;
        let mut duration = 0.;
        for r in &data.records {
            let dt =
                (r.end_unix.min(selected.end_unix) - r.start_unix.max(selected.start_unix)).max(0.);
            if let Some(value) = r.value(column) {
                sum += dt * value;
                duration += dt;
            }
        }
        summary.push(if duration > 0. {
            format!(
                "{title}: {:.2} ({} s covered)",
                sum / duration,
                duration.round()
            )
        } else {
            format!("{title}: unavailable")
        });
    }
    text(
        "fit-instrument-summary",
        &format!("Selected interval · {}", summary.join(" · ")),
    );
}
