//! Experimental independent-pulse echo prediction. Never uses the queried
//! pulse in association, phase/geometry fitting or amplitude estimation.
//! This module is not enabled in production cancellation.
use crate::C;
use anyhow::{ensure, Result};
use nalgebra::{Matrix3, Vector3};
use rustfft::FftPlanner;
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, f64::consts::TAU};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Event {
    pub sample: u64,
    pub delay_samples: f64,
    pub doppler_hz: f64,
    pub amplitude: C,
    pub snr_db: f64,
}
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(default)]
pub struct Settings {
    pub radius_s: f64,
    pub min_pulses: usize,
    pub max_pulses: usize,
    pub sample_rate_hz: f64,
    pub radar_frequency_hz: f64,
    pub max_search_cells: usize,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            radius_s: 0.35,
            min_pulses: 16,
            max_pulses: 48,
            sample_rate_hz: 1e6,
            radar_frequency_hz: 440.2e6,
            max_search_cells: 2_000_000,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Forecast {
    pub supported: bool,
    pub reason: String,
    pub training_samples: Vec<u64>,
    pub training_count: usize,
    pub pulse_clock_samples: u64,
    pub delay_samples: f64,
    pub doppler_hz: f64,
    pub delay_prediction_std_samples: f64,
    pub doppler_prediction_std_hz: f64,
    pub delay_rate_samples_s: f64,
    pub doppler_slope_hz_s: f64,
    pub doppler_from_range_rate_hz: f64,
    pub amplitude: C,
    pub centre_amplitude: C,
    pub phase_frequency_alias_hz: f64,
    pub phase_chirp_hz_s: f64,
    pub explained_fraction: f64,
    pub training_amplitude_residual_rms: f64,
    pub prediction_variance: f64,
    pub search_cells: usize,
    pub phase_fft_length: usize,
}
fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}
fn dt(sample: u64, held: u64, sr: f64) -> f64 {
    (sample as i128 - held as i128) as f64 / sr
}
fn valid(e: &Event) -> bool {
    e.delay_samples.is_finite()
        && e.doppler_hz.is_finite()
        && e.amplitude.re.is_finite()
        && e.amplitude.im.is_finite()
        && e.snr_db.is_finite()
}

/// Loose local association. Both range-rate signs are retained here; physical
/// consistency is reported separately rather than guessed from ion velocity.
pub fn groups(events: &[Event], held: u64, s: &Settings) -> Result<Vec<Vec<Event>>> {
    ensure!(
        s.sample_rate_hz == 1e6
            && s.radar_frequency_hz == 440.2e6
            && s.radius_s.is_finite()
            && s.radius_s > 0.,
        "unsupported_track_layout"
    );
    let rows: Vec<_> = events
        .iter()
        .filter(|e| {
            e.sample != held
                && dt(e.sample, held, s.sample_rate_hz).abs() <= s.radius_s
                && e.snr_db >= 14.
        })
        .cloned()
        .collect();
    ensure!(rows.iter().all(valid), "invalid_track_metadata");
    let mut parent: Vec<_> = (0..rows.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for i in 0..rows.len() {
        for j in 0..i {
            let a = &rows[i];
            let b = &rows[j];
            let time = dt(a.sample, b.sample, s.sample_rate_hz).abs();
            let drift = a.doppler_hz.abs().max(b.doppler_hz.abs()) / s.radar_frequency_hz
                * s.sample_rate_hz
                * time;
            if (a.delay_samples - b.delay_samples).abs() <= 16. + drift
                && (a.doppler_hz - b.doppler_hz).abs() <= 2000.
            {
                let ri = root(&mut parent, i);
                let rj = root(&mut parent, j);
                parent[ri] = rj;
            }
        }
    }
    let mut result: Vec<(usize, Vec<Event>)> = Vec::new();
    for (i, event) in rows.into_iter().enumerate() {
        let r = root(&mut parent, i);
        if let Some((_, group)) = result.iter_mut().find(|(k, _)| *k == r) {
            group.push(event)
        } else {
            result.push((r, vec![event]));
        }
    }
    Ok(result.into_iter().map(|(_, rows)| rows).collect())
}

pub struct Forecaster {
    planner: FftPlanner<f64>,
}
impl Default for Forecaster {
    fn default() -> Self {
        Self {
            planner: FftPlanner::new(),
        }
    }
}
struct LinearPhase<'a> {
    u: &'a [f64],
    b: &'a [C],
    inverse: [[f64; 2]; 2],
    energy: f64,
}
impl LinearPhase<'_> {
    fn solve(&self, f: f64, c: f64) -> ([C; 2], f64) {
        let mut r = [C::default(); 2];
        for (&u, &b) in self.u.iter().zip(self.b) {
            let z = b * C::from_polar(1., -TAU * (f * u + 0.5 * c * u * u));
            r[0] += z;
            r[1] += z * u;
        }
        let coefficients = [
            r[0] * self.inverse[0][0] + r[1] * self.inverse[0][1],
            r[0] * self.inverse[1][0] + r[1] * self.inverse[1][1],
        ];
        let explained = (r[0].conj() * coefficients[0] + r[1].conj() * coefficients[1]).re;
        (coefficients, self.energy - explained)
    }
    fn refine(&self, f: f64, c: f64) -> (f64, f64) {
        // Two-dimensional Nelder-Mead refinement of the variable-projection
        // objective. A shrinking coordinate grid stalls in its narrow valleys.
        let initial = [f, c];
        let initial_cost = self.solve(f, c).1;
        let evaluate = |p: [f64; 2]| (p, self.solve(p[0], p[1]).1);
        let mut simplex = [
            evaluate(initial),
            evaluate([if f != 0. { 1.05 * f } else { 0.00025 }, c]),
            evaluate([f, if c != 0. { 1.05 * c } else { 0.00025 }]),
        ];
        for _ in 1..150 {
            simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
            let mut spread: f64 = 0.;
            for i in 1..3 {
                for j in 0..2 {
                    spread = spread.max((simplex[i].0[j] - simplex[0].0[j]).abs());
                }
            }
            let value_spread = (simplex[2].1 - simplex[0].1)
                .abs()
                .max((simplex[1].1 - simplex[0].1).abs());
            if spread <= 1e-7 && value_spread <= 1e-14 {
                return if simplex[0].1 <= initial_cost {
                    (simplex[0].0[0], simplex[0].0[1])
                } else {
                    (f, c)
                };
            }
            let middle = [
                (simplex[0].0[0] + simplex[1].0[0]) / 2.,
                (simplex[0].0[1] + simplex[1].0[1]) / 2.,
            ];
            let reflected = evaluate([
                2. * middle[0] - simplex[2].0[0],
                2. * middle[1] - simplex[2].0[1],
            ]);
            let mut shrink = false;
            if reflected.1 < simplex[0].1 {
                let expanded = evaluate([
                    3. * middle[0] - 2. * simplex[2].0[0],
                    3. * middle[1] - 2. * simplex[2].0[1],
                ]);
                simplex[2] = if expanded.1 < reflected.1 {
                    expanded
                } else {
                    reflected
                };
            } else if reflected.1 < simplex[1].1 {
                simplex[2] = reflected;
            } else if reflected.1 < simplex[2].1 {
                let contracted = evaluate([
                    1.5 * middle[0] - 0.5 * simplex[2].0[0],
                    1.5 * middle[1] - 0.5 * simplex[2].0[1],
                ]);
                if contracted.1 <= reflected.1 {
                    simplex[2] = contracted
                } else {
                    shrink = true;
                }
            } else {
                let contracted = evaluate([
                    0.5 * middle[0] + 0.5 * simplex[2].0[0],
                    0.5 * middle[1] + 0.5 * simplex[2].0[1],
                ]);
                if contracted.1 < simplex[2].1 {
                    simplex[2] = contracted
                } else {
                    shrink = true;
                }
            }
            if shrink {
                for i in 1..3 {
                    simplex[i] = evaluate([
                        simplex[0].0[0] + 0.5 * (simplex[i].0[0] - simplex[0].0[0]),
                        simplex[0].0[1] + 0.5 * (simplex[i].0[1] - simplex[0].0[1]),
                    ]);
                }
            }
        }
        (f, c)
    }
}
impl Forecaster {
    pub fn forecast(
        &mut self,
        events: &[Event],
        held: u64,
        tx_centre: f64,
        s: &Settings,
    ) -> Result<Forecast> {
        ensure!(
            s.sample_rate_hz == 1e6 && s.radar_frequency_hz == 440.2e6,
            "unsupported_track_layout"
        );
        ensure!(
            s.min_pulses >= 16
                && s.max_pulses >= s.min_pulses
                && s.max_pulses <= 256
                && tx_centre.is_finite(),
            "invalid_track_settings"
        );
        let mut rows: Vec<_> = events
            .iter()
            .filter(|e| e.sample != held)
            .cloned()
            .collect();
        ensure!(rows.iter().all(valid), "invalid_track_metadata");
        let unique: BTreeSet<_> = rows.iter().map(|e| e.sample).collect();
        ensure!(
            unique.len() == rows.len(),
            "multiple_components_in_training_pulse"
        );
        rows.sort_by_key(|e| (e.sample.abs_diff(held), e.sample));
        rows.truncate(s.max_pulses);
        ensure!(rows.len() >= s.min_pulses, "insufficient_other_pulses");
        rows.sort_by_key(|e| e.sample);
        let clock = rows
            .windows(2)
            .map(|w| w[1].sample - w[0].sample)
            .fold(0, gcd);
        ensure!(
            clock > 0 && held.abs_diff(rows[0].sample) % clock == 0,
            "held_epoch_not_on_training_clock"
        );
        ensure!(clock >= 1000, "unsupported_irregular_clock");
        let u: Vec<_> = rows
            .iter()
            .map(|e| dt(e.sample, held, s.sample_rate_hz))
            .collect();
        let span = u[u.len() - 1] - u[0];
        ensure!(span > 0., "zero_training_span");
        let b: Vec<_> = rows
            .iter()
            .map(|e| {
                e.amplitude
                    * C::from_polar(
                        1.,
                        TAU * e.doppler_hz * (tx_centre + e.delay_samples) / s.sample_rate_hz,
                    )
            })
            .collect();
        let n = rows.len() as f64;
        let su = u.iter().sum::<f64>();
        let suu = u.iter().map(|u| u * u).sum::<f64>();
        let det = n * suu - su * su;
        ensure!(det > 0., "singular_phase_design");
        let inverse = [[suu / det, -su / det], [-su / det, n / det]];
        let energy = b.iter().map(|z| z.norm_sqr()).sum::<f64>();
        ensure!(energy > 0., "zero_training_amplitude");
        let dc = 1. / (2. * span * span);
        // Same grid endpoint convention as NumPy arange(-1000,1000+dc/2,dc).
        let nchirp = ((2000. + dc * 0.5) / dc).ceil() as usize;
        let extent = ((rows.last().unwrap().sample - rows[0].sample) / clock) as usize + 1;
        let nfft = (4 * extent).next_power_of_two();
        ensure!(
            nfft.checked_mul(nchirp)
                .is_some_and(|cells| cells <= s.max_search_cells),
            "phase_search_budget_exceeded"
        );
        let fft = self.planner.plan_fft_forward(nfft);
        let mut z0 = vec![C::default(); nfft];
        let mut z1 = z0.clone();
        let mut scratch = vec![C::default(); fft.get_inplace_scratch_len()];
        let mut best = (f64::NEG_INFINITY, 0., 0.);
        for ci in 0..nchirp {
            let chirp = -1000. + ci as f64 * dc;
            z0.fill(C::default());
            z1.fill(C::default());
            for (i, event) in rows.iter().enumerate() {
                let k = ((event.sample - rows[0].sample) / clock) as usize;
                let z = b[i] * C::from_polar(1., -std::f64::consts::PI * chirp * u[i] * u[i]);
                z0[k] = z;
                z1[k] = z * u[i];
            }
            fft.process_with_scratch(&mut z0, &mut scratch);
            fft.process_with_scratch(&mut z1, &mut scratch);
            for i in 0..nfft {
                let explained = inverse[0][0] * z0[i].norm_sqr()
                    + inverse[1][1] * z1[i].norm_sqr()
                    + 2. * inverse[0][1] * (z0[i].conj() * z1[i]).re;
                if explained > best.0 {
                    let bin = if i < nfft.div_ceil(2) {
                        i as f64
                    } else {
                        i as f64 - nfft as f64
                    };
                    best = (
                        explained,
                        bin * s.sample_rate_hz / clock as f64 / nfft as f64,
                        chirp,
                    );
                }
            }
        }
        let phase = LinearPhase {
            u: &u,
            b: &b,
            inverse,
            energy,
        };
        let (f, c) = phase.refine(best.1, best.2);
        let (coefficients, residual) = phase.solve(f, c);
        let mut normal = Matrix3::<f64>::zeros();
        let mut rd = Vector3::<f64>::zeros();
        let mut rf = rd;
        for (&u, event) in u.iter().zip(&rows) {
            let x = Vector3::new(1., u, u * u);
            normal += x * x.transpose();
            rd += x * event.delay_samples;
            rf += x * event.doppler_hz;
        }
        let factor = normal
            .cholesky()
            .ok_or_else(|| anyhow::anyhow!("singular_geometry_design"))?;
        let dfit = factor.solve(&rd);
        let ffit = factor.solve(&rf);
        let covariance = factor.inverse()[(0, 0)];
        let mut dv = 0.;
        let mut fv = 0.;
        for (&u, event) in u.iter().zip(&rows) {
            let x = Vector3::new(1., u, u * u);
            dv += (event.delay_samples - x.dot(&dfit)).powi(2);
            fv += (event.doppler_hz - x.dot(&ffit)).powi(2);
        }
        let amplitude = coefficients[0]
            * C::from_polar(
                1.,
                -TAU * ffit[0] * (tx_centre + dfit[0]) / s.sample_rate_hz,
            );
        let variance = residual.max(0.) / (n - 4.).max(1.);
        Ok(Forecast {
            supported: true,
            reason: "training_forecast".into(),
            training_samples: rows.iter().map(|e| e.sample).collect(),
            training_count: rows.len(),
            pulse_clock_samples: clock,
            delay_samples: dfit[0],
            doppler_hz: ffit[0],
            delay_prediction_std_samples: (dv / (n - 3.).max(1.) * covariance).sqrt(),
            doppler_prediction_std_hz: (fv / (n - 3.).max(1.) * covariance).sqrt(),
            delay_rate_samples_s: dfit[1],
            doppler_slope_hz_s: ffit[1],
            doppler_from_range_rate_hz: -s.radar_frequency_hz * dfit[1] / s.sample_rate_hz,
            amplitude,
            centre_amplitude: coefficients[0],
            phase_frequency_alias_hz: f,
            phase_chirp_hz_s: c,
            explained_fraction: 1. - residual / energy,
            training_amplitude_residual_rms: variance.sqrt(),
            prediction_variance: variance * inverse[0][0],
            search_cells: nfft * nchirp,
            phase_fft_length: nfft,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (Vec<Event>, u64) {
        let mut sample = 1_712_587_630_000_000u64;
        let mut keys = vec![sample];
        for i in 0..40 {
            sample += if i % 2 == 0 { 8910 } else { 17820 };
            keys.push(sample);
        }
        let held = keys[20];
        let rows = keys
            .into_iter()
            .map(|sample| {
                let t = dt(sample, held, 1e6);
                let d = 3700. + 3. * t + 2. * t * t;
                let f = -1000. - 250. * t + 0.5 * t * t;
                let b = C::from_polar(0.1 + 0.02 * t, TAU * (12. * t - 0.5 * 250. * t * t));
                let a = b * C::from_polar(1., -TAU * f * (350. + d) / 1e6);
                Event {
                    sample,
                    delay_samples: d,
                    doppler_hz: f,
                    amplitude: a,
                    snr_db: 30.,
                }
            })
            .collect();
        (rows, held)
    }
    #[test]
    fn unseen_pulse_is_forecast_without_its_measurements() {
        let (mut rows, held) = fixture();
        let s = Settings::default();
        let mut p = Forecaster::default();
        let associated = groups(&rows, held, &s).unwrap();
        assert_eq!(associated.len(), 1);
        let a = p.forecast(&associated[0], held, 350., &s).unwrap();
        assert!(!a.training_samples.contains(&held));
        assert!((a.delay_samples - 3700.).abs() < 1e-7);
        assert!((a.doppler_hz + 1000.).abs() < 1e-7);
        let expected = C::from_polar(0.1, TAU * 1000. * 4050. / 1e6);
        assert!((a.amplitude - expected).norm() < 0.001);
        assert!(a.explained_fraction > 0.999);
        rows[20].amplitude = C::new(0., 1000.);
        rows[20].delay_samples = 99999.;
        rows[20].doppler_hz = 50000.;
        rows[20].snr_db = -100.;
        let b = p
            .forecast(&groups(&rows, held, &s).unwrap()[0], held, 350., &s)
            .unwrap();
        assert_eq!(
            serde_json::to_string(&a).unwrap(),
            serde_json::to_string(&b).unwrap()
        );
    }
    #[test]
    fn duplicate_components_and_unseen_clock_phase_are_explicit() {
        let (mut rows, held) = fixture();
        let mut p = Forecaster::default();
        let s = Settings::default();
        rows.push(rows[0].clone());
        assert!(p
            .forecast(&rows, held, 350., &s)
            .unwrap_err()
            .to_string()
            .contains("multiple_components"));
        rows.pop();
        for (i, e) in rows.iter_mut().enumerate() {
            e.sample = 100000 + i as u64 * 20000;
        }
        assert!(p
            .forecast(&rows, 110000, 350., &s)
            .unwrap_err()
            .to_string()
            .contains("held_epoch_not_on_training_clock"));
    }
}
