//! Coherent complex-voltage cancellation with held-out validation. This is
//! amplitude-domain subtraction, not pulse blanking or post-ACF nuisance fitting.
use crate::{signal::percentile, C};
use anyhow::{ensure, Result};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Detection {
    pub delay_samples: f64,
    pub doppler_hz: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct Settings {
    pub interpolation_half: usize,
    pub validation_block: usize,
    pub max_residual_rms_ratio: f64,
    pub max_edge_rms_ratio: f64,
    pub max_coherent_sigma: f64,
    pub adc_limit: f64,
    pub edge_guard: usize,
    pub refine_delay_samples: f64,
    pub refine_doppler_hz: f64,
    pub min_amplitude_snr_db: f64,
    pub shrink_weak_amplitudes: bool,
    /// Bartlett HAC bandwidth for conditional amplitude uncertainty. Zero keeps
    /// the legacy white-background assumption; positive values are experimental.
    #[serde(skip_serializing_if = "zero_lag")]
    pub correlation_lag_s: f64,
}
fn zero_lag(value: &f64) -> bool {
    *value == 0.
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            interpolation_half: 16,
            validation_block: 16,
            max_residual_rms_ratio: 1.5,
            max_edge_rms_ratio: 2.,
            max_coherent_sigma: 5.,
            adc_limit: 32767.,
            edge_guard: 8,
            refine_delay_samples: 4.,
            refine_doppler_hz: 100.,
            min_amplitude_snr_db: 6.,
            shrink_weak_amplitudes: true,
            correlation_lag_s: 0.,
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Diagnostic {
    pub accepted: bool,
    pub applied: bool,
    pub confirmed: bool,
    pub amplitude_snr_db: Vec<f64>,
    pub amplitude_shrinkage: Vec<f64>,
    #[serde(default)]
    pub amplitude_white_variance: Vec<f64>,
    #[serde(default)]
    pub amplitude_covariance_variance: Vec<f64>,
    pub reason: String,
    pub input_rms_over_noise: f64,
    pub residual_rms_over_noise: f64,
    pub heldout_rms_over_noise: f64,
    pub edge_rms_over_noise: f64,
    pub heldout_coherent_sigma: f64,
    pub input_coherent_sigma: f64,
    pub clip_fraction: f64,
    pub fitted_complex_parameters: usize,
    pub fitted_samples: usize,
    pub projection_fraction: f64,
    pub amplitudes: Vec<[f64; 2]>,
    pub detections: Vec<Detection>,
    pub noise_rms: f64,
}
pub struct Cancellation {
    pub cleaned: Vec<C>,
    pub usable: Vec<bool>,
    pub diagnostic: Diagnostic,
}
fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-12 {
        1.
    } else {
        (PI * x).sin() / (PI * x)
    }
}
/// Noncircular fractional delay with a compact windowed-sinc interpolation of
/// the *measured* complete pulse, preserving all bit edges and receiver ringing.
fn delayed_template(tx: &[C], delay: f64, n: usize, half: usize) -> Vec<C> {
    // Zero samples cannot contribute to interpolation. Preserve the original
    // summation order for every contributing term, including pulse ringing.
    let first = tx
        .iter()
        .position(|z| z.re != 0. || z.im != 0.)
        .map(|i| i as isize)
        .unwrap_or(1);
    let last = tx
        .iter()
        .rposition(|z| z.re != 0. || z.im != 0.)
        .map(|i| i as isize)
        .unwrap_or(0);
    (0..n)
        .map(|i| {
            let t = i as f64 - delay;
            let centre = t.floor() as isize;
            let mut z = C::default();
            for j in (centre - half as isize + 1).max(first)..=(centre + half as isize).min(last) {
                let u = t - j as f64;
                if u.abs() < half as f64 {
                    z += tx[j as usize] * (sinc(u) * sinc(u / half as f64));
                }
            }
            z
        })
        .collect()
}
fn doppler_phase(doppler: f64, sr: f64, n: usize) -> Vec<C> {
    (0..n)
        .map(|i| C::from_polar(1., 2. * PI * doppler * i as f64 / sr))
        .collect()
}
pub fn template(tx: &[C], delay: f64, doppler: f64, sr: f64, n: usize, half: usize) -> Vec<C> {
    delayed_template(tx, delay, n, half)
        .into_iter()
        .enumerate()
        .map(|(i, z)| z * C::from_polar(1., 2. * PI * doppler * i as f64 / sr))
        .collect()
}
fn amplitude(y: &[C], t: &[C], mask: &[bool]) -> Option<C> {
    let mut xy = C::default();
    let mut xx = 0.;
    for ((y, t), use_it) in y.iter().zip(t).zip(mask) {
        if *use_it {
            xy += t.conj() * y;
            xx += t.norm_sqr()
        }
    }
    (xx > 1e-20).then(|| xy / xx)
}
fn rms(y: &[C], mask: &[bool]) -> f64 {
    let n = mask.iter().filter(|v| **v).count();
    if n == 0 {
        return f64::NAN;
    }
    (y.iter()
        .zip(mask)
        .filter(|(_, m)| **m)
        .map(|(z, _)| z.norm_sqr())
        .sum::<f64>()
        / n as f64)
        .sqrt()
}
fn reject_pulse(y: &[C], detections: Vec<Detection>, reason: &str, adc_limit: f64) -> Cancellation {
    Cancellation {
        cleaned: y.to_vec(),
        usable: vec![false; y.len()],
        diagnostic: Diagnostic {
            accepted: false,
            applied: false,
            confirmed: false,
            amplitude_snr_db: vec![],
            amplitude_shrinkage: vec![],
            amplitude_white_variance: vec![],
            amplitude_covariance_variance: vec![],
            reason: reason.into(),
            input_rms_over_noise: f64::NAN,
            residual_rms_over_noise: f64::NAN,
            heldout_rms_over_noise: f64::NAN,
            edge_rms_over_noise: f64::NAN,
            heldout_coherent_sigma: f64::NAN,
            input_coherent_sigma: f64::NAN,
            clip_fraction: y
                .iter()
                .filter(|z| z.re.abs() >= adc_limit || z.im.abs() >= adc_limit)
                .count() as f64
                / y.len() as f64,
            fitted_complex_parameters: 0,
            fitted_samples: 0,
            projection_fraction: 0.,
            amplitudes: vec![],
            detections,
            noise_rms: f64::NAN,
        },
    }
}
pub fn cancel(
    y: &[C],
    tx: &[C],
    detections: &[Detection],
    sr: f64,
    noise_window: [usize; 2],
    settings: &Settings,
) -> Result<Cancellation> {
    ensure!(
        !y.is_empty()
            && y.len() == tx.len()
            && y.iter()
                .chain(tx)
                .all(|z| z.re.is_finite() && z.im.is_finite()),
        "invalid voltage/template arrays"
    );
    ensure!(
        sr.is_finite()
            && sr > 0.
            && noise_window[1] <= y.len()
            && noise_window[0] < noise_window[1]
            && settings.interpolation_half > 0
            && settings.validation_block > 0
            && settings.adc_limit > 0.
            && settings.max_residual_rms_ratio > 0.
            && settings.max_edge_rms_ratio > 0.
            && settings.max_coherent_sigma > 0.
            && settings.correlation_lag_s.is_finite()
            && settings.correlation_lag_s >= 0.,
        "invalid cancellation settings"
    );
    ensure!(
        detections.iter().all(|d| d.delay_samples.is_finite()
            && d.doppler_hz.is_finite()
            && d.delay_samples >= 0.
            && d.doppler_hz.abs() < sr / 2.),
        "invalid detection"
    );
    let bg = &y[noise_window[0]..noise_window[1]];
    ensure!(bg.len() >= 32, "at least 32 background samples required");
    let dc = C::new(
        percentile(&bg.iter().map(|z| z.re).collect::<Vec<_>>(), 0.5),
        percentile(&bg.iter().map(|z| z.im).collect::<Vec<_>>(), 0.5),
    );
    // Median complex Gaussian power / ln 2 is robust against occasional impulsive
    // contamination; its finite-sample uncertainty is retained as a QA limitation.
    let noise2 = percentile(
        &bg.iter().map(|z| (*z - dc).norm_sqr()).collect::<Vec<_>>(),
        0.5,
    ) / 2f64.ln();
    if !noise2.is_finite() || noise2 <= 0. {
        // A flat/dead recorded background is a data-quality failure, not an
        // invented noise floor or a reason to abort unrelated source windows.
        return Ok(reject_pulse(
            y,
            detections.to_vec(),
            "invalid_background_noise",
            settings.adc_limit,
        ));
    }
    let noise = noise2.sqrt();
    let mut refined = Vec::new();
    let mut templates = Vec::new();
    for d in detections {
        let mut centre = d.clone();
        let mut selected = (d.clone(), Vec::new());
        let levels = if settings.refine_delay_samples == 0. && settings.refine_doppler_hz == 0. {
            1
        } else {
            6
        };
        for level in 0..levels {
            let mut best = f64::INFINITY;
            // Delay and Doppler are separable. Reuse the five interpolated
            // pulses and five phase ramps instead of recomputing each 25 times.
            let delays: Vec<_> = (-2..=2)
                .map(|di| {
                    centre.delay_samples
                        + di as f64 * settings.refine_delay_samples / (2. * 4f64.powi(level))
                })
                .collect();
            let frequencies: Vec<_> = (-2..=2)
                .map(|fi| {
                    centre.doppler_hz
                        + fi as f64 * settings.refine_doppler_hz / (2. * 4f64.powi(level))
                })
                .collect();
            let delayed: Vec<_> = delays
                .iter()
                .map(|&dd| delayed_template(tx, dd, y.len(), settings.interpolation_half))
                .collect();
            let phases: Vec<_> = frequencies
                .iter()
                .map(|&ff| doppler_phase(ff, sr, y.len()))
                .collect();
            for di in -2..=2 {
                for fi in -2..=2 {
                    let dd = centre.delay_samples
                        + di as f64 * settings.refine_delay_samples / (2. * 4f64.powi(level));
                    let ff = centre.doppler_hz
                        + fi as f64 * settings.refine_doppler_hz / (2. * 4f64.powi(level));
                    let t: Vec<_> = delayed[(di + 2) as usize]
                        .iter()
                        .zip(&phases[(fi + 2) as usize])
                        .map(|(z, p)| *z * *p)
                        .collect();
                    let peak = t.iter().map(|z| z.norm_sqr()).fold(0., f64::max);
                    let train: Vec<_> = t
                        .iter()
                        .enumerate()
                        .map(|(i, z)| {
                            z.norm_sqr() > peak * 1e-6
                                && (i / settings.validation_block) % 2 == 0
                                && !(noise_window[0]..noise_window[1]).contains(&i)
                        })
                        .collect();
                    if let Some(a) = amplitude(y, &t, &train) {
                        let cost = y
                            .iter()
                            .zip(&t)
                            .zip(&train)
                            .filter(|(_, m)| **m)
                            .map(|((y, t), _)| (*y - a * t).norm_sqr())
                            .sum::<f64>()
                            / train.iter().filter(|x| **x).count().max(1) as f64;
                        if cost < best {
                            best = cost;
                            selected = (
                                Detection {
                                    delay_samples: dd,
                                    doppler_hz: ff,
                                },
                                t,
                            )
                        }
                    }
                }
            }
            centre = selected.0.clone();
        }
        ensure!(
            !selected.1.is_empty(),
            "detection has no measurable template support"
        );
        refined.push(selected.0);
        templates.push(selected.1);
    }
    let mut support = vec![false; y.len()];
    let mut edges = vec![false; y.len()];
    for t in &templates {
        let peak = t.iter().map(|z| z.norm()).fold(0., f64::max);
        for i in 0..t.len() {
            if t[i].norm() > peak * 1e-3 {
                support[i] = true
            }
            if i > 0 && (t[i] - t[i - 1]).norm() > peak * 0.15 {
                for bit in edges
                    .iter_mut()
                    .take((i + settings.edge_guard + 1).min(t.len()))
                    .skip(i.saturating_sub(settings.edge_guard))
                {
                    *bit = true
                }
            }
        }
    }
    if support
        .iter()
        .enumerate()
        .any(|(i, s)| *s && (noise_window[0]..noise_window[1]).contains(&i))
    {
        // No independent background calibration is available for this pulse.
        // Keep original voltages and reject it, including receiver-gain use.
        // Undefined QA ratios are serialized as null, never fabricated values.
        return Ok(reject_pulse(
            y,
            refined,
            "background_gate_contaminated",
            settings.adc_limit,
        ));
    }
    let train: Vec<_> = support
        .iter()
        .enumerate()
        .map(|(i, s)| *s && (i / settings.validation_block) % 2 == 0)
        .collect();
    let hold: Vec<_> = support
        .iter()
        .enumerate()
        .map(|(i, s)| *s && !train[i])
        .collect();
    let ntrain = train.iter().filter(|v| **v).count();
    let nhold = hold.iter().filter(|v| **v).count();
    let ncol = templates.len();
    if ncol == 0 {
        return Ok(Cancellation {
            cleaned: y.to_vec(),
            usable: vec![true; y.len()],
            diagnostic: Diagnostic {
                accepted: true,
                applied: false,
                confirmed: false,
                amplitude_snr_db: vec![],
                amplitude_shrinkage: vec![],
                amplitude_white_variance: vec![],
                amplitude_covariance_variance: vec![],
                reason: "no_detections".into(),
                input_rms_over_noise: 1.,
                residual_rms_over_noise: 1.,
                heldout_rms_over_noise: 1.,
                edge_rms_over_noise: 1.,
                heldout_coherent_sigma: 0.,
                input_coherent_sigma: 0.,
                clip_fraction: 0.,
                fitted_complex_parameters: 0,
                fitted_samples: 0,
                projection_fraction: 0.,
                amplitudes: vec![],
                detections: vec![],
                noise_rms: noise,
            },
        });
    }
    ensure!(
        ntrain > ncol * 8 && nhold > 16,
        "insufficient fit/validation support"
    );
    let scales: Vec<_> = templates
        .iter()
        .map(|t| {
            t.iter()
                .zip(&train)
                .filter(|(_, m)| **m)
                .map(|(z, _)| z.norm_sqr())
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    let a = DMatrix::from_fn(ntrain, ncol, |_, _| C::default());
    let mut a = a;
    let mut b = DVector::zeros(ntrain);
    let mut row = 0;
    for i in 0..y.len() {
        if train[i] {
            b[row] = y[i];
            for j in 0..ncol {
                a[(row, j)] = templates[j][i] / scales[j]
            }
            row += 1
        }
    }
    let svd = a.svd(true, true);
    let singular = svd.singular_values.clone();
    ensure!(
        singular.iter().all(|s| *s > 1e-8),
        "overlapping templates are rank deficient"
    );
    let coefficients = svd.solve(&b, 1e-10).map_err(|e| anyhow::anyhow!(e))?;
    let mut amplitudes: Vec<_> = coefficients
        .iter()
        .zip(&scales)
        .map(|(a, s)| *a / s)
        .collect();
    // Conditional amplitude variance from the joint scaled design. The
    // delay/Doppler search and template uncertainty are not included, so this
    // quantity controls conservative subtraction, not a claimed full posterior.
    let v_t = svd.v_t.as_ref().unwrap();
    let white_variance: Vec<_> = (0..ncol)
        .map(|j| {
            noise2
                * (0..ncol)
                    .map(|k| v_t[(k, j)].norm_sqr() / singular[k].powi(2))
                    .sum::<f64>()
                / scales[j].powi(2)
        })
        .collect();
    let mut variance = white_variance.clone();
    if settings.correlation_lag_s > 0. {
        // Sandwich/HAC covariance of the actual complex least-squares estimator.
        // A+ = V S^-1 U^H, including joint-template conditioning and missing
        // training samples. Correlation is measured after target fitting, before
        // shrinkage. Do not assume alternating training/hold blocks are independent.
        let u = svd.u.as_ref().unwrap();
        let locations: Vec<_> = train
            .iter()
            .enumerate()
            .filter_map(|(i, &yes)| yes.then_some(i))
            .collect();
        let first = *locations.first().unwrap();
        let last = *locations.last().unwrap();
        let residual: Vec<_> = locations
            .iter()
            .map(|&i| {
                y[i] - templates
                    .iter()
                    .zip(&amplitudes)
                    .map(|(t, a)| t[i] * *a)
                    .sum::<C>()
            })
            .collect();
        let bandwidth = ((settings.correlation_lag_s * sr).ceil() as usize).min(last - first);
        for j in 0..ncol {
            let mut influence = vec![C::default(); last - first + 1];
            for (row, (&i, &e)) in locations.iter().zip(&residual).enumerate() {
                let weight = (0..ncol)
                    .map(|k| v_t[(k, j)].conj() * u[(row, k)].conj() / singular[k])
                    .sum::<C>()
                    / scales[j];
                influence[i - first] = weight * e;
            }
            let mut estimate = influence.iter().map(|z| z.norm_sqr()).sum::<f64>();
            for lag in 1..=bandwidth {
                let covariance = influence[..influence.len() - lag]
                    .iter()
                    .zip(&influence[lag..])
                    .map(|(a, b)| (*a * b.conj()).re)
                    .sum::<f64>();
                estimate += 2. * (1. - lag as f64 / (bandwidth + 1) as f64) * covariance;
            }
            estimate *= ntrain as f64 / (ntrain - ncol) as f64;
            variance[j] = variance[j].max(estimate);
        }
    }
    let amplitude_snr_db: Vec<_> = amplitudes
        .iter()
        .zip(&variance)
        .map(|(a, v)| 10. * ((a.norm_sqr() / v - 1.).max(1e-12)).log10())
        .collect();
    let shrinkage: Vec<_> = amplitudes
        .iter()
        .zip(&variance)
        .zip(&amplitude_snr_db)
        .map(|((a, v), snr)| {
            if *snr < settings.min_amplitude_snr_db {
                0.
            } else if settings.shrink_weak_amplitudes && *snr < 14. {
                (1. - v / a.norm_sqr()).clamp(0., 1.)
            } else {
                1.
            }
        })
        .collect();
    for (a, shrink) in amplitudes.iter_mut().zip(&shrinkage) {
        *a *= shrink
    }
    let applied = shrinkage.iter().any(|v| *v > 0.);
    let mut cleaned = y.to_vec();
    for (i, z) in cleaned.iter_mut().enumerate() {
        for (t, a) in templates.iter().zip(&amplitudes) {
            *z -= a * t[i]
        }
    }
    let n = support.iter().filter(|x| **x).count();
    let clipped = y
        .iter()
        .zip(&support)
        .filter(|(z, s)| {
            **s && (z.re.abs() >= settings.adc_limit - 1. || z.im.abs() >= settings.adc_limit - 1.)
        })
        .count();
    let clip_fraction = clipped as f64 / n.max(1) as f64;
    let coherent = templates
        .iter()
        .map(|t| {
            amplitude(&cleaned, t, &hold)
                .map(|a| {
                    a.norm()
                        * t.iter()
                            .zip(&hold)
                            .filter(|(_, m)| **m)
                            .map(|(z, _)| z.norm_sqr())
                            .sum::<f64>()
                            .sqrt()
                        / noise
                })
                .unwrap_or(f64::INFINITY)
        })
        .fold(0., f64::max);
    let input_coherent = templates
        .iter()
        .map(|t| {
            amplitude(y, t, &hold)
                .map(|a| {
                    a.norm()
                        * t.iter()
                            .zip(&hold)
                            .filter(|(_, m)| **m)
                            .map(|(z, _)| z.norm_sqr())
                            .sum::<f64>()
                            .sqrt()
                        / noise
                })
                .unwrap_or(0.)
        })
        .fold(0., f64::max);
    let residual = rms(&cleaned, &support) / noise;
    let heldout = rms(&cleaned, &hold) / noise;
    let edge = rms(&cleaned, &edges) / noise;
    let reason = if clipped > 0 {
        "clipped_adc"
    } else if !applied {
        "skipped_low_snr"
    } else if coherent > settings.max_coherent_sigma {
        "coherent_residual"
    } else if !heldout.is_finite() || heldout > settings.max_residual_rms_ratio {
        "heldout_excess"
    } else if edge.is_finite() && edge > settings.max_edge_rms_ratio {
        "edge_residual"
    } else if input_coherent < 5. {
        "weak_applied"
    } else {
        "accepted"
    };
    let accepted = ["accepted", "weak_applied", "skipped_low_snr"].contains(&reason);
    let confirmed = reason == "accepted";
    let mut usable = vec![true; y.len()];
    if !accepted {
        for i in 0..y.len() {
            if support[i] {
                for v in usable
                    .iter_mut()
                    .take((i + settings.edge_guard + 1).min(y.len()))
                    .skip(i.saturating_sub(settings.edge_guard))
                {
                    *v = false
                }
            }
        }
    }
    Ok(Cancellation {
        cleaned,
        usable,
        diagnostic: Diagnostic {
            accepted,
            applied,
            confirmed,
            amplitude_snr_db,
            amplitude_shrinkage: shrinkage,
            amplitude_white_variance: white_variance,
            amplitude_covariance_variance: variance,
            reason: reason.into(),
            input_rms_over_noise: rms(y, &support) / noise,
            residual_rms_over_noise: residual,
            heldout_rms_over_noise: heldout,
            edge_rms_over_noise: edge,
            heldout_coherent_sigma: coherent,
            input_coherent_sigma: input_coherent,
            clip_fraction,
            fitted_complex_parameters: ncol,
            fitted_samples: ntrain,
            projection_fraction: ncol as f64 / ntrain as f64,
            amplitudes: amplitudes.iter().map(|z| [z.re, z.im]).collect(),
            detections: refined,
            noise_rms: noise,
        },
    })
}
