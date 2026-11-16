//! Unmagnetized Maxwellian ISR spectrum with number-conserving BGK collisions.
//! Spectral density is per Hz per electron; integrate both frequency signs.
use crate::Complex as C;
use anyhow::{ensure, Result};
use errorfunctions::{ComplexErrorFunctions, RealErrorFunctions};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
const KB: f64 = 1.380649e-23;
const ME: f64 = 9.1093837139e-31;
const AMU: f64 = 1.66053906892e-27;
const EPS0: f64 = 8.8541878188e-12;
const E: f64 = 1.602176634e-19;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IonSpecies {
    pub mass_amu: f64,
    /// Number fraction. This initial model supports singly charged ions.
    pub fraction: f64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SpectrumModel {
    pub radar_frequency_hz: f64,
    pub ions: Vec<IonSpecies>,
    pub electron_collision_s: f64,
}
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct PlasmaState {
    pub electron_density_m3: f64,
    pub ion_temperature_k: f64,
    pub electron_temperature_k: f64,
    /// Positive shifts the model spectrum toward positive baseband frequency.
    /// Physical toward/away sign must be established from the receiver convention.
    pub velocity_m_s: f64,
    /// Effective number-conserving BGK rate in s^-1, not cycles per second.
    pub ion_collision_s: f64,
}
impl SpectrumModel {
    pub fn validate(&self) -> Result<()> {
        ensure!(
            self.radar_frequency_hz.is_finite() && self.radar_frequency_hz > 0.,
            "invalid radar frequency"
        );
        ensure!(
            !self.ions.is_empty()
                && self.ions.iter().all(|i| i.mass_amu.is_finite()
                    && i.mass_amu > 0.
                    && i.fraction.is_finite()
                    && i.fraction >= 0.),
            "invalid composition"
        );
        ensure!(
            (self.ions.iter().map(|i| i.fraction).sum::<f64>() - 1.).abs() < 1e-10,
            "ion fractions must sum to one"
        );
        ensure!(
            self.electron_collision_s.is_finite() && self.electron_collision_s >= 0.,
            "invalid electron collision rate"
        );
        Ok(())
    }
    pub fn wavenumber(&self) -> f64 {
        4. * PI * self.radar_frequency_hz / 299792458.
    }
    pub fn evaluate(&self, state: PlasmaState, frequency_hz: &[f64]) -> Result<Vec<f64>> {
        self.validate()?;
        ensure!(
            state.electron_density_m3.is_finite()
                && state.electron_density_m3 > 0.
                && state.ion_temperature_k.is_finite()
                && state.ion_temperature_k > 0.
                && state.electron_temperature_k.is_finite()
                && state.electron_temperature_k > 0.
                && state.velocity_m_s.is_finite()
                && state.ion_collision_s.is_finite()
                && state.ion_collision_s >= 0.,
            "invalid plasma state"
        );
        ensure!(
            frequency_hz.iter().all(|f| f.is_finite()),
            "invalid frequency grid"
        );
        let k = self.wavenumber();
        let alpha = state.electron_density_m3 * E * E / (EPS0 * KB * k * k);
        // Thermal widths and susceptibility factors are constant over a grid.
        // In the collisionless case the Faddeeva argument is real: evaluate
        // its Gaussian and scaled Dawson parts directly, avoiding complex
        // divisions and the general complex-response setup at every sample.
        if self.electron_collision_s == 0. && state.ion_collision_s == 0. {
            let electron = CollisionlessResponse::new(
                k,
                ME,
                state.electron_temperature_k,
                alpha / state.electron_temperature_k,
            );
            if self.ions.len() == 1 {
                let ion = CollisionlessResponse::new(
                    k,
                    self.ions[0].mass_amu * AMU,
                    state.ion_temperature_k,
                    alpha * self.ions[0].fraction / state.ion_temperature_k,
                );
                let shift = k * state.velocity_m_s;
                #[cfg(target_arch = "x86_64")]
                if std::is_x86_feature_detected!("avx2") && frequency_hz.len() >= 4 {
                    // All pointers cover full chunks of four. The scalar tail
                    // handles 0--3 remaining samples without out-of-bounds reads.
                    return unsafe {
                        single_ion_avx2(&electron, &ion, self.ions[0].fraction, shift, frequency_hz)
                    };
                }
                return single_ion_scalar(
                    &electron,
                    &ion,
                    self.ions[0].fraction,
                    shift,
                    frequency_hz,
                );
            }
            let ions: Vec<_> = self
                .ions
                .iter()
                .map(|ion| {
                    (
                        CollisionlessResponse::new(
                            k,
                            ion.mass_amu * AMU,
                            state.ion_temperature_k,
                            alpha * ion.fraction / state.ion_temperature_k,
                        ),
                        ion.fraction,
                    )
                })
                .collect();
            let shift = k * state.velocity_m_s;
            return frequency_hz
                .iter()
                .map(|&f| {
                    let omega = 2. * PI * f - shift;
                    let (ce_re, ce_im, fe) = electron.at(omega);
                    let mut ci_re = 0.;
                    let mut ci_im = 0.;
                    let mut fi = 0.;
                    for (ion, fraction) in &ions {
                        let (re, im, free) = ion.at(omega);
                        ci_re += re;
                        ci_im += im;
                        fi += fraction * free;
                    }
                    let s = (((1. + ci_re).powi(2) + ci_im * ci_im) * fe
                        + (ce_re * ce_re + ce_im * ce_im) * fi)
                        / ((1. + ce_re + ci_re).powi(2) + (ce_im + ci_im).powi(2));
                    ensure!(s.is_finite() && s >= 0., "invalid model spectrum");
                    Ok(s)
                })
                .collect();
        }
        let mut output = Vec::with_capacity(frequency_hz.len());
        for &f in frequency_hz {
            let omega = 2. * PI * f - k * state.velocity_m_s;
            let (ae, fe) = response(
                omega,
                k,
                ME,
                state.electron_temperature_k,
                self.electron_collision_s,
            );
            let ce = ae * (alpha / state.electron_temperature_k);
            let mut ci = C::new(0., 0.);
            let mut fi = 0.;
            for ion in &self.ions {
                let (a, free) = response(
                    omega,
                    k,
                    ion.mass_amu * AMU,
                    state.ion_temperature_k,
                    state.ion_collision_s,
                );
                ci += a * (alpha * ion.fraction / state.ion_temperature_k);
                fi += ion.fraction * free;
            }
            let eps = C::new(1., 0.) + ce + ci;
            let s = ((C::new(1., 0.) + ci).norm_sqr() * fe + ce.norm_sqr() * fi) / eps.norm_sqr();
            ensure!(s.is_finite() && s >= 0., "invalid model spectrum");
            output.push(s);
        }
        Ok(output)
    }
}

fn single_ion_scalar(
    electron: &CollisionlessResponse,
    ion: &CollisionlessResponse,
    fraction: f64,
    shift: f64,
    frequency_hz: &[f64],
) -> Result<Vec<f64>> {
    frequency_hz
        .iter()
        .map(|&f| {
            let omega = 2. * PI * f - shift;
            let (er, ei, fe) = electron.at(omega);
            let (ir, ii, fi) = ion.at(omega);
            let result = (((1. + ir).powi(2) + ii * ii) * fe + (er * er + ei * ei) * fraction * fi)
                / ((1. + er + ir).powi(2) + (ei + ii).powi(2));
            ensure!(result.is_finite() && result >= 0., "invalid model spectrum");
            Ok(result)
        })
        .collect()
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn single_ion_avx2(
    electron: &CollisionlessResponse,
    ion: &CollisionlessResponse,
    fraction: f64,
    shift: f64,
    frequency_hz: &[f64],
) -> Result<Vec<f64>> {
    use std::arch::x86_64::*;
    let mut result = vec![0.; frequency_hz.len()];
    let count = frequency_hz.len() / 4 * 4;
    let one = _mm256_set1_pd(1.);
    let fraction = _mm256_set1_pd(fraction);
    for start in (0..count).step_by(4) {
        let mut er = [0.; 4];
        let mut ei = [0.; 4];
        let mut fe = [0.; 4];
        let mut ir = [0.; 4];
        let mut ii = [0.; 4];
        let mut fi = [0.; 4];
        for lane in 0..4 {
            let omega = 2. * PI * frequency_hz[start + lane] - shift;
            (er[lane], ei[lane], fe[lane]) = electron.at(omega);
            (ir[lane], ii[lane], fi[lane]) = ion.at(omega);
        }
        let er = _mm256_loadu_pd(er.as_ptr());
        let ei = _mm256_loadu_pd(ei.as_ptr());
        let ir = _mm256_loadu_pd(ir.as_ptr());
        let ii = _mm256_loadu_pd(ii.as_ptr());
        let fe = _mm256_loadu_pd(fe.as_ptr());
        let fi = _mm256_loadu_pd(fi.as_ptr());
        let ir1 = _mm256_add_pd(one, ir);
        let en = _mm256_add_pd(_mm256_mul_pd(er, er), _mm256_mul_pd(ei, ei));
        let inn = _mm256_add_pd(_mm256_mul_pd(ir1, ir1), _mm256_mul_pd(ii, ii));
        let dr = _mm256_add_pd(_mm256_add_pd(one, er), ir);
        let di = _mm256_add_pd(ei, ii);
        let num = _mm256_add_pd(
            _mm256_mul_pd(inn, fe),
            _mm256_mul_pd(_mm256_mul_pd(en, fraction), fi),
        );
        let den = _mm256_add_pd(_mm256_mul_pd(dr, dr), _mm256_mul_pd(di, di));
        _mm256_storeu_pd(result.as_mut_ptr().add(start), _mm256_div_pd(num, den));
    }
    if count < frequency_hz.len() {
        // SIMD fraction is a broadcast; use the original model fraction.
        let mut fractions = [0.; 4];
        _mm256_storeu_pd(fractions.as_mut_ptr(), fraction);
        let tail = single_ion_scalar(electron, ion, fractions[0], shift, &frequency_hz[count..])?;
        result[count..].copy_from_slice(&tail);
    }
    ensure!(
        result.iter().all(|s| s.is_finite() && *s >= 0.),
        "invalid model spectrum"
    );
    Ok(result)
}

struct CollisionlessResponse {
    inv_width: f64,
    susceptibility_scale: f64,
    free_scale: f64,
}
impl CollisionlessResponse {
    fn new(k: f64, mass: f64, temperature: f64, susceptibility_scale: f64) -> Self {
        let inv_width = 1. / (k * (2. * KB * temperature / mass).sqrt());
        Self {
            inv_width,
            susceptibility_scale,
            free_scale: 2. * PI.sqrt() * inv_width,
        }
    }
    #[inline]
    fn at(&self, omega: f64) -> (f64, f64, f64) {
        let x = omega * self.inv_width;
        let gaussian = (-x * x).exp();
        let a_re = 1. - PI.sqrt() * x * x.w_im();
        let a_im = PI.sqrt() * x * gaussian;
        (
            a_re * self.susceptibility_scale,
            a_im * self.susceptibility_scale,
            self.free_scale * gaussian,
        )
    }
}
fn response(omega: f64, k: f64, mass: f64, temperature: f64, collision: f64) -> (C, f64) {
    let a = k * (2. * KB * temperature / mass).sqrt();
    let g0 = C::new(omega / a, collision / a).w() * (PI.sqrt() / a);
    let g = g0 / (C::new(1., 0.) - collision * g0);
    (C::new(1., 0.) + C::new(0., omega) * g, 2. * g.re)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn real_axis_response_matches_complex_faddeeva_across_scales() {
        for temperature in [1., 300., 1000., 100000.] {
            for mass in [ME, AMU, 16. * AMU, 30.5 * AMU] {
                let width = 39. * (2. * KB * temperature / mass).sqrt();
                let prepared = CollisionlessResponse::new(39., mass, temperature, 1.);
                for x in [
                    -1000., -45., -10., -3.9, -1., -0.001, 0., 0.001, 1., 3.9, 10., 45., 1000.,
                ] {
                    let (expected, free) = response(x * width, 39., mass, temperature, 0.);
                    let (re, im, actual_free) = prepared.at(x * width);
                    assert!((C::new(re, im) - expected).norm() < 2e-12);
                    assert!((actual_free - free).abs() < 2e-12 / width);
                }
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_matches_scalar_for_short_shifted_grids_and_all_tails() {
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let k = 4. * PI * 440.2e6 / 299792458.;
        let alpha = 1e11 * E * E / (EPS0 * KB * k * k);
        for (te, ti) in [(2000., 1000.), (500., 250.), (10000., 3000.)] {
            let electron = CollisionlessResponse::new(k, ME, te, alpha / te);
            let ion = CollisionlessResponse::new(k, 16. * AMU, ti, alpha / ti);
            let grid: Vec<_> = (-50..=50).map(|i| i as f64 * 10000.).collect();
            for length in 0..=19 {
                for offset in 0..4 {
                    let f = &grid[offset..offset + length];
                    let expected = single_ion_scalar(&electron, &ion, 1., 1234., f).unwrap();
                    let actual = unsafe { single_ion_avx2(&electron, &ion, 1., 1234., f) }.unwrap();
                    assert_eq!(actual.len(), expected.len());
                    for (a, b) in actual.iter().zip(expected) {
                        assert!((a - b).abs() <= 1e-13 * b.abs().max(1e-30), "{a} != {b}");
                    }
                }
            }
        }
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    #[ignore = "manual paired CPU-affinity performance check"]
    fn benchmark_single_ion_scalar_and_avx2() {
        use std::{hint::black_box, time::Instant};
        if !std::is_x86_feature_detected!("avx2") {
            return;
        }
        let k = 4. * PI * 440.2e6 / 299792458.;
        let alpha = 1e11 * E * E / (EPS0 * KB * k * k);
        let electron = CollisionlessResponse::new(k, ME, 2000., alpha / 2000.);
        let ion = CollisionlessResponse::new(k, 16. * AMU, 1000., alpha / 1000.);
        let grid: Vec<_> = (-400..=400).map(|i| i as f64 * 250.).collect();
        let mut scalar = Vec::new();
        let mut simd = Vec::new();
        for trial in 0..31 {
            for vector in if trial % 2 == 0 {
                [true, false]
            } else {
                [false, true]
            } {
                let start = Instant::now();
                for _ in 0..100 {
                    let output = if vector {
                        unsafe { single_ion_avx2(&electron, &ion, 1., 0., black_box(&grid)) }
                    } else {
                        single_ion_scalar(&electron, &ion, 1., 0., black_box(&grid))
                    };
                    black_box(output.unwrap());
                }
                let elapsed = start.elapsed().as_secs_f64() / 100.;
                if vector {
                    simd.push(elapsed);
                } else {
                    scalar.push(elapsed);
                }
            }
        }
        scalar.sort_by(f64::total_cmp);
        simd.sort_by(f64::total_cmp);
        println!(
            "scalar_median_s={} avx2_median_s={} simd_speedup={}",
            scalar[15],
            simd[15],
            scalar[15] / simd[15]
        );
    }
    #[test]
    fn maxwellian_response_matches_independent_time_quadrature() {
        // Independent Simpson integration of the free-particle correlation.
        for nu in [0., 300., 5000.] {
            let k = 39.;
            let t = 900.;
            let mass = 30.5 * AMU;
            let a = k * (2. * KB * t / mass).sqrt();
            for omega in [0., 2000., -8000., 16000.] {
                let n = 10000;
                let dt = 14. / a / n as f64;
                let mut g = C::new(0., 0.);
                for i in 0..=n {
                    let x = i as f64 * dt;
                    let weight = if i == 0 || i == n {
                        1.
                    } else if i % 2 == 0 {
                        2.
                    } else {
                        4.
                    };
                    g += C::from_polar(weight * (-a * a * x * x / 4. - nu * x).exp(), omega * x);
                }
                g *= dt / 3.;
                g /= C::new(1., 0.) - nu * g;
                let (chi, s) = response(omega, k, mass, t, nu);
                assert!((s - 2. * g.re).abs() < 1e-12);
                assert!((chi - (C::new(1., 0.) + C::new(0., omega) * g)).norm() < 1e-8);
            }
        }
    }
    #[test]
    fn spectrum_symmetry_drift_and_composition() {
        let m = SpectrumModel {
            radar_frequency_hz: 930e6,
            ions: vec![
                IonSpecies {
                    mass_amu: 16.,
                    fraction: 0.7,
                },
                IonSpecies {
                    mass_amu: 30.5,
                    fraction: 0.3,
                },
            ],
            electron_collision_s: 0.,
        };
        let p = PlasmaState {
            electron_density_m3: 3e11,
            ion_temperature_k: 900.,
            electron_temperature_k: 1800.,
            velocity_m_s: 0.,
            ion_collision_s: 100.,
        };
        let f = vec![-2000., 0., 2000.];
        let a = m.evaluate(p, &f).unwrap();
        assert!((a[0] - a[2]).abs() < 1e-15);
        let mut moving = p;
        moving.velocity_m_s = 123.;
        let shift = m.wavenumber() * 123. / (2. * PI);
        let b = m
            .evaluate(moving, &f.iter().map(|x| x + shift).collect::<Vec<_>>())
            .unwrap();
        for (a, b) in a.iter().zip(b) {
            assert!((a - b).abs() < 1e-15);
        }
        let mut reordered = m.clone();
        reordered.ions.reverse();
        assert_eq!(
            m.evaluate(p, &f).unwrap(),
            reordered.evaluate(p, &f).unwrap()
        );
    }
}

#[cfg(test)]
mod normalization_tests {
    use super::*;
    #[test]
    fn ion_line_area_matches_equilibrium_debye_limit() {
        let m = SpectrumModel {
            radar_frequency_hz: 930e6,
            ions: vec![IonSpecies {
                mass_amu: 30.5,
                fraction: 1.,
            }],
            electron_collision_s: 0.,
        };
        let p = PlasmaState {
            electron_density_m3: 1e12,
            ion_temperature_k: 1000.,
            electron_temperature_k: 1000.,
            velocity_m_s: 0.,
            ion_collision_s: 0.,
        };
        let f: Vec<_> = (-8000..=8000).map(|i| i as f64 * 10.).collect();
        let s = m.evaluate(p, &f).unwrap();
        let area = s
            .iter()
            .enumerate()
            .map(|(i, value)| value * if i == 0 || i == s.len() - 1 { 5. } else { 10. })
            .sum::<f64>();
        let a = p.electron_density_m3 * E * E
            / (EPS0 * KB * p.electron_temperature_k * m.wavenumber().powi(2));
        let expected = a * a / ((1. + a) * (1. + 2. * a));
        // Finite electron mass and the included electron continuum make the
        // finite ion-band integral differ slightly from the zero-mass limit.
        assert!(
            (area / expected - 1.).abs() < 0.015,
            "area {area}, equilibrium {expected}"
        );
    }
}
