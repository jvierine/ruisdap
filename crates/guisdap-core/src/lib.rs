//! BEATA 2.2 quick-look signal processing, independent of file and network I/O.
pub mod calibration;
pub mod quality;
use anyhow::{ensure, Result};
use rustfft::{num_complex::Complex, FftPlanner};
use serde::{Deserialize, Serialize};
pub type C = Complex<f64>;
pub const DT: f64 = 0.4e-6;
pub const BANDS: usize = 2;
pub const GATES: usize = 90;
pub const NFFT: usize = 1600;
pub const STRIDE: usize = 82232;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FrameMeta {
    pub id: String,
    pub source: String,
    pub source_sha256: String,
    pub start_unix: f64,
    pub end_unix: f64,
    pub elevation_deg: f64,
    pub azimuth_deg: f64,
    pub centers_hz: [f64; 2],
    pub integration_s: f64,
    pub transmitter_w: f64,
    pub schema: String,
}

pub struct Products {
    pub native_long: Vec<f32>, // band, lag (1..800), gate, real/imag; normalized
    pub native_short: Vec<f32>, // band, lag (0..49), gate, real/imag; normalized
    pub native_uncoded: Vec<f32>, // band, lag (0..49), broad gate, real/imag
    pub broad: Vec<f32>,       // band, broad range gate, frequency (2,4,100)
    pub decoded: Vec<f32>,     // band, range gate, frequency (2,90,1600)
    pub decoded_lag0: Vec<f32>, // signed short-block lag-zero diagnostic; no extrapolation
    pub ion_power: Vec<f32>,
    pub ion_short_acf: Vec<f32>, // 431 gates, real/imag, measured 10 us lag
    pub ion_short_slices: Vec<f32>, // 14 subintegrations, 431 gates, real/imag
    pub ion_calibration: f64,
    pub ion_excess: Vec<f32>,
    pub plasma_power: Vec<f32>,
    pub plasma_excess: Vec<f32>,
    pub background: [f64; 3],
}

pub fn altitude_km(range: f64, elevation: f64) -> f64 {
    // Spherical Earth, 6370 km radius; site altitude 0.030 km (legacy site metadata).
    let r = 6370.0 + 0.030;
    (r * r + range * range + 2.0 * r * range * elevation.to_radians().sin()).sqrt() - 6370.0
}
pub fn ion_range_km(i: usize) -> f64 {
    (1548.0 * 0.4 + 10.0 * i as f64) * 0.149896229
}
/// Numerical BEATA response centroid, including all partial-code tails.
pub fn ion_short_range_km(i: usize) -> f64 {
    (810.5 * 0.4 + 10.0 * i as f64) * 0.149896229
}
pub fn plasma_range_km(i: usize) -> f64 {
    (955.0 - 320.0 - 0.4 + 0.4 * i as f64) * 0.149896229
}
pub fn decoded_range_km(i: usize) -> f64 {
    (955.0 - 240.0 - 20.4 + 20.0 * i as f64) * 0.149896229
}

/// Shifted real spectrum of a Hermitian ACF, with an explicit zero at the
/// unmeasured terminal lag. dt scaling is optional (broad GUISDAP convention).
pub fn spectrum(acf: &[C], dt_scale: f64, window: bool) -> Vec<f32> {
    let n = 2 * acf.len();
    let mut v = vec![C::new(0.0, 0.0); n];
    for (l, &a) in acf.iter().enumerate() {
        let w = if window {
            0.54 - 0.46
                * (2.0 * std::f64::consts::PI * (acf.len() + l) as f64 / (n - 1) as f64).cos()
        } else {
            1.0
        };
        v[l] = a * w;
        if l > 0 {
            v[n - l] = v[l].conj();
        }
    }
    v[0].im = 0.0;
    FftPlanner::new().plan_fft_forward(n).process(&mut v);
    (0..n)
        .map(|i| (v[(i + n / 2) % n].re * dt_scale) as f32)
        .collect()
}

// Fractional-lag alternating-code boxcar overlap weighting. See docs/science.md.
fn triangle(k: usize) -> f64 {
    let m = k % 100;
    (if m <= 50 { m } else { 100 - m }) as f64
}
fn filter_weight(k: usize) -> f64 {
    let t = triangle(k);
    if t == 50.0 {
        1.0
    } else {
        t / (50.0 - 1.0 / 3.0)
    }
}
pub fn decoded_weight(k: usize) -> f64 {
    if k == 0 {
        return 31.0 / 3.0;
    }
    let m1 = if k <= 50 {
        31.0 * (100 - k) as f64
    } else {
        triangle(k) * (31.0 - 2.0 * ((k - 1) / 100) as f64)
    };
    let second = if k >= 51 {
        filter_weight(k - 50) * triangle(k - 50) * (30.0 - 2.0 * ((k - 51) / 100) as f64)
    } else {
        0.0
    };
    filter_weight(k) * m1 + second
}

pub fn process_ion(data: &[C], normalization: f64) -> Result<Products> {
    ensure!(
        data.len() == 196493,
        "unsupported BEATA layout: {} elements, expected 196493",
        data.len()
    );
    ensure!(
        normalization.is_finite() && normalization > 0.0,
        "invalid accumulation normalization"
    );
    ensure!(
        data.iter().all(|z| z.re.is_finite() && z.im.is_finite()),
        "non-finite correlator data"
    );
    let mean = |offset: usize| {
        data[offset..offset + 16].iter().map(|z| z.re).sum::<f64>() / 16.0 / normalization
    };
    let background = [2.0 * mean(0), 2.0 * mean(32029), 2.0 * mean(32029 + STRIDE)];
    ensure!(
        background.iter().all(|b| *b > 0.0),
        "non-positive noise reference"
    );
    let ion_power: Vec<f32> = data[48..480]
        .iter()
        .map(|z| (z.re / normalization) as f32)
        .collect();
    let ion_excess = ion_power
        .iter()
        .map(|p| (*p as f64 / background[0] - 1.0) as f32)
        .collect();
    let p = Products {
        native_long: vec![],
        native_short: vec![],
        native_uncoded: vec![],
        broad: vec![],
        decoded: vec![],
        decoded_lag0: vec![],
        ion_power,
        ion_short_acf: data[480..911]
            .iter()
            .flat_map(|z| [(z.re / normalization) as f32, (z.im / normalization) as f32])
            .collect(),
        ion_short_slices: data[25531..31565]
            .iter()
            .flat_map(|z| [(z.re / normalization) as f32, (z.im / normalization) as f32])
            .collect(),
        ion_calibration: 2.0 * mean(16),
        ion_excess,
        plasma_power: vec![],
        plasma_excess: vec![],
        background,
    };
    Ok(p)
}

pub fn process(data: &[C], normalization: f64) -> Result<Products> {
    let mut p = process_ion(data, normalization)?;
    for band in 0..BANDS {
        let offset = STRIDE * band;
        for (dest, start, len) in [
            (&mut p.native_long, 36611, 800 * 90),
            (&mut p.native_short, 32061, 50 * 91),
            (&mut p.native_uncoded, 108611, 50 * 4),
        ] {
            for z in &data[start + offset..start + offset + len] {
                dest.extend([(z.re / normalization) as f32, (z.im / normalization) as f32]);
            }
        }
        for gate in 0..4 {
            let acf: Vec<C> = (0..50)
                .map(|l| data[108611 + offset + l * 4 + gate] / normalization)
                .collect();
            p.broad.extend(spectrum(&acf, 1.0, false));
        }
        for gate in 0..GATES {
            let z = data[32061 + offset + gate] + data[32061 + offset + gate + 1];
            p.decoded_lag0
                .push((z.re / (decoded_weight(0) * normalization)) as f32);
        }
        for i in 0..5450 {
            let power = data[108811 + offset + i].re / normalization;
            p.plasma_power.push(power as f32);
            p.plasma_excess
                .push((power / p.background[band + 1] - 1.0) as f32);
        }
    }
    p.decoded = native_spectra(&p.native_long, true)?;
    Ok(p)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn positive_frequency_sign_and_parseval() {
        let acf: Vec<C> = (0..50)
            .map(|l| C::from_polar(1.0, 2.0 * std::f64::consts::PI * 0.1 * l as f64))
            .collect();
        let s = spectrum(&acf, 1.0, false);
        assert_eq!(
            s.iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap()
                .0,
            60
        );
        assert!((s.iter().sum::<f32>() / 100.0 - 1.0).abs() < 1e-5);
    }
    #[test]
    fn geometry_is_not_slant_range() {
        assert!((altitude_km(100.0, 90.0) - 100.03).abs() < 1e-8);
        assert!(altitude_km(100.0, 77.8) < 100.0);
        assert!(altitude_km(100.0, 0.0) > 0.03);
    }
    #[test]
    fn weight_is_positive_and_periodic_overlap_matches() {
        for k in 0..801 {
            assert!(decoded_weight(k) > 0.0);
        }
        assert_eq!(decoded_weight(50), 1550.0);
        assert_eq!(decoded_weight(100), 1500.0);
        assert_eq!(decoded_weight(800), 800.0);
    }
    #[test]
    fn unsupported_layout_fails() {
        assert!(process(&[], 1.0).is_err());
    }
}

/// Direct FFT of every recorded long lag. The absent R(0) is explicitly zero,
/// so the signed output is spectral contrast, not absolute power.
/// Optional positive sqrt-overlap taper follows the user's ipypl quick-look.
pub fn native_spectra(data: &[f32], weighted: bool) -> Result<Vec<f32>> {
    ensure!(
        data.len() == 2 * 800 * 90 * 2,
        "invalid native ACF dimensions"
    );
    let mut output = Vec::with_capacity(2 * 90 * NFFT);
    for band in 0..2 {
        for gate in 0..90 {
            let mut acf = vec![C::new(0., 0.); 801];
            for (lag, z) in acf.iter_mut().enumerate().skip(1) {
                let i = ((band * 800 + lag - 1) * 90 + gate) * 2;
                let w = if weighted {
                    ((640. - 0.4 * (lag - 1) as f64) / 640.).sqrt()
                } else {
                    1.
                };
                *z = C::new(data[i] as f64, data[i + 1] as f64) * w;
            }
            // plwin performs a circular 1600-point inverse FFT. R(800) is
            // its self-conjugate midpoint, so include it once, not twice.
            // Preserve its recorded imaginary residual in native_long, but
            // a real spectrum uses only its real component at this endpoint.
            let mut circular = vec![C::new(0., 0.); NFFT];
            circular[1..=800].copy_from_slice(&acf[1..]);
            circular[800].im = 0.;
            for lag in 1..800 {
                circular[NFFT - lag] = acf[lag].conj();
            }
            FftPlanner::new()
                .plan_fft_forward(NFFT)
                .process(&mut circular);
            output.extend((0..NFFT).map(|i| (circular[(i + NFFT / 2) % NFFT].re * DT) as f32));
        }
    }
    Ok(output)
}

/// Merge two band/row/frequency arrays onto a common increasing frequency grid.
/// Linear interpolation within each receiver, equal-weight mean in overlap.
/// Unsupported frequencies are NaN, never extrapolated.
pub fn merge_bands(
    data: &[f32],
    rows: usize,
    bins: usize,
    centers: [f64; 2],
    df: f64,
) -> Result<(Vec<f32>, usize, f64)> {
    ensure!(
        bins > 1 && rows > 0 && data.len() == 2 * rows * bins && df > 0.,
        "invalid band grid"
    );
    let starts = centers.map(|c| c - (bins / 2) as f64 * df);
    let low = starts[0].min(starts[1]);
    let high = starts[0].max(starts[1]) + (bins - 1) as f64 * df;
    let width = ((high - low) / df).floor() as usize + 1;
    let mut out = vec![f32::NAN; rows * width];
    for row in 0..rows {
        for k in 0..width {
            let f = low + k as f64 * df;
            let mut sum = 0.;
            let mut count = 0;
            for (b, start) in starts.iter().enumerate() {
                let q = (f - start) / df;
                if q >= -1e-8 && q <= (bins - 1) as f64 + 1e-8 {
                    let q = q.clamp(0., (bins - 1) as f64);
                    let i = q.floor() as usize;
                    let j = (i + 1).min(bins - 1);
                    let t = (q - i as f64) as f32;
                    let offset = (b * rows + row) * bins;
                    sum += data[offset + i] * (1. - t) + data[offset + j] * t;
                    count += 1;
                }
            }
            if count > 0 {
                out[row * width + k] = sum / count as f32;
            }
        }
    }
    Ok((out, width, low))
}

#[cfg(test)]
mod direct_tests {
    use super::*;
    #[test]
    fn native_lag_index_sign_and_last_lag_are_preserved() {
        let mut d = vec![0.; 2 * 800 * 90 * 2];
        for l in 1..=800 {
            let z = C::from_polar(
                1.,
                2. * std::f64::consts::PI * 100. * l as f64 / NFFT as f64,
            );
            let i = (l - 1) * 90 * 2;
            d[i] = z.re as f32;
            d[i + 1] = z.im as f32;
        }
        let s = native_spectra(&d, false).unwrap();
        assert_eq!(
            s[..NFFT]
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.total_cmp(b.1))
                .unwrap()
                .0,
            900
        );
        assert!(s[..NFFT].iter().sum::<f32>().abs() < 1e-6);
        d.fill(0.);
        d[(799 * 90 + 89) * 2] = 1.;
        let s = native_spectra(&d, true).unwrap();
        assert!(s[89 * NFFT..90 * NFFT].iter().any(|v| v.abs() > 1e-8));
        assert!(s[..89 * NFFT].iter().all(|v| *v == 0.));
    }
    #[test]
    fn circular_midpoint_is_counted_once() {
        let mut d = vec![0.; 2 * 800 * 90 * 2];
        d[799 * 90 * 2] = 3.;
        let s = native_spectra(&d, false).unwrap();
        for (i, x) in s[..NFFT].iter().enumerate() {
            let expected = if i % 2 == 0 { 3. * DT } else { -3. * DT };
            assert!((*x as f64 - expected).abs() < 1e-12);
        }
    }
    #[test]
    fn short_and_uncoded_blocks_never_enter_long_lag_spectrum() {
        let mut d = vec![C::new(1., 0.); 196493];
        let baseline = process(&d, 1.).unwrap();
        for band in 0..2 {
            for z in &mut d[32061 + band * STRIDE..36611 + band * STRIDE] {
                *z = C::new(137., -53.);
            }
            for z in &mut d[108611 + band * STRIDE..108811 + band * STRIDE] {
                *z = C::new(293., 71.);
            }
        }
        let changed = process(&d, 1.).unwrap();
        assert_eq!(baseline.native_long, changed.native_long);
        assert_eq!(baseline.decoded, changed.decoded);
        assert_ne!(baseline.native_short, changed.native_short);
        assert_ne!(baseline.native_uncoded, changed.native_uncoded);
        assert_ne!(baseline.broad, changed.broad);
        assert_eq!(changed.native_uncoded.len(), 2 * 50 * 4 * 2);
    }
    #[test]
    fn merging_preserves_frequency_positions_and_averages_overlap() {
        let (v, w, lo) =
            merge_bands(&[10., 20., 30., 40., 1., 2., 3., 4.], 1, 4, [2., 0.], 1.).unwrap();
        assert_eq!((w, lo), (6, -2.));
        assert_eq!(v, vec![1., 2., 6.5, 12., 30., 40.]);
        let (v, _, _) =
            merge_bands(&[0., 1., 2., 3., 0., 1., 2., 3.], 1, 4, [0.5, 0.], 1.).unwrap();
        assert!((v[1] - 0.75).abs() < 1e-6);
    }
}

/// Zero-lag raw-Ne estimate from the ion-line receiver with its own noise and calibration.
/// Raw Ne means Te/Ti=1 and no Debye correction, as in GUISDAP power_prof.
/// Signed subtraction is linear conditional on the estimated calibration gain;
/// this is a response-weighted profile, not an inversion of its 96-km kernel.
pub fn ion_raw_ne(
    p: &Products,
    cal_k: f64,
    tx_w: f64,
    pulses: f64,
) -> Result<(Vec<f32>, Vec<f32>)> {
    ensure!(
        cal_k > 0. && tx_w > 0. && pulses > 1.,
        "invalid ion calibration metadata"
    );
    let b = p.background[0];
    let c = p.ion_calibration;
    let gain = c - b;
    ensure!(gain > 0., "non-positive calibration increment");
    let effective_length = 7851.068411222712;
    let radar = ion_radar_factor(tx_w);
    // Sixteen background/calibration range samples, half the pulses in each.
    // Account for covariance between samples caused by the receiver filter.
    let corr = [
        0.18246371104456985f64,
        0.0005802602186175348,
        1.2516622566827527e-10,
    ];
    let corr_factor = 1.
        + 2. * corr
            .iter()
            .enumerate()
            .map(|(i, r)| (1. - (i + 1) as f64 / 16.) * r * r)
            .sum::<f64>();
    let reference_samples = (pulses / 2.) * 16. / corr_factor;
    let vb = b * b / (reference_samples + 1.);
    let vc = c * c / (reference_samples + 1.);
    let mut ne = Vec::with_capacity(432);
    let mut se = Vec::with_capacity(432);
    for (i, &power) in p.ion_power.iter().enumerate() {
        let total = power as f64;
        let net = total - b;
        let r = ion_range_km(i) * 1000.;
        let factor = 2. * cal_k * r * r / (radar * effective_length);
        ne.push((factor * net / gain) as f32);
        // GUISDAP var_prof approach: estimate variance from nearby range
        // samples after range normalization, instead of guessing the unknown
        // online clutter-filter degrees of freedom. Six samples, with finite
        // sample and receiver-correlation corrections. Real spatial structure
        // can inflate this estimate; this limitation is carried in provenance.
        let begin = i.saturating_sub(2).min(p.ion_power.len() - 6);
        let values: Vec<f64> = (begin..begin + 6)
            .map(|j| p.ion_power[j] as f64 * (ion_range_km(j) / ion_range_km(i)).powi(2))
            .collect();
        let mean = values.iter().sum::<f64>() / 6.;
        let dof = 5.
            - 2. * corr
                .iter()
                .enumerate()
                .map(|(k, rho)| (6. - (k + 1) as f64) / 6. * rho * rho)
                .sum::<f64>();
        let vs = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / dof;
        // Delta-method propagation includes shared B in numerator/denominator.
        let var = vs / gain.powi(2)
            + (total - c).powi(2) * vb / gain.powi(4)
            + net.powi(2) * vc / gain.powi(4);
        se.push((factor * var.sqrt()) as f32);
    }
    Ok((ne, se))
}

#[cfg(test)]
mod power_tests {
    use super::*;
    #[test]
    fn raw_ne_preserves_signed_power_and_range_squared() {
        let mut d = vec![C::new(1., 0.); 196493];
        for z in &mut d[16..32] {
            z.re = 2.;
        }
        for z in &mut d[48..480] {
            z.re = 3.;
        }
        let p = process(&d, 1.).unwrap();
        let (ne, se) = ion_raw_ne(&p, 177.8, 1.2e6, 896.).unwrap();
        assert!(
            (ne[100] as f64 / ne[0] as f64 - (ion_range_km(100) / ion_range_km(0)).powi(2)).abs()
                < 1e-6
        );
        assert!(se.iter().all(|s| *s > 0.));
        let mut p = p;
        p.ion_power[0] = 1.;
        p.ion_power[1] = 2.;
        let (ne, _) = ion_raw_ne(&p, 177.8, 1.2e6, 896.).unwrap();
        assert!(ne[0] < 0.);
        assert_eq!(ne[1], 0.);
    }
}

/// Monostatic radar coefficient multiplied by range squared.
fn ion_radar_factor(tx_w: f64) -> f64 {
    // Numerical response constants from beataT_2init.mat, SHA in docs/science.md.
    let ap0 = 0.18370696050680618;
    let receiver_bandwidth = ap0 / DT;
    let scale_m = 299792458. * DT / 2.;
    let lambda = 299792458. / 930e6;
    let radar_gain = 10f64.powf(4.81);
    let electron_radius: f64 = 2.8179403262e-15;
    // 0.7 is the documented calibration multiplier in this event's GUISDAP run.
    electron_radius.powi(2) * radar_gain * lambda * lambda * 0.460 * 0.66 * scale_m * tx_w
        / 0.7
        / (4. * std::f64::consts::PI * 1.380649e-23 * receiver_bandwidth)
}

/// Fourteen consecutive 10-us-lag subintegrations, in raw electron-density units.
/// Their mean has the same calibration as the recorded five-second sum.
/// Equal-duration slice timing is assigned by the exporter from record bounds.
pub fn ion_short_slice_raw_ne(slices: &[f32], gain: f64, cal_k: f64, tx_w: f64) -> Result<Vec<f32>> {
    ensure!(slices.len() == 14 * 431 * 2 && gain > 0. && cal_k > 0. && tx_w > 0., "invalid short-lag slice calibration");
    let radar = ion_radar_factor(tx_w);
    Ok((0..14).flat_map(|j| (0..431).map(move |i| {
        let pairs = (i / 2 + 1).min(31) as f64;
        let r = ion_short_range_km(i) * 1000.;
        let factor = 2. * cal_k * r * r / (radar * 156.249999713498 * pairs);
        (14. * factor * slices[(j * 431 + i) * 2] as f64 / gain) as f32
    })).collect())
}

/// Signed 10-us raw-Ne equivalent using GUISDAP's peak spectral response.
/// Linear in measured lag conditional on calibration. This is not a measured
/// zero lag: decorrelation/Doppler and the absolute calibration remain model errors.
pub fn ion_short_raw_ne(
    p: &Products,
    cal_k: f64,
    tx_w: f64,
    pulses: f64,
) -> Result<(Vec<f32>, Vec<f32>)> {
    ensure!(
        cal_k > 0. && tx_w > 0. && pulses > 1.,
        "invalid ion calibration metadata"
    );
    ensure!(
        p.ion_short_acf.len() == 862 && p.ion_short_slices.len() == 14 * 862,
        "invalid ion short-lag dimensions"
    );
    let b = p.background[0];
    let c = p.ion_calibration;
    let gain = c - b;
    ensure!(gain > 0., "non-positive calibration increment");
    let rho = [
        0.18246371104456985f64,
        0.0005802602186175348,
        1.2516622566827527e-10,
    ];
    let corr = 1.
        + 2. * rho
            .iter()
            .enumerate()
            .map(|(i, r)| (1. - (i + 1) as f64 / 16.) * r * r)
            .sum::<f64>();
    let nref = pulses / 2. * 16. / corr;
    let gain_var = (b * b + c * c) / (nref + 1.);
    let radar = ion_radar_factor(tx_w);
    let mut ne = Vec::with_capacity(431);
    let mut sigma = Vec::with_capacity(431);
    for i in 0..431 {
        // lpg_ND = 64 * min(1 + floor(gate/2), 31). Calibration
        // references above have already been put on the ND=64 scale.
        let code_pairs = (i / 2 + 1).min(31) as f64;
        let value = p.ion_short_acf[2 * i] as f64;
        let slices: Vec<f64> = (0..14)
            .map(|j| p.ion_short_slices[(j * 431 + i) * 2] as f64)
            .collect();
        let mean = slices.iter().sum::<f64>() / 14.;
        // Variance of the SUM, not the mean. Approximate: online clutter
        // filtering may correlate subintegrations; no covariance is invented.
        let var_sum = 14. / 13. * slices.iter().map(|x| (x - mean).powi(2)).sum::<f64>();
        let r = ion_short_range_km(i) * 1000.;
        let factor = 2. * cal_k * r * r / (radar * 156.249999713498 * code_pairs);
        // lpg_bac=0 for this decoded nonzero lag: do not subtract raw noise.
        ne.push((factor * value / gain) as f32);
        sigma.push(
            (factor * (var_sum / gain.powi(2) + value * value * gain_var / gain.powi(4)).sqrt())
                as f32,
        );
    }
    Ok((ne, sigma))
}

#[cfg(test)]
mod short_power_tests {
    use super::*;
    #[test]
    fn partial_gate_normalization_signed_sum_and_uncertainty() {
        let mut d = vec![C::new(1., 0.); 196493];
        d[16..32].fill(C::new(2., 0.));
        for i in 0..431 {
            let pairs = (i / 2 + 1).min(31) as f64;
            let scale = pairs / ion_short_range_km(i).powi(2);
            for j in 0..14 {
                d[25531 + j * 431 + i] = C::new((j as f64 - 6.) * scale, 0.);
            }
            d[480 + i] = (0..14).map(|j| d[25531 + j * 431 + i]).sum();
        }
        let mut p = process(&d, 1.).unwrap();
        let (ne, se) = ion_short_raw_ne(&p, 177.8, 1.2e6, 896.).unwrap();
        let native = ion_short_slice_raw_ne(&p.ion_short_slices, p.ion_calibration - p.background[0], 177.8, 1.2e6).unwrap();
        for i in 0..431 {
            let mean = (0..14).map(|j| native[j * 431 + i] as f64).sum::<f64>() / 14.;
            assert!((mean / ne[i] as f64 - 1.).abs() < 2e-6);
        }
        for v in &ne {
            assert!((*v / ne[0] - 1.).abs() < 1e-6);
        }
        assert!(se.iter().all(|v| *v > 0.));
        let h = altitude_km(ion_short_range_km(0), 77.8);
        assert!(h > 47. && h < 48.);
        p.ion_short_acf[0] *= -1.;
        p.ion_short_acf[2] = 0.;
        let (signed, _) = ion_short_raw_ne(&p, 177.8, 1.2e6, 896.).unwrap();
        assert_eq!(signed[0], -ne[0]);
        assert_eq!(signed[1], 0.);
    }
}
