//! Weighted lag-profile inversion using f64 complex arithmetic. The unknowns
//! are narrow range gates and an optional range-independent background column.
use crate::{
    signal::{lag_products, lowpass, percentile},
    C,
};
use anyhow::{ensure, Result};
use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Settings {
    pub sample_rate_hz: f64,
    pub decimation: usize,
    pub lag_samples: Vec<usize>,
    pub range_start: usize,
    pub range_end: usize,
    pub measurement_start: usize,
    pub measurement_end: usize,
    pub tx_gate: [usize; 2],
    pub echo_gate: [usize; 2],
    pub noise_gate: [usize; 2],
    pub injection_gate: [usize; 2],
    pub filter_half: usize,
    pub passband_hz: f64,
    #[serde(default = "default_threshold")]
    pub outlier_threshold: f64,
    #[serde(default)]
    pub background_column: bool,
}
fn default_threshold() -> f64 {
    10.
}
pub struct Pulse {
    pub echo: Vec<C>,
    pub transmit: Vec<C>,
    pub pair: usize,
    pub usable: Vec<bool>,
}
pub struct Estimate {
    pub acf: Vec<C>,
    pub ground_clutter_acf: Vec<C>,
    pub variance: Vec<f64>,
    pub retained_fraction: Vec<f64>,
    pub condition_proxy: Vec<f64>,
    pub noise_power: f64,
    pub injection_power: f64,
    pub receiver_gain: f64,
    pub system_temperature_k: f64,
}
/// Column-scaled Hermitian WLS. Exactly unconstrained columns remain NaN;
/// non-positive-definite systems fail rather than silently regularizing.
pub fn solve(
    rows: &[Vec<(usize, C)>],
    y: &[C],
    sigma: &[f64],
    n: usize,
) -> Result<(Vec<C>, Vec<f64>, f64)> {
    ensure!(
        rows.len() == y.len() && y.len() == sigma.len() && n > 0,
        "WLS dimensions"
    );
    let mut scale = vec![0.; n];
    for (row, &s) in rows.iter().zip(sigma) {
        ensure!(s.is_finite() && s > 0., "nonpositive row sigma");
        for &(j, a) in row {
            ensure!(j < n, "column outside design");
            scale[j] += a.norm_sqr() / (s * s)
        }
    }
    let active: Vec<_> = (0..n).filter(|&j| scale[j] > 0.).collect();
    ensure!(!active.is_empty(), "unconstrained design");
    for s in &mut scale {
        *s = s.sqrt()
    }
    let mut map = vec![usize::MAX; n];
    for (i, &j) in active.iter().enumerate() {
        map[j] = i
    }
    let mut normal = DMatrix::<C>::zeros(active.len(), active.len());
    let mut rhs = DVector::<C>::zeros(active.len());
    for ((row, &y), &s) in rows.iter().zip(y).zip(sigma) {
        ensure!(
            y.re.is_finite() && y.im.is_finite(),
            "nonfinite measurement"
        );
        for &(j, a) in row {
            let a = a / (s * scale[j]);
            rhs[map[j]] += a.conj() * y / s;
            for &(k, b) in row {
                normal[(map[j], map[k])] += a.conj() * b / (s * scale[k]);
            }
        }
    }
    let factor = normal
        .cholesky()
        .ok_or_else(|| anyhow::anyhow!("rank-deficient lag design"))?;
    let diag: Vec<_> = factor.l().diagonal().iter().map(|v| v.norm_sqr()).collect();
    let ratio = diag.iter().copied().fold(f64::INFINITY, f64::min)
        / diag.iter().copied().fold(0., f64::max);
    ensure!(ratio > 1e-12, "ill-conditioned lag design");
    let answer = factor.solve(&rhs);
    let cov = factor.inverse();
    let mut x = vec![C::new(f64::NAN, f64::NAN); n];
    let mut variance = vec![f64::NAN; n];
    for (i, &j) in active.iter().enumerate() {
        x[j] = answer[i] / scale[j];
        variance[j] = cov[(i, i)].re / (scale[j] * scale[j])
    }
    Ok((x, variance, ratio))
}
fn mean_power(z: &[C], window: [usize; 2], dc: C) -> f64 {
    z[window[0]..window[1]]
        .iter()
        .map(|z| (*z - dc).norm_sqr())
        .sum::<f64>()
        / (window[1] - window[0]) as f64
}
pub fn estimate(pulses: &[Pulse], s: &Settings, injection_temperature_k: f64) -> Result<Estimate> {
    estimate_with_progress(pulses, s, injection_temperature_k, |_, _, _| Ok(()))
}
pub fn estimate_with_progress(
    pulses: &[Pulse],
    s: &Settings,
    injection_temperature_k: f64,
    mut progress: impl FnMut(usize, usize, &Estimate) -> Result<()>,
) -> Result<Estimate> {
    ensure!(
        pulses.len() >= 16
            && s.decimation > 0
            && s.range_end > s.range_start
            && s.measurement_end > s.measurement_start
            && !s.lag_samples.is_empty(),
        "insufficient pulses or invalid LPI geometry"
    );
    let n = pulses[0].echo.len();
    ensure!(
        pulses.iter().all(|p| p.echo.len() == n
            && p.transmit.len() == n
            && p.usable.len() == n
            && p.pair < pulses.len()),
        "pulse dimensions/pair index mismatch"
    );
    ensure!(
        s.lag_samples.iter().all(|&l| l < n)
            && [s.tx_gate, s.echo_gate, s.noise_gate, s.injection_gate]
                .iter()
                .all(|g| g[0] < g[1] && g[1] <= n)
            && s.measurement_end * s.decimation + s.lag_samples.iter().max().unwrap() <= n,
        "invalid gate timing"
    );
    let nr = s.range_end - s.range_start;
    let nl = s.lag_samples.len();
    let cols = nr + usize::from(s.background_column);
    let mut filtered = Vec::new();
    let mut txs = Vec::new();
    let mut masks = Vec::new();
    for p in pulses {
        filtered.push(lowpass(
            &p.echo,
            s.sample_rate_hz,
            1.2 * s.passband_hz,
            s.filter_half,
        )?);
        let mut tx = p.transmit.clone();
        for (i, z) in tx.iter_mut().enumerate() {
            if !(s.tx_gate[0]..s.tx_gate[1]).contains(&i) {
                *z = C::default()
            }
        }
        let norm = tx.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt();
        ensure!(norm > 0., "zero transmit energy");
        for z in &mut tx {
            *z /= norm
        }
        txs.push(tx);
        let mut mask = p.usable.clone();
        for i in 0..n {
            if !p.usable[i] {
                for m in mask
                    .iter_mut()
                    .take((i + s.filter_half + 1).min(n))
                    .skip(i.saturating_sub(s.filter_half))
                {
                    *m = false
                }
            }
        }
        for (i, m) in mask.iter_mut().enumerate() {
            if !(s.echo_gate[0] + s.filter_half..s.echo_gate[1].saturating_sub(s.filter_half))
                .contains(&i)
            {
                *m = false
            }
        }
        masks.push(mask);
    }
    let means: Vec<_> = filtered
        .iter()
        .map(|z| {
            z[s.noise_gate[0]..s.noise_gate[1]].iter().sum::<C>()
                / (s.noise_gate[1] - s.noise_gate[0]) as f64
        })
        .collect();
    let dc = C::new(
        percentile(&means.iter().map(|z| z.re).collect::<Vec<_>>(), 0.5),
        percentile(&means.iter().map(|z| z.im).collect::<Vec<_>>(), 0.5),
    );
    let noise = percentile(
        &filtered
            .iter()
            .map(|z| mean_power(z, s.noise_gate, dc))
            .collect::<Vec<_>>(),
        0.5,
    );
    let injection = percentile(
        &filtered
            .iter()
            .map(|z| mean_power(z, s.injection_gate, dc))
            .collect::<Vec<_>>(),
        0.5,
    );
    let gain = (injection - noise) / injection_temperature_k;
    ensure!(
        gain.is_finite() && gain > 0.,
        "noise injection does not establish positive receiver gain"
    );
    let mut out = Estimate {
        acf: vec![C::new(f64::NAN, f64::NAN); nr * nl],
        ground_clutter_acf: vec![C::new(f64::NAN, f64::NAN); nr * nl],
        variance: vec![f64::NAN; nr * nl],
        retained_fraction: vec![0.; nl],
        condition_proxy: vec![0.; nl],
        noise_power: noise,
        injection_power: injection,
        receiver_gain: gain,
        system_temperature_k: noise / gain,
    };
    for (li, &lag) in s.lag_samples.iter().enumerate() {
        let mut products = Vec::new();
        let mut difference = Vec::new();
        let mut ambiguities = Vec::new();
        for (pi, p) in pulses.iter().enumerate() {
            products.push(lag_products(&filtered[pi], lag, s.decimation)?);
            let zd: Vec<_> = filtered[pi]
                .iter()
                .zip(&filtered[p.pair])
                .map(|(a, b)| a - b)
                .collect();
            difference.push(lag_products(&zd, lag, s.decimation)?);
            ambiguities.push(lag_products(&txs[pi], lag, s.decimation)?);
        }
        let nm = s.measurement_end - s.measurement_start;
        let robust: Vec<_> = (s.measurement_start..s.measurement_end)
            .map(|m| {
                (2. * percentile(
                    &products.iter().map(|p| p[m].norm_sqr()).collect::<Vec<_>>(),
                    0.34,
                ))
                .sqrt()
            })
            .collect();
        let robust_g: Vec<_> = (s.measurement_start..s.measurement_end)
            .map(|m| {
                (2. * percentile(
                    &difference
                        .iter()
                        .map(|p| p[m].norm_sqr())
                        .collect::<Vec<_>>(),
                    0.34,
                ))
                .sqrt()
            })
            .collect();
        // Direct f64 moving power sum: never FFT a single-precision weight matrix
        // dominated by a bright satellite. Preserve the legacy ten-pulse window.
        let local: Vec<Vec<f64>> = (0..pulses.len())
            .map(|p| {
                (s.measurement_start..s.measurement_end)
                    .map(|m| {
                        ((0..10)
                            .map(|j| {
                                products[(p + pulses.len() + 5 - j % pulses.len()) % pulses.len()]
                                    [m]
                                    .norm_sqr()
                            })
                            .sum::<f64>()
                            / 10.)
                            .sqrt()
                    })
                    .collect()
            })
            .collect();
        let median = percentile(&local.iter().flatten().copied().collect::<Vec<_>>(), 0.5);
        ensure!(
            median > 0. && median.is_finite(),
            "invalid product noise weights"
        );
        let mut rows = Vec::new();
        let mut y = Vec::new();
        let mut yg = Vec::new();
        let mut sigma = Vec::new();
        for (pi, p) in pulses.iter().enumerate() {
            for (mi, m) in (s.measurement_start..s.measurement_end).enumerate() {
                let start = m * s.decimation;
                let end = (m + 1) * s.decimation;
                let valid = (start..end).all(|i| {
                    masks[pi][i] && masks[pi][i + lag] && masks[p.pair][i] && masks[p.pair][i + lag]
                });
                let weight = local[pi][mi].max(median);
                if !valid
                    || products[pi][m].norm() > s.outlier_threshold * robust[mi]
                    || difference[pi][m].norm() > s.outlier_threshold * robust_g[mi]
                    || weight > 100. * median
                {
                    continue;
                }
                let mut row = Vec::new();
                for ri in s.range_start..s.range_end {
                    if m >= ri {
                        let j = m - ri;
                        if j < ambiguities[pi].len() {
                            let a = ambiguities[pi][j];
                            if a.norm_sqr() > 1e-30 {
                                row.push((ri - s.range_start, a))
                            }
                        }
                    }
                }
                if s.background_column {
                    row.push((nr, C::new(1., 0.)))
                }
                if row.is_empty() {
                    continue;
                }
                rows.push(row);
                y.push(products[pi][m]);
                yg.push(difference[pi][m]);
                sigma.push(weight);
            }
        }
        out.retained_fraction[li] = rows.len() as f64 / (pulses.len() * nm) as f64;
        if rows.len() < cols {
            continue;
        }
        let result = solve(&rows, &y, &sigma, cols);
        let result_g = solve(&rows, &yg, &sigma, cols);
        if let (Ok((x, var, condition)), Ok((xg, _, _))) = (result, result_g) {
            out.condition_proxy[li] = condition;
            for ri in 0..nr {
                out.acf[ri * nl + li] = x[ri] / gain;
                out.ground_clutter_acf[ri * nl + li] = xg[ri] / (2. * gain);
                out.variance[ri * nl + li] = var[ri] / gain.powi(2);
            }
        }
        progress(li + 1, nl, &out)?;
    }
    Ok(out)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scaled_wls_and_unconstrained() {
        let rows = vec![
            vec![(0, C::new(1e-9, 0.)), (1, C::new(0., 1e8))],
            vec![(0, C::new(2e-9, 0.)), (1, C::new(1e8, 0.))],
            vec![(0, C::new(3e-9, 0.)), (1, C::new(-1e8, 0.))],
        ];
        let truth = [C::new(2., -3.), C::new(4., 1.)];
        let y: Vec<_> = rows
            .iter()
            .map(|r| r.iter().map(|(j, a)| truth[*j] * a).sum())
            .collect();
        let (x, v, _) = solve(&rows, &y, &[1.; 3], 3).unwrap();
        assert!((x[1] - truth[1]).norm() < 1e-12);
        assert!(x[2].re.is_nan() && v[2].is_nan());
    }
    #[test]
    fn singular_design_is_rejected() {
        let rows = vec![vec![(0, C::new(1., 0.)), (1, C::new(1., 0.))]; 4];
        assert!(solve(&rows, &[C::new(1., 0.); 4], &[1.; 4], 2).is_err());
    }
}
