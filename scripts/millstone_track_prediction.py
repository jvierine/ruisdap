"""Experimental training-only coherent echo forecast from amplitude metadata.

The queried pulse is removed before association, geometry fitting or phase
search. This module does not read RF or change production processing.
"""
import math

import numpy as np
from scipy.optimize import minimize


def fft_grid_search(u, b, gram_inv, clock, chirps):
    """Exact sparse amplitude-sum search on the measured pulse clock.

    Missing epochs have zero *weight* in these metadata sums. No RF samples
    are synthesized and the normalization uses only actual training events.
    """
    integer = np.rint((u-u.min())*1e6/clock).astype(int)
    np.testing.assert_allclose(u-u.min(), integer*clock/1e6, atol=1e-12, rtol=0)
    nfft = 1 << int(np.ceil(np.log2(4*(integer.max()+1))))
    frequencies = np.fft.fftfreq(nfft, clock/1e6)
    best = (-math.inf, 0., 0.)
    for start in range(0, len(chirps), 64):
        values = chirps[start:start+64]
        modulated = b*np.exp(-1j*np.pi*values[:, None]*u*u)
        z0 = np.zeros((len(values), nfft), complex)
        z1 = np.zeros_like(z0)
        z0[:, integer] = modulated
        z1[:, integer] = modulated*u
        r0, r1 = np.fft.fft(z0, axis=1), np.fft.fft(z1, axis=1)
        explained = (gram_inv[0, 0]*abs(r0)**2+gram_inv[1, 1]*abs(r1)**2
                     +2*gram_inv[0, 1]*np.real(np.conj(r0)*r1))
        ci, fi = np.unravel_index(np.argmax(explained), explained.shape)
        if explained[ci, fi] > best[0]:
            best = (float(explained[ci, fi]), float(frequencies[fi]), float(values[ci]))
    return best, nfft


def groups(events, held_sample, radius_s=.35):
    """Associate other pulses locally without looking at the queried pulse."""
    rows = [e for e in events if e['sample'] != held_sample and
            abs(e['sample']-held_sample) <= radius_s*1e6 and e['snr_db'] >= 14]
    parent = list(range(len(rows)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, a in enumerate(rows):
        for j in range(i):
            b = rows[j]
            dt = abs(a['sample']-b['sample'])/1e6
            # Both signs are allowed here: the stored physical velocity sign
            # convention is still under review. This is association only.
            drift = max(abs(a['doppler_hz']), abs(b['doppler_hz']))/440.2e6*1e6*dt
            if abs(a['delay_samples']-b['delay_samples']) <= 16+drift and abs(a['doppler_hz']-b['doppler_hz']) <= 2000:
                parent[root(i)] = root(j)
    found = {}
    for i, row in enumerate(rows):
        found.setdefault(root(i), []).append(row)
    return list(found.values())


def forecast(rows, held_sample, tx_centre, min_pulses=16, max_pulses=48):
    # Defence in depth: never rely on the caller to exclude the held pulse.
    rows = [r for r in rows if r['sample'] != held_sample]
    # One component per training pulse; unresolved associations are reported as
    # unsupported instead of silently choosing with held-pulse information.
    if len({r['sample'] for r in rows}) != len(rows):
        return dict(supported=False, reason='multiple_components_in_training_pulse')
    rows = sorted(rows, key=lambda r: (abs(r['sample']-held_sample), r['sample']))[:max_pulses]
    if len(rows) < min_pulses:
        return dict(supported=False, reason='insufficient_other_pulses', training_count=len(rows))
    rows.sort(key=lambda r: r['sample'])
    keys = np.array([r['sample'] for r in rows], np.int64)
    clock = int(np.gcd.reduce(np.diff(keys)))
    if clock <= 0 or (held_sample-int(keys[0])) % clock:
        return dict(supported=False, reason='held_epoch_not_on_training_clock')
    # Avoid constructing a huge search if the experiment's timestamps do not
    # provide a stable pulse clock. No unsupported timing is guessed.
    if clock < 1000:
        return dict(supported=False, reason='unsupported_irregular_clock')
    u = (keys-held_sample)/1e6
    span = np.ptp(u)
    if span <= 0:
        return dict(supported=False, reason='zero_training_span')
    delay = np.array([r['delay_samples'] for r in rows])
    frequency = np.array([r['doppler_hz'] for r in rows])
    amplitude = np.array([r['amplitude'] for r in rows], complex)
    # Refer each fit to its echo centre. Fit phase at the pulse epoch rather
    # than treating a noisy single-pulse Doppler estimate as an exact phase.
    b = amplitude*np.exp(2j*np.pi*frequency*(tx_centre+delay)/1e6)
    x = np.column_stack([np.ones(len(u)), u])
    gram_inv = np.linalg.inv(x.T@x)
    dc = 1/(2*span**2)
    chirps = np.arange(-1000., 1000.+dc/2, dc)
    nfft = 1 << int(np.ceil(np.log2(4*(int(round(span*1e6/clock))+1))))
    if nfft*len(chirps) > 2_000_000:
        return dict(supported=False, reason='phase_search_budget_exceeded')
    best, nfft = fft_grid_search(u, b, gram_inv, clock, chirps)

    def solve(parameters):
        theta = 2*np.pi*(parameters[0]*u+.5*parameters[1]*u*u)
        rhs = x.T@(b*np.exp(-1j*theta))
        coefficients = gram_inv@rhs
        return coefficients, float(np.vdot(b, b).real-np.vdot(rhs, coefficients).real)

    result = minimize(lambda z: solve(z)[1], [best[1], best[2]], method='Nelder-Mead',
                      options=dict(maxiter=150, xatol=1e-7, fatol=1e-14))
    parameters = result.x if result.success and solve(result.x)[1] <= solve([best[1], best[2]])[1] else np.array([best[1], best[2]])
    coefficients, residual = solve(parameters)
    geometry = np.column_stack([np.ones(len(u)), u, u*u])
    delay_fit = np.linalg.lstsq(geometry, delay, rcond=None)[0]
    frequency_fit = np.linalg.lstsq(geometry, frequency, rcond=None)[0]
    predicted_delay = float(delay_fit[0])
    predicted_frequency = float(frequency_fit[0])
    geometry_cov = np.linalg.inv(geometry.T@geometry)
    delay_std = float(np.sqrt(np.sum((delay-geometry@delay_fit)**2)/max(1, len(rows)-3)*geometry_cov[0, 0]))
    frequency_std = float(np.sqrt(np.sum((frequency-geometry@frequency_fit)**2)/max(1, len(rows)-3)*geometry_cov[0, 0]))
    predicted_amplitude = coefficients[0]*np.exp(-2j*np.pi*predicted_frequency*(tx_centre+predicted_delay)/1e6)
    energy = float(np.vdot(b, b).real)
    return dict(supported=True, reason='training_forecast', training_samples=keys.tolist(), training_count=len(rows),
                pulse_clock_samples=clock, delay_samples=predicted_delay, doppler_hz=predicted_frequency,
                delay_prediction_std_samples=delay_std, doppler_prediction_std_hz=frequency_std,
                delay_rate_samples_s=float(delay_fit[1]), doppler_slope_hz_s=float(frequency_fit[1]),
                doppler_from_range_rate_hz=float(-440.2e6*delay_fit[1]/1e6),
                amplitude=[float(predicted_amplitude.real), float(predicted_amplitude.imag)],
                centre_amplitude=[float(coefficients[0].real), float(coefficients[0].imag)],
                phase_frequency_alias_hz=float(parameters[0]), phase_chirp_hz_s=float(parameters[1]),
                explained_fraction=float(1-residual/energy) if energy else 0.,
                training_amplitude_residual_rms=float(np.sqrt(max(0., residual)/max(1, len(rows)-4))),
                prediction_variance=float(max(0., residual)/max(1, len(rows)-4)*gram_inv[0, 0]),
                search_cells=nfft*len(chirps), phase_fft_length=nfft)


def predict(events, held_sample, tx_centre, radius_s=.35):
    return [forecast(group, held_sample, tx_centre) for group in groups(events, held_sample, radius_s)]
