#!/usr/bin/env python3
"""Oracle control: estimate a known coherent target's amplitude in other pulses.

Known injected delay, Doppler and inter-pulse phase are supplied. This isolates
same-pulse projection bias; it is not a working target tracker or a campaign
false-alarm test. Synthetic voltages stay in memory and anonymous pipes only.
"""
import argparse
import hashlib
import json
import pathlib

import h5py
import numpy as np


def leave_one_out(coefficients):
    coefficients = np.asarray(coefficients)
    assert coefficients.ndim == 2 and len(coefficients) >= 2 and np.isfinite(coefficients).all()
    return (coefficients.sum(axis=0)-coefficients)/(len(coefficients)-1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=pathlib.Path, required=True)
    p.add_argument('--spectrum', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--pulses', type=int, default=16)
    p.add_argument('--groups', type=int, default=8)
    a = p.parse_args()
    assert 4 <= a.pulses <= 256 and 2 <= a.groups <= 256 and not a.output.exists()
    import importlib.util
    path = pathlib.Path(__file__).with_name('check-millstone-residual-controls.py')
    spec = importlib.util.spec_from_file_location('independent_control', path)
    control = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(control)
    interpolation = control.load('independent_template', path.with_name('check-millstone-residual-cfar.py'))
    summary = control.load('independent_summary', path.with_name('summarize-millstone-residual-controls.py'))
    with h5py.File(a.spectrum) as h:
        psd = h['ion_line_psd'][:]
        assert h.attrs['sample_rate_hz'] == 1e6
    native = control.NativeCancellation(a.binary.resolve())
    n, delay = 10000, round(550e3*2/299792458*1e6)
    lags = np.arange(0, 451, 10)
    lag10 = 1
    try:
        with h5py.File(a.output, 'x') as report:
            report.attrs['complete'] = False
            report.attrs['scope'] = __doc__
            report.attrs['phase_assumption'] = 'All injected pulses have the same already-aligned complex target amplitudes. Geometry is known; no CFAR selection or parameter refinement. Leave-one-out estimator uses every other independent pulse, without conditioning on validation acceptance.'
            report.attrs['interval_definition'] = 'Paired bootstrap across independent groups of pulses; pulses within a leave-one-out group are not treated as independent bootstrap units.'
            for name, filename in [('script', pathlib.Path(__file__)), ('native_binary', a.binary), ('spectrum', a.spectrum), ('pipe_helper', path), ('interpolation', pathlib.Path(interpolation.__file__)), ('summary', pathlib.Path(summary.__file__))]:
                report.attrs[name+'_sha256'] = hashlib.sha256(filename.read_bytes()).hexdigest()
            report['lags_samples'] = lags
            for family, tx_end, noise_gate in [('coded', 624, [7700, 8200]), ('long', 645, [7200, 7700])]:
                tx = np.zeros(n, complex)
                body = np.zeros(tx_end-76, complex)
                code = np.array([1, 1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1])
                body[:480] = np.repeat(code, 30) if family == 'coded' else 1
                tx[76:tx_end] = body
                for width in [0, 240]:
                    offsets = np.array([0]) if width == 0 else np.round(np.linspace(-3*width, 3*width, 49)).astype(int)
                    weights = np.ones(1) if width == 0 else np.exp(-.5*(offsets/width)**2)
                    weights /= weights.sum()
                    profile = np.zeros(n)
                    for offset, weight in zip(offsets, weights):
                        start = delay+76+int(offset)
                        profile[start:start+len(body)] += weight*abs(body)**2
                    support = profile > profile.max()*1e-6
                    for kind in ['single', 'double']:
                        detections = [dict(delay_samples=float(delay), doppler_hz=1000.)]
                        true_amplitudes = [3*np.exp(.3j)]
                        if kind == 'double':
                            detections.append(dict(delay_samples=delay+60.25, doppler_hz=2500.))
                            true_amplitudes.append(.6*np.exp(-.7j))
                        target = sum(amplitude*interpolation.template(tx, d['delay_samples'], d['doppler_hz'], 1e6, 16) for amplitude, d in zip(true_amplitudes, detections))
                        for power in [.03, .1, .2]:
                            group = report.create_group(f'{family}/width_{width}/{kind}_{power}')
                            group.attrs['nominal_diffuse_receiver_power'] = power
                            group.attrs['pulses'] = a.pulses
                            group.attrs['groups'] = a.groups
                            base_means, independent_means, truth_means = [], [], []
                            for replicate in range(a.groups):
                                subgroup = group.create_group(str(replicate+1))
                                frames = []
                                coefficients = []
                                for trial in range(a.pulses):
                                    rng = np.random.default_rng(replicate*a.pulses+trial+1)
                                    receiver = (rng.standard_normal(n)+1j*rng.standard_normal(n))/np.sqrt(2)
                                    diffuse = np.zeros(n, complex)
                                    for offset, weight in zip(offsets, weights):
                                        z = (rng.standard_normal(len(psd))+1j*rng.standard_normal(len(psd)))/np.sqrt(2)
                                        z = np.fft.ifft(z*np.sqrt(psd))*len(psd)/np.sqrt(psd.sum())
                                        start = delay+76+int(offset)
                                        diffuse[start:start+len(body)] += body*z[np.arange(start, start+len(body)) % len(z)]*np.sqrt(power*weight)
                                    ideal = receiver+diffuse+100*tx
                                    raw = ideal+target
                                    measured = np.zeros(n, complex)
                                    measured[76:tx_end] = raw[76:tx_end]
                                    fixed = dict(refine_delay_samples=0., refine_doppler_hz=0., min_amplitude_snr_db=0., shrink_weak_amplitudes=False)
                                    clean, mask, qa = native.cancel(raw, measured, detections, noise_gate, fixed)
                                    assert 'error' not in qa, qa
                                    coefficients.append([complex(*v) for v in qa['amplitudes']])
                                    frames.append((raw, ideal, diffuse, measured, clean, qa))
                                coefficients = np.asarray(coefficients)
                                predictions = leave_one_out(coefficients)
                                base, independent, truth = [], [], []
                                for trial, (raw, ideal, diffuse, measured, clean, qa) in enumerate(frames):
                                    predicted = raw.copy()
                                    for coefficient, d in zip(predictions[trial], detections):
                                        predicted -= coefficient*interpolation.template(measured, d['delay_samples'], d['doppler_hz'], 1e6, 16)
                                    row = subgroup.create_group(str(trial+1))
                                    reference = control.moments(ideal, support, lags)
                                    row['same_pulse_acf_error'] = control.moments(clean, support, lags)-reference
                                    row['independent_acf_error'] = control.moments(predicted, support, lags)-reference
                                    row['diffuse_acf'] = control.moments(diffuse, support, lags)
                                    row['predicted_amplitudes'] = predictions[trial]
                                    row['same_pulse_diagnostic_json'] = np.frombuffer(json.dumps(qa).encode(), np.uint8)
                                    base.append(row['same_pulse_acf_error'][lag10].real)
                                    independent.append(row['independent_acf_error'][lag10].real)
                                    truth.append(row['diffuse_acf'][lag10].real)
                                base_means.append(np.mean(base))
                                independent_means.append(np.mean(independent))
                                truth_means.append(np.mean(truth))
                            rng = np.random.default_rng(72123)
                            group['same_pulse_lag10_relative_error_interval'] = summary.paired_ratio(np.asarray(base_means), np.asarray(truth_means), rng)
                            group['independent_lag10_relative_error_interval'] = summary.paired_ratio(np.asarray(independent_means), np.asarray(truth_means), rng)
                            report.flush()
                            print(group.name, 'same-pulse', group['same_pulse_lag10_relative_error_interval'][:].round(4).tolist(), 'independent', group['independent_lag10_relative_error_interval'][:].round(4).tolist(), flush=True)
            report.attrs['complete'] = True
    finally:
        native.close()


if __name__ == '__main__':
    main()
