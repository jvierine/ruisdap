#!/usr/bin/env python3
"""CFAR -> native amplitude fit -> native held-pulse forecast controls.

Unlike the oracle control, no injected geometry or phase is given to the
forecaster. All voltages remain in memory/pipes; HDF5 stores paired lag moments
and diagnostics. This is not yet a full lag-profile inversion or campaign test.
"""
import argparse
import hashlib
import importlib.util
import json
import pathlib
import sys
import types

import h5py
import numpy as np


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--canceller', type=pathlib.Path, required=True)
    p.add_argument('--forecaster', type=pathlib.Path, required=True)
    p.add_argument('--detector-script', type=pathlib.Path, required=True)
    p.add_argument('--spectrum', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--groups', type=int, default=4)
    p.add_argument('--pulses', type=int, default=32)
    a = p.parse_args()
    assert a.groups >= 2 and 24 <= a.pulses <= 128 and not a.output.exists()
    folder = pathlib.Path(__file__).parent
    control = load('predictive_cancel', folder/'check-millstone-residual-controls.py')
    native = load('predictive_track', folder/'check-millstone-native-track.py')
    waveform = load('predictive_waveform', folder/'check-millstone-residual-cfar.py')
    summary = load('predictive_summary', folder/'summarize-millstone-residual-controls.py')
    d = load('predictive_detector', a.detector_script)
    with h5py.File(a.spectrum) as h:
        psd = h['ion_line_psd'][:]
        assert h.attrs['sample_rate_hz'] == d.FS == 1e6
    cancellation = control.NativeCancellation(a.canceller.resolve())
    prediction = native.NativePredictor(a.forecaster.resolve())
    fft = d.FFTBank(128)
    parameters = types.SimpleNamespace(decimation_factor=8, fft_padding=4., max_doppler_hz=60000., pfa=1e-10,
                                       max_echoes_per_pulse=32, receiver_delay_samples=0.)
    lags = np.arange(0, 451, 10)

    class MemoryRF:
        def __init__(self, y, tx):
            self.values = {'zenith-l': y, 'tx-h': tx}

        def read_vector(self, sample, length, channel):
            return self.values[channel][:length]

    try:
        with h5py.File(a.output, 'x') as report:
            report.attrs['complete'] = False
            report.attrs['scope'] = __doc__
            report.attrs['selection'] = 'Native training-only model explained_fraction >=0.9. Neighbourhood 0.35 s coded / 0.5 s long. No held-voltage acceptance conditioning or residual-seed addition.'
            report.attrs['intervals'] = 'Paired bootstrap across independent pulse groups; no resampling of dependent pulses within a group'
            for name, path in [('script', pathlib.Path(__file__)), ('canceller', a.canceller), ('forecaster', a.forecaster), ('detector', a.detector_script), ('spectrum', a.spectrum)]:
                report.attrs[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            report['lags_samples'] = lags
            for family, sweep in [('coded', 1), ('long', 300)]:
                mode = d.MODES[sweep]
                body = np.zeros(mode.tx1-mode.tx0, complex)
                code = np.array([1, 1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1])
                body[:480] = np.repeat(code, 30) if family == 'coded' else 1
                tx = np.zeros(10000, complex);tx[mode.tx0:mode.tx1] = body
                keys = (1712587630000000 + (np.r_[0, np.cumsum(np.resize([8910, 17820], a.pulses-1))] if family == 'coded' else np.arange(a.pulses)*26730)).astype(np.int64)
                noise_gate = [7700, 8200] if family == 'coded' else [7200, 7700]
                for width in [0, 240]:
                    offsets = np.array([0]) if width == 0 else np.round(np.linspace(-3*width, 3*width, 49)).astype(int)
                    weights = np.ones(1) if width == 0 else np.exp(-.5*(offsets/width)**2);weights/=weights.sum()
                    for kind in ['null', 'single']:
                        case = report.create_group(f'{family}/width_{width}/{kind}')
                        case.attrs['nominal_diffuse_receiver_power'] = .2
                        case.attrs['groups'] = a.groups;case.attrs['pulses_per_group'] = a.pulses
                        averages = []
                        total_models = total_false = 0
                        for group_index in range(a.groups):
                            group = case.create_group(str(group_index))
                            frames, events = [], []
                            for i, key in enumerate(keys):
                                rng = np.random.default_rng(1+group_index*a.pulses+i)
                                t = (int(key)-int(keys[0]))/1e6
                                # The true trajectory is used only to generate
                                # voltages and the paired target-free reference.
                                tau = 550e3*2/d.C*1e6-1e6/440.2e6*(1000*t-.5*250*t*t)
                                frequency = 1000-250*t
                                phase = 2*np.pi*(1000*t-.5*250*t*t)+.3
                                target = (3*(1+.1*t)*np.exp(1j*phase)*waveform.template(tx, tau, frequency, 1e6, 16)) if kind == 'single' else np.zeros(len(tx), complex)
                                receiver = (rng.standard_normal(len(tx))+1j*rng.standard_normal(len(tx)))/np.sqrt(2)
                                diffuse = np.zeros(len(tx), complex)
                                profile = np.zeros(len(tx))
                                for offset, weight in zip(offsets, weights):
                                    z = (rng.standard_normal(len(psd))+1j*rng.standard_normal(len(psd)))/np.sqrt(2)
                                    z = np.fft.ifft(z*np.sqrt(psd))*len(psd)/np.sqrt(psd.sum())
                                    start = round(550e3*2/d.C*1e6)+mode.tx0+int(offset)
                                    diffuse[start:start+len(body)] += body*z[np.arange(start,start+len(body)) % len(z)]*np.sqrt(.2*weight)
                                    profile[start:start+len(body)] += weight*abs(body)**2
                                ideal = receiver+diffuse+100*tx;raw = ideal+target
                                measured = np.zeros(len(tx), complex);measured[mode.tx0:mode.tx1] = raw[mode.tx0:mode.tx1]
                                found = d.detect_pulse(MemoryRF(raw, tx), int(key), sweep, 'zenith-l', parameters, fft)
                                seeds = [dict(delay_samples=v['corrected_delay_samples'], doppler_hz=v['doppler_hz']) for v in found]
                                baseline, mask, qa = cancellation.cancel(raw, measured, seeds, noise_gate)
                                assert 'error' not in qa, qa
                                for v, amplitude, snr in zip(qa['detections'], qa['amplitudes'], qa['amplitude_snr_db']):
                                    if 250 <= v['delay_samples']*.149896229 <= 575:
                                        events.append(dict(sample=int(key), delay_samples=v['delay_samples'], doppler_hz=v['doppler_hz'], amplitude=complex(*amplitude), snr_db=snr))
                                frames.append((raw, ideal, diffuse, measured, baseline, qa, profile > profile.max()*1e-6))
                            values = []
                            for key, (raw, ideal, diffuse, measured, baseline, qa, support) in zip(keys, frames):
                                models = prediction.predict(events, int(key), (mode.tx0+mode.tx1)/2, radius_s=.5 if family == 'long' else .35)
                                models = [m for m in models if m['supported'] and m['explained_fraction'] >= .9]
                                assert all(int(key) not in m['training_samples'] for m in models)
                                corrected = raw.copy()
                                for model in models:
                                    corrected -= complex(*model['amplitude'])*waveform.template(measured, model['delay_samples'], model['doppler_hz'], 1e6, 16)
                                total_models += bool(models)
                                total_false += kind == 'null' and bool(models)
                                row = group.create_group(str(int(key)))
                                ref = control.moments(ideal, support, lags)
                                row['baseline_acf_error'] = control.moments(baseline, support, lags)-ref
                                row['predictive_acf_error'] = control.moments(corrected, support, lags)-ref
                                row['diffuse_acf'] = control.moments(diffuse, support, lags)
                                row['models_json'] = np.frombuffer(json.dumps(models).encode(), np.uint8)
                                row['baseline_qa_json'] = np.frombuffer(json.dumps(qa).encode(), np.uint8)
                                values.append([row['baseline_acf_error'][1].real,row['predictive_acf_error'][1].real,row['diffuse_acf'][1].real])
                            averages.append(np.mean(values, axis=0))
                            report.flush()
                        averages = np.asarray(averages);rng = np.random.default_rng(2935)
                        case['baseline_lag10_relative_error_interval'] = summary.paired_ratio(averages[:,0], averages[:,2], rng)
                        case['predictive_lag10_relative_error_interval'] = summary.paired_ratio(averages[:,1], averages[:,2], rng)
                        case.attrs['predicted_subtraction_pulses'] = total_models;case.attrs['false_subtraction_pulses'] = total_false
                        print(case.name, 'predicted', total_models, 'false', total_false, 'baseline', case['baseline_lag10_relative_error_interval'][:].round(4).tolist(), 'predictive', case['predictive_lag10_relative_error_interval'][:].round(4).tolist(), flush=True)
            report.attrs['complete'] = True
    finally:
        cancellation.close();prediction.close()


if __name__ == '__main__':
    main()
