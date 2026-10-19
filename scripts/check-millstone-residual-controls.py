#!/usr/bin/env python3
"""Synthetic CFAR/native-cancellation controls for one residual-seed joint refit.

All synthetic voltages remain in RAM and anonymous process pipes. HDF5 reports
contain diagnostics and paired ACF moments only. This tests a stated stochastic
model, not the empirical campaign false-alarm rate or the full LPI inversion.
"""
import argparse
import hashlib
import importlib.util
import json
import pathlib
import struct
import subprocess
import sys
import types

import h5py
import numpy as np


class NativeCancellation:
    def __init__(self, binary):
        self.process = subprocess.Popen([str(binary)], stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    def read(self, n):
        data = bytearray()
        while len(data) < n:
            part = self.process.stdout.read(n-len(data))
            if not part:
                raise RuntimeError('Native diagnostic helper ended inside a response')
            data.extend(part)
        return bytes(data)

    def cancel(self, y, tx, detections, noise_gate, settings=None):
        header = dict(n=len(y), sample_rate_hz=1e6, noise_gate=noise_gate, detections=detections)
        if settings is not None:
            header['settings'] = settings
        request = json.dumps(header).encode()
        self.process.stdin.write(struct.pack('<Q', len(request)))
        self.process.stdin.write(request)
        for array in [y, tx]:
            self.process.stdin.write(np.asarray(array, '<c16').tobytes())
        self.process.stdin.flush()
        size, = struct.unpack('<Q', self.read(8))
        assert 0 < size <= 1_000_000
        qa = json.loads(self.read(size))
        if 'error' in qa:
            return None, None, qa
        clean = np.frombuffer(self.read(len(y)*16), '<c16').copy()
        mask = np.frombuffer(self.read(len(y)), np.uint8).astype(bool)
        return clean, mask, qa

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)
        self.process.stdout.close()
        assert self.process.returncode == 0


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def moments(y, mask, lags):
    values = []
    for lag in lags:
        use = mask[:len(y)-lag] & mask[lag:]
        values.append(np.mean(y[:len(y)-lag][use] * np.conj(y[lag:][use])) if use.any() else complex(np.nan, np.nan))
    return np.asarray(values)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=pathlib.Path, required=True)
    p.add_argument('--detector-script', type=pathlib.Path, required=True)
    p.add_argument('--spectrum', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--trials', type=int, default=50)
    a = p.parse_args()
    assert a.trials > 0 and not a.output.exists()
    d = load('control_detector', a.detector_script)
    interpolation = load('control_interpolation', pathlib.Path(__file__).with_name('check-millstone-residual-cfar.py'))
    with h5py.File(a.spectrum) as h:
        psd = h['ion_line_psd'][:]
        assert h.attrs['sample_rate_hz'] == d.FS == 1e6
    native = NativeCancellation(a.binary.resolve())
    fft = d.FFTBank(128)
    settings = types.SimpleNamespace(decimation_factor=8, fft_padding=4., max_doppler_hz=60000., pfa=1e-10,
                                     max_echoes_per_pulse=32, receiver_delay_samples=0.)

    class MemoryRF:
        def __init__(self, y, tx):
            self.values = {'zenith-l': y, 'tx-h': tx}

        def read_vector(self, sample, length, channel):
            return self.values[channel][:length]

    def detect(y, tx, sweep):
        return d.detect_pulse(MemoryRF(y, tx), 0, sweep, 'zenith-l', settings, fft)

    def seed(v):
        return dict(delay_samples=v['corrected_delay_samples'], doppler_hz=v['doppler_hz'])

    lags = np.arange(0, 451, 10)
    try:
        with h5py.File(a.output, 'x') as report:
            report.attrs['complete'] = False
            report.attrs['scope'] = __doc__
            report.attrs['algorithm'] = 'Unchanged CFAR; native cancellation; if accepted/applied/confirmed, detect residual and add nonduplicate seeds (>2 samples or >300 Hz from existing component); one joint refit from original voltages; no further iteration.'
            for name, path in [('native_binary', a.binary), ('detector', a.detector_script), ('spectrum', a.spectrum), ('script', pathlib.Path(__file__)), ('interpolation', pathlib.Path(interpolation.__file__))]:
                report.attrs[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            report['lags_samples'] = lags
            for family, sweep in [('coded', 1), ('long', 300)]:
                mode = d.MODES[sweep]
                n = 10000
                delay = round(550e3*2/d.C*1e6)
                code = np.array([1, 1, -1, 1, -1, -1, 1, -1, 1, 1, 1, -1, -1, 1, -1, -1])
                body = np.zeros(mode.tx1-mode.tx0, complex)
                body[:480] = 1 if family == 'long' else np.repeat(code, 30)
                tx = np.zeros(n, complex)
                tx[mode.tx0:mode.tx1] = body
                noise_gate = [7700, 8200] if family == 'coded' else [7200, 7700]
                for width in [0, 240]:
                    offsets = np.array([0]) if width == 0 else np.round(np.linspace(-3*width, 3*width, 49)).astype(int)
                    weights = np.ones(1) if width == 0 else np.exp(-.5*(offsets/width)**2)
                    weights /= weights.sum()
                    profile = np.zeros(n)
                    for offset, weight in zip(offsets, weights):
                        start = delay+mode.tx0+int(offset)
                        profile[start:start+len(body)] += weight*abs(body)**2
                    support = profile > profile.max()*1e-6
                    cases = [('null', power) for power in [0., .03, .1, .2, .3, 1.]] + [(kind, power) for kind in ['single', 'double'] for power in [.03, .1, .2]]
                    for kind, power in cases:
                        group = report.create_group(f'{family}/width_{width}/{kind}_{power}')
                        group.attrs['nominal_diffuse_receiver_power'] = power
                        group.attrs['maximum_expected_diffuse_receiver_ratio'] = power*profile.max()
                        counts = dict(trials=0, kernel_errors=0, initial_layer_candidates=0, initial_accepted_subtractions=0, joint_refits=0, final_accepted_subtractions=0, final_masked=0, residual_layer=0)
                        target_specs = [] if kind == 'null' else [(delay, 1000., 3*np.exp(.3j))]
                        if kind == 'double':
                            target_specs.append((delay+60.25, 2500., .6*np.exp(-.7j)))
                        target = np.zeros(n, complex)
                        for tau, frequency, amplitude in target_specs:
                            target += amplitude*interpolation.template(tx, tau, frequency, 1e6, 16)
                        group['target_parameters_json'] = np.frombuffer(json.dumps([[tau, frequency, amplitude.real, amplitude.imag] for tau, frequency, amplitude in target_specs]).encode(), np.uint8)
                        for trial in range(a.trials):
                            rng = np.random.default_rng(trial+1)
                            receiver = (rng.standard_normal(n)+1j*rng.standard_normal(n))/np.sqrt(2)
                            diffuse = np.zeros(n, complex)
                            if power:
                                for offset, weight in zip(offsets, weights):
                                    z = (rng.standard_normal(len(psd))+1j*rng.standard_normal(len(psd)))/np.sqrt(2)
                                    z = np.fft.ifft(z*np.sqrt(psd))*len(psd)/np.sqrt(psd.sum())
                                    start = delay+mode.tx0+int(offset)
                                    diffuse[start:start+len(body)] += body*z[np.arange(start, start+len(body)) % len(z)]*np.sqrt(power*weight)
                            ideal = receiver+diffuse+100*tx
                            raw = ideal+target
                            measured = np.zeros(n, complex)
                            measured[mode.tx0:mode.tx1] = raw[mode.tx0:mode.tx1]
                            initial = detect(raw, tx, sweep)
                            seeds = [seed(v) for v in initial]
                            clean, usable, qa = native.cancel(raw, measured, seeds, noise_gate)
                            row = group.create_group(str(trial+1))
                            counts['trials'] += 1
                            counts['initial_layer_candidates'] += any(525 <= v['range_km'] <= 575 for v in initial)
                            if 'error' in qa:
                                counts['kernel_errors'] += 1
                                row['kernel_error_json'] = np.frombuffer(json.dumps(qa).encode(), np.uint8)
                                continue
                            base, base_mask, base_qa = clean, usable, qa
                            counts['initial_accepted_subtractions'] += qa['accepted'] and qa['applied']
                            residual = detect(clean, tx, sweep)
                            extra = []
                            if qa['accepted'] and qa['applied'] and qa['confirmed']:
                                for v in residual:
                                    candidate = seed(v)
                                    if all(abs(candidate['delay_samples']-old['delay_samples']) > 2 or abs(candidate['doppler_hz']-old['doppler_hz']) > 300 for old in qa['detections']):
                                        extra.append(candidate)
                            if extra:
                                counts['joint_refits'] += 1
                                clean, usable, qa = native.cancel(raw, measured, seeds+extra, noise_gate)
                                if 'error' in qa:
                                    counts['kernel_errors'] += 1
                                    row['kernel_error_json'] = np.frombuffer(json.dumps(qa).encode(), np.uint8)
                                    continue
                            final = detect(clean, tx, sweep) if extra else residual
                            counts['final_accepted_subtractions'] += qa['accepted'] and qa['applied']
                            counts['final_masked'] += not qa['accepted']
                            counts['residual_layer'] += any(525 <= v['range_km'] <= 575 for v in final)
                            # Paired moments use identical availability masks for
                            # each estimate and its ideal target-free reference.
                            for label, y, mask in [('baseline', base, base_mask), ('joint', clean, usable)]:
                                valid = support & mask
                                row[label+'_acf_error'] = moments(y, valid, lags)-moments(ideal, valid, lags)
                                row.attrs[label+'_masked_fraction'] = 1-valid.sum()/support.sum()
                                row.attrs[label+'_distortion_power'] = np.mean(abs(y[support]-ideal[support])**2)
                            row['diffuse_acf'] = moments(diffuse, support, lags)
                            row['raw_acf_error'] = moments(raw, support, lags)-moments(ideal, support, lags)
                            for label, value in [('initial', initial), ('residual', residual), ('extra_seeds', extra), ('final', final), ('baseline_qa', base_qa), ('joint_qa', qa)]:
                                row[label+'_json'] = np.frombuffer(json.dumps(value).encode(), np.uint8)
                        for name, value in counts.items():
                            group.attrs[name] = value
                        report.flush()
                        print(f'{family} width={width} {kind} power={power}: {json.dumps(counts)}', flush=True)
            report.attrs['complete'] = True
    finally:
        native.close()


if __name__ == '__main__':
    main()
