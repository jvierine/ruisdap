#!/usr/bin/env python3
"""Re-run the unchanged CFAR detector before/after recorded native subtraction.

Read original DigitalRF directly. Reconstruct the exact recorded subtraction in
RAM, checking its support RMS against native diagnostics before any inference.
Only detections and scalar QA enter the HDF5 report; no voltages are written.
This is a targeted residual check, not a campaign false-alarm measurement.
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


def template(tx, delay, doppler, sample_rate, half):
    """Noncircular Lanczos interpolation, in native contributing-sum order."""
    t = np.arange(len(tx)) - delay
    centre = np.floor(t).astype(np.int64)
    z = np.zeros(len(tx), np.complex128)
    for offset in range(-half + 1, half + 1):
        j = centre + offset
        u = t - j
        use = (j >= 0) & (j < len(tx)) & (abs(u) < half)
        z[use] += tx[j[use]] * np.sinc(u[use]) * np.sinc(u[use] / half)
    return z * np.exp(2j * np.pi * doppler * np.arange(len(tx)) / sample_rate)


def reconstruct(raw, configuration, diagnostic):
    assert len(diagnostic['detections']) == len(diagnostic['amplitudes'])
    dc = complex(*configuration['instrument_metadata']['dc_offset'])
    y = np.asarray(raw, np.complex128) - dc
    tx = np.zeros(len(y), np.complex128)
    lo, hi = configuration['lpi']['tx_gate']
    tx[lo:hi] = y[lo:hi]
    cleaned = y.copy()
    support = np.zeros(len(y), bool)
    for detection, amplitude in zip(diagnostic['detections'], diagnostic['amplitudes']):
        v = template(tx, detection['delay_samples'], detection['doppler_hz'],
                     configuration['lpi']['sample_rate_hz'],
                     configuration['cancellation']['interpolation_half'])
        support |= abs(v) > abs(v).max() * 1e-3
        cleaned -= complex(*amplitude) * v
    if diagnostic['detections']:
        rms = np.sqrt(np.mean(abs(cleaned[support]) ** 2)) / diagnostic['noise_rms']
        # Reject an incorrect reconstruction instead of studying a surrogate.
        np.testing.assert_allclose(rms, diagnostic['residual_rms_over_noise'], rtol=1e-9, atol=1e-10)
    else:
        rms = np.nan
        np.testing.assert_array_equal(cleaned, y)
    return cleaned, rms


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--detector-script', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--product', type=pathlib.Path, action='append', required=True)
    p.add_argument('--limit', type=int, default=32, help='evenly spaced eligible pulses per product')
    p.add_argument('--controls', type=int, default=4, help='no-seed science pulses per product')
    p.add_argument('--failed-limit', type=int, default=16, help='failed-validation layer pulses per product')
    a = p.parse_args()
    assert a.limit > 0 and a.controls >= 0 and a.failed_limit >= 0 and not a.output.exists()
    spec = importlib.util.spec_from_file_location('residual_cfar_reference', a.detector_script)
    d = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = d
    spec.loader.exec_module(d)
    from digital_rf import DigitalRFReader

    class MemoryRF:
        def __init__(self, y, tx):
            self.arrays = {'zenith-l': y, 'tx-h': tx}

        def read_vector(self, sample, length, channel):
            return self.arrays[channel][:length]

    fft = d.FFTBank(128)
    readers = {}
    with h5py.File(a.output, 'x') as report:
        report.attrs['complete'] = False
        report.attrs['detector_sha256'] = hashlib.sha256(a.detector_script.read_bytes()).hexdigest()
        report.attrs['script_sha256'] = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
        report.attrs['scope'] = 'Targeted 525-575 km accepted/applied/confirmed and failed-validation science pulses, evenly spaced within each category and product; no-seed controls. Not a campaign false-alarm estimate.'
        for product in a.product:
            with h5py.File(product) as h:
                for record_id, g in h['records'].items():
                    c = json.loads(bytes(g['configuration_json'][:]))
                    qa = json.loads(bytes(g['debris/diagnostics_json'][:]))
                    keys = g['pulse_sample'][:]
                    assert len(keys) == len(qa) == len(c['lpi']['science_pulse'])
                    eligible, controls, failed = [], [], []
                    for i, q in enumerate(qa):
                        if not c['lpi']['science_pulse'][i]:
                            continue
                        if not q['detections']:
                            controls.append(i)
                        elif q['amplitudes'] and any(525 <= v['delay_samples'] * d.C / 2e9 <= 575 for v in q['detections']):
                            if q['accepted'] and q['applied'] and q['confirmed']:
                                eligible.append(i)
                            elif not q['accepted']:
                                failed.append(i)
                    select = lambda v, n: [v[i] for i in np.unique(np.linspace(0, len(v)-1, min(n, len(v))).astype(int))] if v and n else []
                    selected = ([(i, 'layer') for i in select(eligible, a.limit)]
                                + [(i, 'no_seed') for i in select(controls, a.controls)]
                                + [(i, 'masked') for i in select(failed, a.failed_limit)])
                    roots = tuple(c['instrument_metadata']['digital_rf_roots'])
                    if roots not in readers:
                        readers[roots] = DigitalRFReader([str(pathlib.Path(r) / 'rf_data') for r in roots],
                            rdcc_nbytes=64 * 1024 * 1024)
                    rf = readers[roots]
                    out = report.create_group(product.stem + '/' + record_id)
                    out.attrs['product_path'] = str(product)
                    out.attrs['product_sha256'] = hashlib.sha256(product.read_bytes()).hexdigest()
                    out.attrs['eligible_layer_pulses'] = len(eligible)
                    out.attrs['eligible_masked_layer_pulses'] = len(failed)
                    out['configuration_json'] = np.frombuffer(json.dumps(c).encode(), np.uint8)
                    n_residual, n_other = 0, 0
                    for i, kind in selected:
                        sample = int(keys[i])
                        raw = np.asarray(rf.read_vector(sample, 10000, 'zenith-l')).reshape(-1)
                        tx = np.asarray(rf.read_vector(sample, 10000, 'tx-h')).reshape(-1)
                        assert len(raw) == len(tx) == 10000
                        cleaned, rms = reconstruct(raw, c, qa[i])
                        settings = types.SimpleNamespace(decimation_factor=8, fft_padding=4., max_doppler_hz=60000., pfa=1e-10, max_echoes_per_pulse=32,
                            receiver_delay_samples=c['instrument_metadata']['tx_to_receiver_delay_samples'][i])
                        sweep = 1 if c['mode'].startswith('alternating') else 300
                        before = d.detect_pulse(MemoryRF(raw, tx), sample, sweep, 'zenith-l', settings, fft)
                        after = d.detect_pulse(MemoryRF(cleaned, tx), sample, sweep, 'zenith-l', settings, fft)
                        layer = [v for v in after if 525 <= v['range_km'] <= 575]
                        other = [v for v in layer if all(abs(v['doppler_hz'] - old['doppler_hz']) > 4000 for old in qa[i]['detections'])]
                        n_residual += bool(layer)
                        n_other += bool(other)
                        row = out.create_group(str(sample))
                        row.attrs['selection'] = kind
                        row.attrs['native_rms_reconstruction_error'] = rms - qa[i]['residual_rms_over_noise'] if kind != 'no_seed' else 0.
                        row.attrs['residual_layer_count'] = len(layer)
                        row.attrs['residual_other_doppler_count'] = len(other)
                        for name, value in [('native_diagnostic', qa[i]), ('before', before), ('after', after)]:
                            row[name + '_json'] = np.frombuffer(json.dumps(value).encode(), np.uint8)
                        report.flush()
                    print(product.name, 'eligible', len(eligible), 'tested', len(selected), 'residual_layer', n_residual, 'different_doppler', n_other, flush=True)
        report.attrs['complete'] = True


if __name__ == '__main__':
    main()
