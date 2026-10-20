#!/usr/bin/env python3
"""Test training-only forecasts on original held RF, without storing voltages.

Forecasts are completed before opening held-pulse RF. A fixed training-only
explained-fraction threshold selects models; no coefficient is fitted to the
queried voltage. CFAR residuals are diagnostics, not target identifications.
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
from digital_rf import DigitalRFReader


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--forecasts', type=pathlib.Path, required=True)
    p.add_argument('--detector-script', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    reconstruct = load('forecast_voltage_reconstruction', pathlib.Path(__file__).with_name('check-millstone-residual-cfar.py'))
    control = load('forecast_voltage_moments', pathlib.Path(__file__).with_name('check-millstone-residual-controls.py'))
    d = load('forecast_voltage_detector', a.detector_script)
    fft = d.FFTBank(128)
    readers = {}
    lags = np.arange(0, 451, 10)

    class MemoryRF:
        def __init__(self, y, tx):
            self.values = {'zenith-l': y, 'tx-h': tx}

        def read_vector(self, sample, length, channel):
            return self.values[channel][:length]

    with h5py.File(a.forecasts) as source, h5py.File(a.output, 'x') as report:
        assert source.attrs['complete']
        report.attrs['complete'] = False
        report.attrs['scope'] = __doc__
        report.attrs['model_selection'] = 'supported training-only forecasts with explained_fraction >=0.9; no held-voltage acceptance conditioning, amplitude refinement or residual-seed addition'
        for name, path in [('forecasts', a.forecasts), ('detector', a.detector_script), ('script', pathlib.Path(__file__)), ('reconstruction', pathlib.Path(reconstruct.__file__)), ('moments', pathlib.Path(control.__file__))]:
            report.attrs[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        report['lag_samples'] = lags
        for product_name, product in source.items():
            for record_id, record in product.items():
                original = pathlib.Path(record.attrs['source_path'])
                assert hashlib.sha256(original.read_bytes()).hexdigest() == record.attrs['source_sha256']
                c = json.loads(bytes(record['configuration_json'][:]))
                assert c['mode'] in ['alternating-code-16', 'long-pulse-480us']
                assert c['lpi']['sample_rate_hz'] == d.FS == 1e6
                by_pulse = {}
                for region in record.values():
                    if not isinstance(region, h5py.Group):
                        continue
                    for sample, g in region.items():
                        models = json.loads(bytes(g['predictions_json'][:]))
                        qa = json.loads(bytes(g['held_diagnostic_json'][:]))
                        entry = by_pulse.setdefault(int(sample), dict(models=[], qa=qa))
                        entry['models'].extend({k: value for k, value in v.items() if k != 'assessment'}
                                               for v in models if v['supported'] and v['explained_fraction'] >= .9)
                roots = tuple(c['instrument_metadata']['digital_rf_roots'])
                if roots not in readers:
                    readers[roots] = DigitalRFReader([str(pathlib.Path(r)/'rf_data') for r in roots], rdcc_nbytes=64*1024*1024)
                rf = readers[roots]
                out = report.create_group(product_name+'/'+record_id)
                out.attrs['source_product_sha256'] = record.attrs['source_sha256']
                with h5py.File(original) as h:
                    keys = h['records'][record_id]['pulse_sample'][:]
                count, applied, base_hits, independent_hits, empty_new_hits = 0, 0, 0, 0, 0
                applied_base_hits = applied_independent_hits = 0
                for sample, entry in sorted(by_pulse.items()):
                    models, qa = entry['models'], entry['qa']
                    assert all(sample not in m['training_samples'] for m in models)
                    if qa['reason'] in ['invalid_background_noise', 'background_gate_contaminated']:
                        row = out.create_group(str(sample))
                        row.attrs['not_evaluable_reason'] = qa['reason']
                        continue
                    index = int(np.flatnonzero(keys == sample)[0])
                    raw = np.asarray(rf.read_vector(sample, 10000, 'zenith-l')).reshape(-1)
                    tx_h = np.asarray(rf.read_vector(sample, 10000, 'tx-h')).reshape(-1)
                    assert len(raw) == len(tx_h) == 10000
                    y = np.asarray(raw, complex)-complex(*c['instrument_metadata']['dc_offset'])
                    baseline, rms = reconstruct.reconstruct(raw, c, qa)
                    measured = np.zeros(len(y), complex)
                    lo, hi = c['lpi']['tx_gate']
                    measured[lo:hi] = y[lo:hi]
                    predicted = y.copy()
                    support = np.zeros(len(y), bool)
                    for model in models:
                        t = reconstruct.template(measured, model['delay_samples'], model['doppler_hz'], 1e6, c['cancellation']['interpolation_half'])
                        support |= abs(t) > abs(t).max()*1e-3
                        predicted -= complex(*model['amplitude'])*t
                    settings = types.SimpleNamespace(decimation_factor=8, fft_padding=4., max_doppler_hz=60000., pfa=1e-10, max_echoes_per_pulse=32,
                        receiver_delay_samples=c['instrument_metadata']['tx_to_receiver_delay_samples'][index])
                    sweep = 1 if c['mode'].startswith('alternating') else 300
                    before = d.detect_pulse(MemoryRF(raw, tx_h), sample, sweep, 'zenith-l', settings, fft)
                    base = d.detect_pulse(MemoryRF(baseline, tx_h), sample, sweep, 'zenith-l', settings, fft)
                    after = d.detect_pulse(MemoryRF(predicted, tx_h), sample, sweep, 'zenith-l', settings, fft)
                    relevant = lambda values: [v for v in values if 250 <= v['range_km'] <= 575]
                    count += 1
                    applied += bool(models)
                    base_hits += bool(relevant(base))
                    independent_hits += bool(relevant(after))
                    empty_new_hits += not bool(relevant(before)) and bool(relevant(after))
                    if models:
                        applied_base_hits += bool(relevant(base))
                        applied_independent_hits += bool(relevant(after))
                    row = out.create_group(str(sample))
                    row.attrs['applied_model_count'] = len(models)
                    row.attrs['native_accepted'] = qa['accepted']
                    row.attrs['native_rms_reconstruction_error'] = rms-qa['residual_rms_over_noise'] if qa['detections'] else 0.
                    for name, values in [('before', before), ('baseline', base), ('independent', after), ('training_models', models)]:
                        row[name+'_json'] = np.frombuffer(json.dumps(values).encode(), np.uint8)
                    if support.any():
                        for name, values in [('raw', y), ('baseline', baseline), ('independent', predicted)]:
                            row[name+'_support_acf'] = control.moments(values, support, lags)
                    report.flush()
                print(product_name, dict(tested=count, predicted_subtraction=applied, baseline_residual=base_hits,
                      independent_residual=independent_hits, applied_baseline_residual=applied_base_hits,
                      applied_independent_residual=applied_independent_hits, newly_detectable_after_independent=empty_new_hits), flush=True)
        report.attrs['complete'] = True


if __name__ == '__main__':
    main()
