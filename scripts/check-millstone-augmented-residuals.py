#!/usr/bin/env python3
"""Check native residual-candidate joint refits without writing RF voltages."""
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
    p.add_argument('--augmented', type=pathlib.Path, required=True)
    p.add_argument('--products', type=pathlib.Path, required=True)
    p.add_argument('--detector-script', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    reconstruction = load('residual_reconstruction', pathlib.Path(__file__).with_name('check-millstone-residual-cfar.py'))
    d = load('augmentation_detector', a.detector_script)
    fft = d.FFTBank(128)
    readers = {}
    counts = dict(tested=0, augmentation_errors=0, native_accepted=0, residual_layer=0, accepted_with_residual=0)

    class MemoryRF:
        def __init__(self, y, tx):
            self.values = {'zenith-l': y, 'tx-h': tx}

        def read_vector(self, sample, length, channel):
            return self.values[channel][:length]

    with h5py.File(a.augmented) as source, h5py.File(a.output, 'x') as report:
        assert source.attrs['complete']
        report.attrs['complete'] = False
        report.attrs['scope'] = 'Targeted residual-seed experiment only. No production seeds changed. A smaller coherent residual does not validate removal of diffuse ISR.'
        for name, path in [('augmentation', a.augmented), ('detector', a.detector_script), ('script', pathlib.Path(__file__))]:
            report.attrs[name + '_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        for name, g in source.items():
            if not isinstance(g, h5py.Group):
                continue
            counts['tested'] += 1
            if 'augmentation_error_json' in g:
                counts['augmentation_errors'] += 1
                continue
            sample = int(name)
            c = json.loads(bytes(g['configuration_json'][:]))
            assert c['mode'] in ['alternating-code-16', 'long-pulse-480us']
            assert c['lpi']['sample_rate_hz'] == d.FS == 1e6
            qa = json.loads(bytes(g['augmented_diagnostic_json'][:]))
            family = 'coded' if c['mode'].startswith('alternating') else 'long'
            product = a.products / f'{sample // 10000000 * 10}-{family}.h5'
            with h5py.File(product) as h:
                records = [r for r in h['records'].values() if sample in r['pulse_sample'][:]]
                assert len(records) == 1
                index = int(np.flatnonzero(records[0]['pulse_sample'][:] == sample)[0])
            roots = tuple(c['instrument_metadata']['digital_rf_roots'])
            if roots not in readers:
                readers[roots] = DigitalRFReader([str(pathlib.Path(r)/'rf_data') for r in roots], rdcc_nbytes=64*1024*1024)
            rf = readers[roots]
            raw = np.asarray(rf.read_vector(sample, 10000, 'zenith-l')).reshape(-1)
            tx = np.asarray(rf.read_vector(sample, 10000, 'tx-h')).reshape(-1)
            assert len(raw) == len(tx) == 10000
            clean, rms = reconstruction.reconstruct(raw, c, qa)
            settings = types.SimpleNamespace(decimation_factor=8, fft_padding=4., max_doppler_hz=60000., pfa=1e-10, max_echoes_per_pulse=32,
                receiver_delay_samples=c['instrument_metadata']['tx_to_receiver_delay_samples'][index])
            found = d.detect_pulse(MemoryRF(clean, tx), sample, 1 if family == 'coded' else 300, 'zenith-l', settings, fft)
            layer = [v for v in found if 525 <= v['range_km'] <= 575]
            counts['native_accepted'] += qa['accepted']
            counts['residual_layer'] += bool(layer)
            counts['accepted_with_residual'] += qa['accepted'] and bool(layer)
            row = report.create_group(name)
            row.attrs['native_accepted'] = qa['accepted']
            row.attrs['residual_layer_count'] = len(layer)
            row.attrs['native_rms_reconstruction_error'] = rms - qa['residual_rms_over_noise']
            row['residual_detections_json'] = np.frombuffer(json.dumps(found).encode(), np.uint8)
            report.flush()
        for k, v in counts.items():
            report.attrs[k] = v
        report.attrs['complete'] = True
    print(json.dumps(counts), flush=True)


if __name__ == '__main__':
    main()
