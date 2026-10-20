#!/usr/bin/env python3
"""Compare native held-pulse forecasts with the independent Python reference.

Only amplitude metadata crosses an anonymous pipe; no RF is opened or stored.
The output also retains native forecasts for subsequent original-RF checks.
"""
import argparse
import hashlib
import importlib.util
import json
import pathlib
import subprocess
import sys
import time

import h5py
import numpy as np


class NativePredictor:
    def __init__(self, binary):
        self.process = subprocess.Popen([str(binary)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)

    def predict(self, events, held, centre, radius_s=None):
        values = [dict(e, amplitude=[e['amplitude'].real, e['amplitude'].imag]) for e in events]
        request = dict(events=values, held_sample=int(held), tx_centre=float(centre))
        if radius_s is not None:
            request['settings'] = dict(radius_s=radius_s)
        self.process.stdin.write(json.dumps(request)+'\n')
        self.process.stdin.flush()
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError('Native track predictor ended without a response')
        return json.loads(line)

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)
        self.process.stdout.close()
        assert self.process.returncode == 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=pathlib.Path, required=True)
    p.add_argument('--reference', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    path = pathlib.Path(__file__).with_name('check-millstone-track-prediction.py')
    spec = importlib.util.spec_from_file_location('native_track_reference', path)
    reference = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = reference
    spec.loader.exec_module(reference)
    native = NativePredictor(a.binary.resolve())
    total, matched, worst_amplitude = 0, 0, 0.
    try:
        with h5py.File(a.reference) as source, h5py.File(a.output, 'x') as report:
            assert source.attrs['complete']
            report.attrs['complete'] = False
            report.attrs['scope'] = __doc__
            report.attrs['tolerances'] = 'Training membership and support must match; geometry <=1e-6 sample/Hz; complex amplitude relative difference <=0.005; native explained fraction no worse by >1e-6'
            for name, path in [('binary', a.binary), ('reference', a.reference), ('script', pathlib.Path(__file__))]:
                report.attrs[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
            for product_name, product in source.items():
                for record_id, record in product.items():
                    original = pathlib.Path(record.attrs['source_path'])
                    assert hashlib.sha256(original.read_bytes()).hexdigest() == record.attrs['source_sha256']
                    with h5py.File(original) as h:
                        c = json.loads(bytes(h['records'][record_id]['configuration_json'][:]))
                        _, _, events = reference.events_from_record(h['records'][record_id], c)
                    out = report.create_group(product_name+'/'+record_id)
                    out.attrs.update(record.attrs)
                    out['configuration_json'] = record['configuration_json'][:]
                    for region_name, region in record.items():
                        if not isinstance(region, h5py.Group):
                            continue
                        lo, hi = map(float, region_name.split('_')[1:])
                        rows = [e for e in events if lo <= e['delay_samples']*.149896229 <= hi]
                        group = out.create_group(region_name)
                        for sample, row in region.items():
                            before = time.monotonic()
                            predictions = native.predict(rows, int(sample), np.mean(c['lpi']['tx_gate']), radius_s=float(record.attrs.get('prediction_radius_s', .35)))
                            elapsed = time.monotonic()-before
                            py = json.loads(bytes(row['predictions_json'][:]))
                            equal = len(py) == len(predictions)
                            metrics = []
                            for x, y in zip(predictions, py):
                                equal &= x['supported'] == y['supported']
                                if not x['supported'] or not y['supported']:
                                    equal &= x['reason'] == y['reason']
                                    continue
                                assert int(sample) not in x['training_samples']
                                equal &= x['training_samples'] == y['training_samples']
                                relative = abs(complex(*x['amplitude'])-complex(*y['amplitude']))/max(1e-12, abs(complex(*y['amplitude'])))
                                delta = [abs(x['delay_samples']-y['delay_samples']), abs(x['doppler_hz']-y['doppler_hz']), relative, x['explained_fraction']-y['explained_fraction']]
                                metrics.append(delta)
                                worst_amplitude = max(worst_amplitude, relative)
                                equal &= delta[0] <= 1e-6 and delta[1] <= 1e-6 and delta[2] <= .005 and delta[3] >= -1e-6
                            result = group.create_group(sample)
                            result.attrs['equivalent_within_tolerance'] = equal
                            result.attrs['native_seconds'] = elapsed
                            result['comparison'] = np.asarray(metrics).reshape(-1, 4)
                            result['predictions_json'] = np.frombuffer(json.dumps(predictions).encode(), np.uint8)
                            result['held_diagnostic_json'] = row['held_diagnostic_json'][:]
                            total += 1
                            matched += bool(equal)
                            report.flush()
            report.attrs['complete'] = True
            report.attrs['compared_queries'] = total
            report.attrs['matched_queries'] = matched
            report.attrs['worst_relative_amplitude_difference'] = worst_amplitude
            report.attrs['equivalent_within_tolerance'] = total == matched
    finally:
        native.close()
    print(dict(compared=total, matched=matched, worst_relative_amplitude_difference=worst_amplitude), flush=True)
    if total != matched:
        raise SystemExit('Native/reference forecast differences require investigation')


if __name__ == '__main__':
    main()
