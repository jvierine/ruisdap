#!/usr/bin/env python3
"""Fit original-RF native LPI injection controls with the frozen campaign fitter.

Only reduced signed ACFs and derived fit diagnostics are saved. This is a
ten-second matched control, not a campaign result or an empirical CFAR test.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import pathlib
import shutil
import subprocess

import h5py
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def data_json(group, name, value):
    group[name] = np.frombuffer(json.dumps(value).encode(), np.uint8)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', type=pathlib.Path, required=True)
    p.add_argument('--controls', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    with h5py.File(a.campaign/'kernel.h5') as ledger:
        binary = a.campaign/str(ledger.attrs['fit_binary_name'])
        assert digest(binary) == ledger.attrs['fit_kernel_sha256']
    day = a.campaign/'2024-04-08'
    template = json.loads((day/'stream-jobs-legacy300-v3/1712599200/fit.json').read_text())
    jobs = []
    source_ledger = []
    with h5py.File(a.output/'verification.h5', 'x') as verification:
        verification.attrs['complete'] = False
        verification.attrs['scope'] = __doc__
        verification.attrs['fit_kernel_sha256'] = digest(binary)
        for target in [300, 550]:
            report_path = a.controls/f'track-lpi-injection-{target}-v1.h5'
            with h5py.File(report_path) as report:
                assert report.attrs['complete']
                source_ledger.append(dict(path=str(report_path), sha256=digest(report_path)))
                for branch in ['original', 'injected', 'same_pulse', 'independent']:
                    folder = a.output/str(target)/branch
                    folder.mkdir(parents=True)
                    inputs = []
                    for mode in ['coded', 'long']:
                        source = day/'products'/f'1712599210-{mode}.h5'
                        dest = folder/source.name
                        shutil.copy2(source, dest)
                        g = report[mode]
                        with h5py.File(dest, 'r+') as product:
                            record = product['records'][next(iter(product['records']))]
                            for name in ['acf', 'ground_clutter_acf', 'variance']:
                                original = record['ion_line/'+name][:]
                                assert np.array_equal(original, g['original/'+name][:], equal_nan=True), 'Original native LPI differs from the canonical source'
                                record['ion_line/'+name][...] = g[branch+'/'+name][:]
                            np.testing.assert_array_equal(record['pulse_sample'][:], g['pulse_sample'][:])
                            np.testing.assert_allclose(record['range_m'][:]/1000, g['slant_range_km'][:], rtol=1e-14, atol=1e-12)
                            np.testing.assert_allclose(record['ion_line/lag_s'][:]*1e6, g['lag_samples'][:], rtol=1e-14, atol=1e-12)
                            for name in ['receiver_gain', 'system_temperature_k']:
                                np.testing.assert_equal(record['ion_line/'+name][0], g[branch].attrs[name])
                            record['ion_line/retained_fraction'][...] = g[branch+'/retained_fraction'][:]
                            control = record.create_group('injection_control')
                            data_json(control, 'provenance_json', dict(report=str(report_path), report_sha256=digest(report_path),
                                source=str(source), source_sha256=digest(source), branch=branch, target_slant_range_km=target,
                                note='Original RF plus declared synthetic target. Reduced ACF replacement only; never publish as observations.'))
                            for name in ['injection_json', 'independent_forecasts_json', 'prediction_failures_json', 'legacy_cancellation_diagnostics_json']:
                                control[name] = g[name][:]
                        inputs.append(str(dest))
                    plan = dict(template, inputs=inputs, output=str(folder/'fits'), range_km=[200., 650.])
                    path = folder/'fit.json'
                    path.write_text(json.dumps(plan))
                    jobs.append((target, branch, path))
        data_json(verification, 'sources_json', source_ledger)
    env = {**os.environ, 'OMP_NUM_THREADS':'1', 'OPENBLAS_NUM_THREADS':'1'}
    env.pop('DISPLAY', None)
    def fit(job):
        target, branch, plan = job
        with (plan.parent/'fit.log').open('x') as log:
            subprocess.run([str(binary), 'fit', str(plan)], check=True, env=env, stdout=log, stderr=log)
        print('FIT READY', target, branch, flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(fit, jobs))
    with h5py.File(a.output/'verification.h5', 'a') as verification:
        for target in [300, 550]:
            fits = {}
            for branch in ['original', 'injected', 'same_pulse', 'independent']:
                path = a.output/str(target)/branch/'fits/fit-000000.h5'
                with h5py.File(path) as product:
                    record = product['records'][next(iter(product['records']))]
                    group = verification.create_group(f'{target}/{branch}')
                    group.attrs['fit_path'] = str(path)
                    group.attrs['fit_sha256'] = digest(path)
                    for name in ['range_m', 'fit/parameters', 'fit/standard_error', 'ion_fit/results_json']:
                        group[name.split('/')[-1]] = record[name][:]
                    fits[branch] = (record['range_m'][:], record['fit/parameters'][:])
            r, reference = fits['original']
            region = abs(r/1000-target) < 30
            for branch in ['injected', 'same_pulse', 'independent']:
                np.testing.assert_array_equal(fits[branch][0], r)
                delta = fits[branch][1]-reference
                print(target, branch, 'median change [ne fraction,Ti K,Te K,vi m/s]',
                      [np.nanmedian(delta[region, 0]/reference[region, 0]), *np.nanmedian(delta[region, 1:], axis=0)], flush=True)
        verification.attrs['complete'] = True
    print('MATCHED FIT CONTROL COMPLETE', flush=True)


if __name__ == '__main__':
    main()
