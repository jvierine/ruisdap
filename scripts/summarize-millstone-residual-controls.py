#!/usr/bin/env python3
"""Summarize paired synthetic residual controls, retaining bootstrap intervals."""
import argparse
import hashlib
import pathlib

import h5py
import numpy as np


def paired_ratio(error, reference, rng, draws=2000):
    use = np.isfinite(error) & np.isfinite(reference)
    error, reference = error[use], reference[use]
    if len(error) < 2 or reference.mean() <= 0:
        return np.full(3, np.nan)
    index = rng.integers(0, len(error), (draws, len(error)))
    den = reference[index].mean(axis=1)
    ratios = np.divide(error[index].mean(axis=1), den, out=np.full(draws, np.nan), where=den > 0)
    return np.r_[error.mean()/reference.mean(), np.nanpercentile(ratios, [2.5, 97.5])]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    a = p.parse_args()
    assert not a.output.exists()
    rng = np.random.default_rng(514312)
    with h5py.File(a.input) as source, h5py.File(a.output, 'x') as report:
        assert source.attrs['complete']
        report.attrs['source_sha256'] = hashlib.sha256(a.input.read_bytes()).hexdigest()
        report.attrs['script_sha256'] = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
        report.attrs['scope'] = 'Paired direct-voltage lag moments with identical masks for processed and ideal signals. Not a full LPI inversion, empirical campaign false-alarm estimate, or proof of target identity.'
        report.attrs['interval_definition'] = 'Estimate and 2.5/97.5 percentiles of 2000 paired trial bootstrap ratios; conditional on synthetic model.'
        report['lag_samples'] = source['lags_samples'][:]
        lag_index = int(np.flatnonzero(source['lags_samples'][:] == 10)[0])
        for family in source:
            if not isinstance(source[family], h5py.Group):
                continue
            for width in source[family]:
                for case in source[family][width]:
                    path = f'{family}/{width}/{case}'
                    g = source[path]
                    out = report.create_group(path)
                    out.attrs.update(g.attrs)
                    rows = [r for r in g.values() if 'diffuse_acf' in r]
                    if not rows:
                        continue
                    reference = np.array([r['diffuse_acf'][:].real for r in rows])
                    for kind in ['baseline', 'joint']:
                        errors = np.array([r[kind+'_acf_error'][:] for r in rows])
                        out[kind+'_mean_acf_error'] = np.nanmean(errors, axis=0)
                        out[kind+'_lag10_relative_error_interval'] = paired_ratio(errors[:, lag_index].real, reference[:, lag_index], rng)
                        out[kind+'_masked_fraction'] = [r.attrs[kind+'_masked_fraction'] for r in rows]
                        out[kind+'_distortion_power'] = [r.attrs[kind+'_distortion_power'] for r in rows]
                    out['mean_diffuse_acf'] = reference.mean(axis=0)
                    # Directly paired difference isolates the proposed change
                    # from distortion already present in baseline subtraction.
                    delta = np.array([(r['joint_acf_error'][lag_index]-r['baseline_acf_error'][lag_index]).real for r in rows])
                    out['joint_minus_baseline_lag10_interval'] = paired_ratio(delta, reference[:, lag_index], rng)
                    b = out['baseline_lag10_relative_error_interval'][0]
                    j = out['joint_lag10_relative_error_interval'][0]
                    print(path, f'false/true applied {g.attrs["initial_accepted_subtractions"]}->{g.attrs["final_accepted_subtractions"]}; joint refits {g.attrs["joint_refits"]}; residual {g.attrs["residual_layer"]}; lag10 error/reference {b:.4g}->{j:.4g}', flush=True)
        report.attrs['complete'] = True


if __name__ == '__main__':
    main()
