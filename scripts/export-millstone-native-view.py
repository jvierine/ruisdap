#!/usr/bin/env python3
"""Export native fits into the existing Millstone viewer; HDF5 remains authoritative."""
import argparse, hashlib, json, pathlib, shutil
import h5py
import numpy as np


def read_json(dataset):
    return json.loads(bytes(dataset[:]))


def clean(value):
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    return value


def write_json(path, value):
    path.write_text(json.dumps(clean(value), allow_nan=False, separators=(',', ':')))
    path.chmod(0o644)


def export(historical, fits, output, public):
    assert not output.exists(), 'Choose a new immutable output directory'
    original = json.loads((historical / 'manifest.json').read_text())
    manifest = dict(original)
    manifest.update(analysis_kind='native', analysis_id=output.name,
                    analysis_label='Rust experimental · 180 s ZENITH fits',
                    density_semantics='native_calibrated_ne_m3', fits=[], fit_times=[],
                    analysis_note='Experimental Rust fits of sparsely sampled integrations; O+ composition, 200–540 km. Transmit power is removed per input before averaging. Raw Ne and plasma-line panels retain the historical calibrated context. Velocity toward/away sign is not yet validated.')
    for key in ['raw_url', 'plasma_url']:
        url = original[key]
        manifest[key] = url if url.startswith('/') else '/guisdap/' + (historical / url).relative_to(public).as_posix()
    manifest['parameters_url'] = 'parameters.f32'
    output.mkdir(parents=True)
    output.chmod(0o755)
    all_parameters = []
    ranges = None
    for path in sorted(fits.glob('fit-*.h5')):
        with h5py.File(path) as h:
            r = h['records'][next(iter(h['records']))]
            config = read_json(r['configuration_json'])
            assert config['antenna'] == original['antenna']
            rg = r['range_m'][:]
            if ranges is None:
                ranges = rg
            assert np.array_equal(ranges, rg), 'Range grids differ'
            parameters, errors = r['fit/parameters'][:], r['fit/standard_error'][:]
            assert parameters.shape == errors.shape == (len(rg), 4)
            sources = read_json(h['sources/combined/records_json'])
            gates = read_json(r['ion_fit/results_json'])
            assumptions = read_json(r['ion_fit/assumptions_json'])
            tsys = []
            for source in sources:
                with h5py.File(source['path']) as s:
                    sr = s['records'][next(iter(s['records']))]
                    tsys.append(float(sr['ion_line/system_temperature_k'][0]))
            name = path.stem
            lag = r['ion_line/lag_s'][:]
            observed, variance = r['ion_line/acf'][:], r['ion_line/variance'][:]
            diagnostics = []
            statuses = []
            for i, range_m in enumerate(rg):
                gate = next((g for g in gates if g['range_m'] == range_m), None)
                result = gate['result'] if gate else None
                curve = result.get('model_curve', []) if result else []
                # Use exact fitted measurement rows, including any omitted/masked lags.
                model = [next((p['value'] for p in curve if abs(p['lag_s'] - t) < 1e-12), [None, None]) for t in lag]
                status = result['status'] if result else (gate['error'] if gate else 'outside fitted range')
                statuses.append(status)
                diagnostics.append(dict(range_m=float(range_m), observed=observed[i] * 1e6,
                    sigma=np.sqrt(np.where(variance[i] > 0, variance[i] / 2, np.nan)) * 1e6,
                    model=[[None if v is None else v * 1e6 for v in row] for row in model],
                    status=status, reduced_chi_square=result['reduced_chi_square'] if result else None))
            start, end = config['start_unix_s'], config['end_unix_s']
            exposure = sum(s['end_unix_s'] - s['start_unix_s'] for s in sources)
            write_json(output / (name + '-acf.json'), dict(lag_s=lag, units='µK W⁻¹ km⁻¹', negative_lag_convention=True, gates=diagnostics))
            sha = hashlib.sha256(path.read_bytes()).hexdigest()
            shutil.copyfile(path, output / path.name)
            (output / path.name).chmod(0o644)
            manifest['fits'].append(dict(id=name, start_unix_s=start, end_unix_s=end,
                parameters=parameters, standard_error=errors, fit_status=statuses,
                transmitter_w=float(np.mean([s['transmitter_w'] for s in sources])),
                system_temperature_k=float(np.mean(tsys)),
                azimuth_deg=float(np.mean([s['azimuth_deg'] for s in sources])),
                elevation_deg=float(np.mean([s['elevation_deg'] for s in sources])),
                source_integrations=len(sources), exposure_s=exposure,
                diagnostics_url=name + '-acf.json', hdf5_url=path.name, sha256=sha))
            manifest['fit_times'].append((start + end) / 2)
            all_parameters.append(parameters)
            manifest['native_calibration'] = dict(gain=assumptions['calibration'], units=assumptions['calibration_units'])
    assert all_parameters, 'No native fits'
    manifest['range_m'] = ranges
    manifest['analysis_default_window'] = [manifest['fits'][0]['start_unix_s'], manifest['fits'][-1]['end_unix_s']]
    # Viewer packing is [parameter, range, time], the same as historical exports.
    np.stack(all_parameters).transpose(2, 1, 0).astype('<f4').tofile(output / 'parameters.f32')
    (output / 'parameters.f32').chmod(0o644)
    write_json(output / 'manifest.json', manifest)
    print(f'Exported {len(all_parameters)} native fits to {output}')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--historical', type=pathlib.Path, required=True)
    p.add_argument('--fits', type=pathlib.Path, required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--public', type=pathlib.Path, default=pathlib.Path('/srv/guisdap'))
    a = p.parse_args()
    from isr_product import validate
    for path in sorted(a.fits.glob('fit-*.h5')):
        with h5py.File(path) as h:
            validate(h)
    export(a.historical, a.fits, a.output, a.public)
