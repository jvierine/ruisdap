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



def debris_events(record):
    """Per-science-pulse CFAR/refined detections, never repeated clutter partners.
    Refined delay is relative to the measured transmit leakage and has had the
    measured transmitter-to-receiver timing offset removed in the archive reader.
    """
    if 'debris/diagnostics_json' not in record:return []
    config=read_json(record['configuration_json'])
    if not config.get('cancel_debris',True):return []
    diagnostics=read_json(record['debris/diagnostics_json'])
    samples=record['pulse_sample'][:]
    science=config['lpi'].get('science_pulse') or [True]*len(samples)
    assert len(diagnostics)==len(samples)==len(science),'Debris/pulse identity dimensions differ'
    sr=config['lpi']['sample_rate_hz'];assert sr>0
    events=[]
    for sample,use,diagnostic in zip(samples,science,diagnostics):
        if not use:continue
        # A failed joint validation masks the contaminated support even if a
        # subtraction was attempted. Do not call that an accepted subtraction.
        status=(2 if not diagnostic['accepted'] else
                0 if diagnostic['applied'] and diagnostic['confirmed'] else
                1 if diagnostic['applied'] else 3)
        shrinkage=diagnostic.get('amplitude_shrinkage')
        if shrinkage is not None:assert len(shrinkage) in [0,len(diagnostic['detections'])], 'Debris amplitude dimensions differ'
        for index,detection in enumerate(diagnostic['detections']):
            event_status=3 if diagnostic['accepted'] and shrinkage and shrinkage[index]==0 else status
            # Joint pulse acceptance does not confirm every candidate. A weak
            # amplitude must remain weak even alongside a strong fitted echo.
            snr=diagnostic.get('amplitude_snr_db',[])
            if event_status==0 and index<len(snr) and snr[index]<14:
                event_status=1
            delay=detection['delay_samples'];doppler=detection['doppler_hz']
            assert np.isfinite(delay) and delay>=0 and np.isfinite(doppler),'Invalid debris coordinates'
            events.append((int(sample)/sr,delay/sr*299792.458/2,event_status,doppler))
    return events

def export(historical, fits, output, public):
    assert not output.exists(), 'Choose a new immutable output directory'
    original = json.loads((historical / 'manifest.json').read_text())
    manifest = dict(original)
    manifest.update(analysis_kind='native', analysis_id=output.name,
                    analysis_label='Rust experimental · 180 s ZENITH fits',
                    density_semantics='native_calibrated_ne_m3', fits=[], fit_times=[],
                    analysis_note='Rust ACF fits. Recorded transmit power is removed per input before averaging; pulse energy and range measure are explicit. Raw Ne and plasma-line panels retain the source context. Velocity toward/away sign is not yet validated.')
    for key in ['raw_url', 'plasma_url']:
        url = original[key]
        manifest[key] = url if url.startswith('/') else '/guisdap/' + (historical / url).relative_to(public).as_posix()
    manifest['parameters_url'] = 'parameters.f32'
    output.mkdir(parents=True)
    output.chmod(0o755)
    all_parameters = []
    ranges = None
    debris=[];debris_sources=set()
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
                    if source['path'] not in debris_sources:
                        debris.extend(debris_events(sr));debris_sources.add(source['path'])
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
            # Coded and long-pulse observations share time windows. Sum their
            # union, not their overlapping spans, when reporting wall exposure.
            intervals=sorted((s['start_unix_s'],s['end_unix_s'])for s in sources)
            exposure=0.;last_end=-np.inf
            for lo,hi in intervals:
                exposure+=max(0.,hi-max(lo,last_end));last_end=max(last_end,hi)
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
                source_integrations=len(sources), exposure_s=exposure,science_pulses=sum(s.get('science_pulse_count')or 0 for s in sources),
                diagnostics_url=name + '-acf.json', hdf5_url=path.name, sha256=sha))
            manifest['fit_times'].append((start + end) / 2)
            all_parameters.append(parameters)
            manifest['native_integration_s']=assumptions.get('integration_s',180)
            manifest['native_historical_settings']=assumptions.get('historical_settings',False)
            manifest['native_calibration'] = dict(gain=assumptions['calibration'], units=assumptions['calibration_units'])
    assert all_parameters, 'No native fits'
    manifest['range_m'] = ranges
    manifest['analysis_default_window'] = [manifest['fits'][0]['start_unix_s'], manifest['fits'][-1]['end_unix_s']]
    # Viewer packing is [parameter, range, time], the same as historical exports.
    np.stack(all_parameters).transpose(2, 1, 0).astype('<f4').tofile(output / 'parameters.f32')
    (output / 'parameters.f32').chmod(0o644)
    # UTC-relative float32 coordinates keep sub-second precision on the GPU;
    # authoritative pulse keys and Doppler/quality diagnostics remain in HDF5.
    if 'date' in original:
        import datetime as dt
        origin=dt.datetime.fromisoformat(original['date']).replace(tzinfo=dt.timezone.utc).timestamp()
        events=np.asarray(debris,dtype=np.float64).reshape(-1,4)
        if len(events):events[:,0]-=origin
        events.astype('<f4').tofile(output/'debris.f32');(output/'debris.f32').chmod(0o644)
        manifest.update(debris_url='debris.f32',debris_count=len(events),debris_time_origin_unix_s=origin,
            debris_columns=['time_offset_s','slant_range_km','status','doppler_hz'],
            debris_status=['Subtracted, confirmed','Weak subtraction','Masked after failed validation','Low-SNR candidate, retained'],
            debris_definition='CFAR-seeded detections in science pulses, excluding duplicate clutter partners. Refined delay is corrected for measured TX/RX channel timing. Candidate detections are not independent confirmations of an orbiting object.')
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
