#!/usr/bin/env python
"""Versioned ISR HDF5 interchange, source adapters and static browser export.

Run with conda run -n base python scripts/isr_product.py --help.
No processing-chain modules are imported and no raw data are modified.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import h5py
import numpy as np

SCHEMA = 'ruisdap-isr-v1'
PARAMETERS = ['density', 'ion_temperature', 'electron_temperature', 'velocity']
UNITS = ['source_proxy', 'K', 'K', 'm s^-1']


def text(value):
    if isinstance(value, bytes):
        return value.decode('utf-8')
    if isinstance(value, np.ndarray) and value.dtype == np.uint8:
        return bytes(value).decode('utf-8')
    return str(value)


def scalar(g, key, default=None):
    return g[key][()].item() if key in g else default


def metadata(g):
    return json.loads(text(g.attrs['metadata']))


def dataset(g, name, values, units, dimensions):
    a = np.asarray(values)
    if a.dtype.kind in 'fc' and np.isinf(a).any():
        raise ValueError(f'{name}: infinite values are unsupported; use NaN for missing')
    d = g.create_dataset(name, data=a, **({'compression': 'gzip', 'shuffle': True} if a.size > 128 else {}))
    d.attrs['units'] = units
    d.attrs['dimensions'] = dimensions
    return d


def source_record(h, path, key, meta):
    # Hash once per input file, and retain exact source bytes as a dataset tree.
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    sid = digest[:20]
    if sid not in h['sources']:
        s = h['sources'].create_group(sid)
        s.attrs['path'] = str(Path(path).resolve())
        s.attrs['sha256'] = digest
        with h5py.File(path, 'r') as original:
            for name in original:
                original.copy(name, s)
            s.attrs['original_attributes_json'] = json.dumps({k: text(v) for k, v in original.attrs.items()})
    rid = hashlib.sha256(f'{digest}:{key}'.encode()).hexdigest()[:24]
    if rid in h['records']:
        raise ValueError(f'duplicate source record: {path}:{key}')
    r = h['records'].create_group(rid)
    meta.update(id=rid, source_id=sid, source_key=key)
    r.attrs['metadata'] = json.dumps(meta, allow_nan=False)
    return r


def identity(station, antenna, channel, mode, frequency, kind, start, end, **extra):
    return dict(station=station, antenna=antenna, channel=channel, mode=str(mode),
                radar_frequency_hz=float(frequency), kind=kind,
                start_unix_s=float(start), end_unix_s=float(end), **extra)


def geometry(r, ranges_m, elevation=None, azimuth=None, altitude_m=None, site_height_m=0.):
    ranges_m = np.asarray(ranges_m, dtype=np.float64)
    dataset(r, 'range_m', ranges_m, 'm', 'range')
    if elevation is not None:
        r.attrs['elevation_deg'] = float(elevation)
        # Explicit spherical-Earth approximation; never silently rename range.
        if altitude_m is None:
            earth = 6370000. + site_height_m
            altitude_m = np.sqrt(earth**2 + ranges_m**2 + 2*earth*ranges_m*np.sin(np.deg2rad(elevation))) - 6370000.
            r.attrs['altitude_model'] = 'spherical Earth, 6370000 m radius'
    if azimuth is not None:
        r.attrs['azimuth_deg'] = float(azimuth)
    if altitude_m is not None:
        dataset(r, 'altitude_m', np.asarray(altitude_m, np.float64), 'm', 'range')


def millstone(h, path, args):
    with h5py.File(path) as src:
        if all(k in src for k in ['ne', 'Te', 'Ti', 'rgs', 't0', 't1']):
            if not args.antenna or not args.mode:
                raise ValueError('Millstone fit products require --antenna and --mode (not stored reliably in legacy files)')
            start, end = scalar(src, 't0'), scalar(src, 't1')
            r = source_record(h, path, '/', identity('Millstone Hill', args.antenna, args.channel or args.antenna.lower()+'-l', args.mode, args.radar_frequency_hz or 440.2e6, 'ion_fit', start, end,
                density_semantics='uncalibrated_ne_proxy', density_units='source_proxy',
                density_definition='fitted zero-lag power * (1+Te/Ti) * range_km^2 / transmitter_W; no magic constant applied',
                velocity_convention='source convention; sign not independently verified'))
            geometry(r, src['rgs'][()]*1000., scalar(src, 'el'), scalar(src, 'az'), site_height_m=146.)
            f = r.create_group('fit')
            vals = np.column_stack([src[k][()] for k in ['ne', 'Ti', 'Te', 'vi']])
            errs = np.full_like(vals, np.nan, dtype=np.float64)
            for i, k in enumerate(['dne', 'dTi', 'dTe', 'dvi']):
                if k in src and isinstance(src[k], h5py.Dataset):
                    errs[:, i] = src[k][()]
            dataset(f, 'parameters', vals, ','.join(UNITS), 'range,parameter')
            dataset(f, 'standard_error', errs, ','.join(UNITS), 'range,parameter')
            f.attrs['parameter_names'] = json.dumps(PARAMETERS)
            f.attrs['uncertainty_note'] = 'dTe/Ti is ratio uncertainty, not dTe; retained separately. No covariance invented.'
            if 'dTe/Ti' in src:
                dataset(f, 'temperature_ratio_standard_error', src['dTe/Ti'][()], '1', 'range')
            for k, unit in [('P_tx','W'), ('T_sys','K'), ('space_object_count','1'), ('range_avg_window_km','km'), ('range_avg_limits_km','km')]:
                if k in src:
                    dataset(f, k, src[k][()], unit, 'source-defined')
            return
        if all(k in src for k in ['acfs_e', 'lags', 'rgs_km', 'i0']):
            if not args.antenna or not args.mode or args.integration_s is None:
                raise ValueError('Millstone LPI requires --antenna, --mode and --integration-s (source stores only i0)')
            start = scalar(src, 'i0')
            r = source_record(h, path, '/', identity('Millstone Hill', args.antenna, args.channel or args.antenna.lower()+'-l', args.mode, args.radar_frequency_hz or 440.2e6, 'ion_acf', start, start+args.integration_s,
                signal_units='source_ACF', complex_convention='source complex ACF; no conjugation applied'))
            geometry(r, src['rgs_km'][()]*1000.)
            a = r.create_group('ion_line')
            dataset(a, 'lag_s', src['lags'][()], 's', 'lag')
            for key in ['acfs_e', 'acfs_g']:
                if key in src:
                    z = src[key][()]
                    dataset(a, key, np.stack([z.real,z.imag], axis=-1), 'source_ACF', 'range,lag,complex_component')
            if 'acfs_var' in src:
                v = src['acfs_var'][()]
                if np.iscomplexobj(v):
                    # Preserve the original in sources; do not mislabel complex variance.
                    a.attrs['variance_note'] = 'complex source variance retained in sources; no diagonal sigma inferred'
                else:
                    dataset(a, 'variance', v, 'source_ACF^2', 'range,lag')
            return
        groups = [(k,g) for k,g in src.items() if isinstance(g,h5py.Group) and 'spec' in g]
        if not groups:
            raise ValueError(f'{path}: unsupported Millstone layout')
        for key,g in groups:
            if not all(k in g for k in ['freqs','ranges','antenna','t1']):
                raise ValueError(f'{path}:{key}: missing spectrum axes or identity')
            antenna = text(g['antenna'][()])
            start, end = float(key), scalar(g,'t1')
            freq, ranges = g['freqs'][()], g['ranges'][()]*1000.
            # This Digital Metadata writer stores MHz, despite calling its
            # dataset freqs. Verify against receiver sampling before conversion.
            sr = scalar(g, 'sr')
            if sr is None:
                raise ValueError('Millstone plasma spectrum requires receiver sample rate to verify frequency units')
            if np.isclose(np.ptp(freq)*1e6, sr, rtol=0.002):
                freq = freq*1e6
            elif not np.isclose(np.ptp(freq), sr, rtol=0.002):
                raise ValueError('frequency span disagrees with receiver sample rate')
            spec = g['spec'][()]
            if spec.shape != (len(freq),len(ranges)):
                raise ValueError('spectrum dimensions do not match measured axes')
            order = np.argsort(freq)
            r = source_record(h, path, key, identity('Millstone Hill', antenna, text(g['channel'][()]) if 'channel' in g else antenna.lower()+'-h',
                f"pl-{scalar(g,'pl','unknown')}", args.radar_frequency_hz or 440.2e6, 'plasma_spectrum', start, end,
                signal_units='source_spectral_power', frequency_convention='signed receiver frequency offset; source bins retained and sorted',
                frequency_sign_verified=False))
            el = 90. if antenna == 'ZENITH' else scalar(g,'misa_elevation')
            az = 0. if antenna == 'ZENITH' else scalar(g,'misa_azimuth')
            geometry(r, ranges, el, az, site_height_m=146.)
            p = r.create_group('plasma_line')
            dataset(p, 'frequency_hz', freq[order], 'Hz', 'frequency')
            dataset(p, 'spectrum', spec[order].T, 'source_spectral_power', 'range,frequency')
            p.attrs['normalization'] = 'unchanged source power; no background subtraction or clipping'


def eiscat(h, path, args):
    with h5py.File(path) as src:
        if 'metadata_json' not in src:
            raise ValueError('EISCAT adapter expects processed BEATA HDF5; raw MAT remains handled by guisdap-worker')
        m = json.loads(bytes(src['metadata_json'][()]))
        if not m['schema'].startswith(('beata22-quicklook-', 'beata22-ion-input-')):
            raise ValueError(f"unsupported EISCAT layout {m['schema']}")
        r = source_record(h, path, '/', identity('EISCAT Tromsø UHF', 'UHF', 'UHF', 'BEATA 2.2', args.radar_frequency_hz or 930e6, 'legacy_eiscat', m['start_unix'], m['end_unix'],
            density_semantics='signed_raw_ne_proxy', density_units='m^-3',
            density_definition='source raw Ne equivalent; retain stored lag/zero-lag interpretation in source dataset units',
            frequency_convention='signed receiver offset; no sign reversal'))
        if 'ion_range_km' in src:
            geometry(r, src['ion_range_km'][()]*1000., m['elevation_deg'],m['azimuth_deg'],src['ion_altitude_km'][()]*1000.)
        else:
            # Minimal indexed input contains altitude only: range cannot be reconstructed honestly.
            dataset(r, 'altitude_m', src['ion_altitude_km'][()]*1000., 'm', 'range')
        ion = r.create_group('ion_line')
        for old,new,unit in [('ion_raw_ne','raw_density','m^-3'),('ion_raw_ne_sigma','raw_density_standard_error','m^-3')]:
            if old in src:
                dataset(ion,new,src[old][()],unit,'range')
        for prefix in ['long','short','uncoded']:
            old = f'native_{prefix}_acf'
            if old in src:
                a = ion.create_group(prefix)
                dataset(a,'acf',src[old][()],'source_accumulator / normalization','band,lag,range,complex_component')
                # Long and short sample intervals differ: copy sample indices and
                # preserve source timebase definition without guessing seconds.
                dataset(a,'lag_samples',src[f'native_{prefix}_lag_samples'][()],'sample','lag')
                a.attrs['timebase_note'] = 'sample interval requires BEATA response definition; retained in source experiment metadata'
        if 'decoded_spectrum' in src:
            for band in range(src['decoded_spectrum'].shape[0]):
                p = r.create_group(f'plasma_line_band_{band}')
                dataset(p,'frequency_hz',src['decoded_frequency_hz'][band],'Hz','frequency')
                dataset(p,'range_m',src['decoded_range_km'][()]*1000.,'m','range')
                dataset(p,'altitude_m',src['decoded_altitude_km'][()]*1000.,'m','range')
                dataset(p,'spectrum',src['decoded_spectrum'][band],'source_signed_fft_contrast','range,frequency')
                p.attrs['normalization'] = 'source signed FFT contrast; missing zero lag, not calibrated spectral power'


def validate(h):
    if text(h.attrs.get('schema')) != SCHEMA or 'records' not in h or 'sources' not in h:
        raise ValueError('unsupported ISR schema')
    if not len(h['records']):
        raise ValueError('empty product')
    for rid,r in h['records'].items():
        m = metadata(r)
        if not (np.isfinite([m['start_unix_s'],m['end_unix_s'],m['radar_frequency_hz']]).all() and m['end_unix_s'] > m['start_unix_s'] and m['radar_frequency_hz']>0):
            raise ValueError(f'{rid}: invalid time/frequency')
        if not all(m.get(k) for k in ['station','antenna','channel','mode','kind']):
            raise ValueError(f'{rid}: incomplete identity')
        if m['source_id'] not in h['sources']:
            raise ValueError('missing source provenance')
        for g in [r]+[x for x in r.values() if isinstance(x,h5py.Group)]:
            for axis in ['range_m','altitude_m','frequency_hz','lag_s']:
                if axis in g:
                    a = g[axis][()]
                    if a.ndim != 1 or not np.isfinite(a).all() or (np.diff(a)<=0).any():
                        raise ValueError(f'{rid}/{axis}: axis must be finite and strictly increasing')
            if 'spectrum' in g:
                rg = g['range_m'] if 'range_m' in g else r['range_m']
                if g['spectrum'].shape != (len(rg),len(g['frequency_hz'])):
                    raise ValueError('spectrum axes mismatch')
        if 'fit' in r:
            n = len(r['range_m'])
            if r['fit/parameters'].shape != (n,4) or r['fit/standard_error'].shape != (n,4):
                raise ValueError('fit shape mismatch')
            e = r['fit/standard_error'][()]
            if (e[np.isfinite(e)]<0).any():
                raise ValueError('negative standard errors')
    return len(h['records'])


def atomic_hdf5(destination, callback):
    destination = Path(destination)
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():
        raise ValueError(f'{destination} exists; choose a new product version')
    fd,tmp = tempfile.mkstemp(dir=destination.parent,suffix='.h5')
    os.close(fd)
    try:
        with h5py.File(tmp,'w') as h:
            h.attrs['schema'] = SCHEMA
            h.attrs['time_convention'] = 'UTC Unix seconds; half-open integration intervals [start,end)'
            h.attrs['missing_values'] = 'NaN; signed values retained; absent uncertainty is not zero'
            h.create_group('records'); h.create_group('sources'); h.create_group('calibration')
            callback(h)
            count = validate(h)
        os.replace(tmp,destination)
        return count
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def clean(value):
    if isinstance(value,np.ndarray): return clean(value.tolist())
    if isinstance(value,np.generic): return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value): return None
    if isinstance(value,list): return [clean(v) for v in value]
    if isinstance(value,dict): return {k:clean(v) for k,v in value.items()}
    return value


def export(input_path, output):
    output = Path(output)
    if output.exists(): raise ValueError('export destination exists; choose a new version')
    output.mkdir(parents=True)
    entries=[]
    with h5py.File(input_path) as h:
        validate(h)
        for rid,r in h['records'].items():
            m = metadata(r)
            detail = dict(m, elevation_deg=r.attrs.get('elevation_deg'), azimuth_deg=r.attrs.get('azimuth_deg'))
            for name in ['range_m','altitude_m']:
                if name in r: detail[name] = r[name][()]
            if 'fit' in r:
                f = r['fit']
                detail['fit'] = dict(parameters=f['parameters'][()],standard_error=f['standard_error'][()])
                for k in ['P_tx','T_sys','space_object_count','temperature_ratio_standard_error']:
                    if k in f: detail['fit'][k] = f[k][()]
            if 'ion_line/raw_density' in r:
                detail['raw_density'] = r['ion_line/raw_density'][()]
            spectra=[]
            for name,g in r.items():
                if name.startswith('plasma_line') and isinstance(g,h5py.Group) and 'spectrum' in g:
                    values = g['spectrum'][()]
                    # Full-resolution slices remain in HDF5. Browser mesh is
                    # averaged in frequency only, with its actual bin centres.
                    step = max(1,int(np.ceil(values.shape[1]/2400)))
                    trim = values.shape[1]//step*step
                    a = values[:,:trim].reshape(len(values),-1,step).mean(axis=2)
                    freq = g['frequency_hz'][:trim].reshape(-1,step).mean(axis=1)
                    filename=f'{rid}-{name}.f32'
                    a.astype('<f4').tofile(output/filename)
                    ranges = g['range_m'][()] if 'range_m' in g else r['range_m'][()]
                    alt = g['altitude_m'][()] if 'altitude_m' in g else (r['altitude_m'][()] if 'altitude_m' in r else None)
                    spectra.append(dict(name=name,url=filename,shape=list(a.shape),frequency_hz=freq,range_m=ranges,altitude_m=alt,units=text(g['spectrum'].attrs['units']),frequency_average_bins=step))
            detail['spectra'] = spectra
            (output/f'{rid}.json').write_text(json.dumps(clean(detail),allow_nan=False))
            # Source paths/build details remain in canonical HDF5, not visible UI.
            entries.append({k:v for k,v in m.items() if k not in ['source_id','source_key']} | dict(url=f'{rid}.json',has_spectra=bool(spectra),has_fit='fit' in r,has_raw_density='ion_line/raw_density' in r))
    entries.sort(key=lambda e:(e['start_unix_s'],e['station'],e['antenna'],e['kind']))
    (output/'index.json').write_text(json.dumps(dict(schema=SCHEMA,records=entries),allow_nan=False))
    return len(entries)


def calibrate(products, review_path, destination):
    """Recompute reviewed points from immutable HDF5 sources; never trust a
    browser-supplied gain, density, temperature, or cached quality result."""
    review=json.loads(Path(review_path).read_text())
    if review.get('schema')!='ruisdap-calibration-review-v1':raise ValueError('unsupported review schema')
    rows=[];records={}
    for path in products:
        with h5py.File(path)as h:
            validate(h)
            for rid,r in h['records'].items():
                if rid in records:raise ValueError('duplicate record identity in product set')
                records[rid]=(path,metadata(r))
    for point in review['points']:
        if not point.get('included'):continue
        oid,fid=point['observation_id'],point['fit_id']
        if oid not in records or fid not in records:raise ValueError('review references a record outside the supplied products')
        op,om=records[oid];fp,fm=records[fid]
        if (om['station'],om['antenna'])!=(fm['station'],fm['antenna']):raise ValueError('cross-instrument calibration pairing')
        if fm.get('density_semantics')!='uncalibrated_ne_proxy':raise ValueError('unsupported density normalization')
        dt=abs((om['start_unix_s']+om['end_unix_s']-fm['start_unix_s']-fm['end_unix_s'])/2)
        if dt>300:raise ValueError('calibration pairing exceeds 300 seconds')
        frequency=float(point['frequency_hz']);sigma=float(point['frequency_standard_error_hz']);gate=int(point['gate']);selected=float(point['range_m'])
        if not np.isfinite([frequency,sigma,selected]).all()or sigma<=0:raise ValueError('invalid frequency/range uncertainty')
        with h5py.File(op)as h:
            r=h['records'][oid];spectra=[g for g in r.values()if isinstance(g,h5py.Group)and 'frequency_hz'in g and 'spectrum'in g]
            if not any(g['frequency_hz'][0]<=frequency<=g['frequency_hz'][-1] and r['range_m'][0]<=selected<=r['range_m'][-1]for g in spectra):raise ValueError('selection outside plasma spectrum')
        with h5py.File(fp)as h:
            r=h['records'][fid];ranges=r['range_m'][()]
            if gate<0 or gate>=len(ranges)or not np.isclose(selected,ranges[gate],atol=0.01):raise ValueError('review gate/range mismatch')
            raw,ti,te,velocity=r['fit/parameters'][gate]
            if not np.isfinite([raw,ti,te]).all()or min(raw,ti,te)<=0:raise ValueError('nonphysical fit')
            eps=8.8541878128e-12;kb=1.380649e-23;me=9.1093837015e-31;e=1.602176634e-19
            k=4*np.pi*om['radar_frequency_hz']/299792458.
            density=((2*np.pi*frequency)**2-3*k*k*kb*te/me)*eps*me/e**2
            if density<=0:raise ValueError('thermal term exceeds selected frequency squared')
            q=k*k*eps*kb*te/(density*e**2);denom=(1+q)*(1+q+te/ti)
            multiplier=density*(1+te/ti)/denom/raw
            rows.append(dict(station=om['station'],antenna=om['antenna'],mode=fm['mode'],observation_id=oid,fit_id=fid,gate=gate,range_m=selected,frequency_hz=frequency,frequency_standard_error_hz=sigma,plasma_density_m3=density,k2_debye2=q,multiplier=multiplier,fit_start_unix_s=fm['start_unix_s']))
    if not rows:raise ValueError('no included reviewed calibration points')
    groups={}
    for row in rows:groups.setdefault((row['station'],row['antenna'],row['mode']),[]).append(row)
    def write(h):
        for path in products:
            with h5py.File(path)as original:
                for rid in original['records']:original.copy(original['records'][rid],h['records'],name=rid)
                for sid in original['sources']:
                    if sid not in h['sources']:original.copy(original['sources'][sid],h['sources'],name=sid)
        h['calibration'].attrs['method']='Bohm-Gross density; full finite-Debye integrated ion power; median multiplier, separately per station/antenna/fit-mode'
        h['calibration'].attrs['review_sha256']=hashlib.sha256(Path(review_path).read_bytes()).hexdigest()
        h['calibration'].attrs['application']='record only; no automatic rescaling of input fits'
        for i,(key,points)in enumerate(groups.items()):
            g=h['calibration'].create_group(f'group_{i:03}')
            for name,value in zip(['station','antenna','mode'],key):g.attrs[name]=value
            g.attrs['multiplier']=float(np.median([p['multiplier']for p in points]))
            g.attrs['independent_fit_count']=len(set(p['fit_id']for p in points))
            g.attrs['uncertainty']='frequency errors retained; no independent estimate of total calibration uncertainty supplied'
            g.attrs['points_json']=json.dumps(points,allow_nan=False)
    return atomic_hdf5(destination,write)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    imp=sub.add_parser('import'); imp.add_argument('--adapter',choices=['millstone','eiscat'],required=True)
    imp.add_argument('--output',required=True); imp.add_argument('--antenna'); imp.add_argument('--channel'); imp.add_argument('--mode')
    imp.add_argument('--radar-frequency-hz',type=float); imp.add_argument('--integration-s',type=float)
    imp.add_argument('inputs',nargs='+')
    v=sub.add_parser('validate'); v.add_argument('input')
    e=sub.add_parser('export'); e.add_argument('input'); e.add_argument('--output',required=True)
    c=sub.add_parser('calibrate');c.add_argument('--products',nargs='+',required=True);c.add_argument('--review',required=True);c.add_argument('--output',required=True)
    a=p.parse_args()
    if a.command=='import':
        adapter={'millstone':millstone,'eiscat':eiscat}[a.adapter]
        def run(h):
            for path in a.inputs: adapter(h,path,a)
        print(f'{atomic_hdf5(a.output,run)} records imported')
    elif a.command=='calibrate':print(f'{calibrate(a.products,a.review,a.output)} source records retained with reviewed calibration')
    elif a.command=='validate':
        with h5py.File(a.input) as h: print(f'{validate(h)} valid records')
    else: print(f'{export(a.input,a.output)} records exported')

if __name__=='__main__': main()
