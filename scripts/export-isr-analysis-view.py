#!/usr/bin/env python3
"""Export one day of isr_analysis fits into the Millstone viewer as an analysis.

isr_analysis (github.com/ionodev/isr_analysis, FYS-3730 project) fits the
10-second lag-profile-inversion ACFs of one run with fit_lpi.py into
lpi_<postfix>/<channel>/pp-<t0>.h5. Its runs differ in how satellite echoes and
impulse noise are handled (satellite echoes as extra unknowns in the
inversion, pass tracking, removal of impulse-hit rows). This script writes one
run's day in the viewer's quicklook schema, so that runs can be compared in
the same maps as the historical and native analyses:

    parameters.f32   Ne, Ti, Te, velocity, [parameter][range][time]
    raw-ne.f32       raw Ne from the run's own 10 microsecond lag, [range][time]
    manifest.json    the viewer's manifest; fits carry their 1-sigma errors
    quicklook.h5     the same arrays with units and source hashes

Density. The fitted zero-lag power P0 (receiver gain removed, in K) gives the
proxy p = P0 (1 + Te/Ti) r_km^2 / P_tx, stored by fit_lpi.py as "ne". The
viewer computes Ne from C p with its finite-Debye correction. C is the run's
own density constant (W m^-3 K^-1 km^-2): it depends on the processing (gate
summation, filter), so isr_analysis's lpi_60 runs use their digisonde value
(7.15e8, isr_analysis memo 43), not the constant of the 30 us historical
products.

Raw Ne is 2 C Re ACF(10 us) / alpha * r_km^2 / P_tx per 10-second inversion
(alpha the receiver gain, Te/Ti = 1), with the ground-clutter-cancelled ACF
below 300 km as fit_lpi.py uses it. The plasma-line panel is the day's
existing export (the plasma lines do not depend on the ion-line processing);
its manifest fields are copied and its spectrum file is referenced, not copied.

Usage:
    python3 scripts/export-isr-analysis-view.py RUN OUTPUT --date 2024-04-08 \
        --plasma-manifest BASE/manifest.json --plasma-url /guisdap/data/BASE/plasma.f32 \
        --label "..." [--magic-constant 7.15e8] [--postfix lpi_60] [--channel zenith-l]
"""
import argparse, datetime as dt, hashlib, json
from pathlib import Path

import h5py
import numpy as np

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('run', type=Path, help='isr_analysis output root (contains lpi_<postfix>/<channel>)')
p.add_argument('output', type=Path, help='new export directory')
p.add_argument('--date', required=True)
p.add_argument('--postfix', default='lpi_60')
p.add_argument('--channel', default='zenith-l')
p.add_argument('--antenna', default='ZENITH')
p.add_argument('--magic-constant', type=float, default=7.15e8)
p.add_argument('--calibration-method', default='Digisonde foF2 against the fitted proxy at the F2 peak '
               '(isr_analysis memo 43), for lpi_60 products with a 100 kHz pass band and filter length 20')
p.add_argument('--plasma-manifest', type=Path, required=True, help="manifest.json of the day's existing export")
p.add_argument('--plasma-url', required=True, help='URL of that export\'s plasma spectrum file')
p.add_argument('--label', required=True, help='what the run is, for the manifest')
p.add_argument('--source-commit', default='', help='isr_analysis commit of the run')
a = p.parse_args()

if a.output.exists():
    raise SystemExit('choose a new export version: %s exists' % a.output)
start = dt.datetime.fromisoformat(a.date).replace(tzinfo=dt.timezone.utc).timestamp()
end = start + 86400
folder = a.run / a.postfix / a.channel
C = a.magic_constant


def finite(x):
    if isinstance(x, np.ndarray):
        return finite(x.tolist())
    if isinstance(x, (list, tuple)):
        return [finite(v) for v in x]
    if isinstance(x, dict):
        return {k: finite(v) for k, v in x.items()}
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.generic):
        return x.item()
    return x


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


fits, sources = [], []
for path in sorted(folder.glob('pp-*.h5')):
    t = float(path.stem.split('-')[1])
    if not start <= t < end:
        continue
    with h5py.File(path, 'r') as h:
        if not all(k in h for k in ['ne', 'Ti', 'Te', 'vi', 'rgs', 't0', 't1', 'P_tx']):
            continue
        proxy = h['ne'][()]
        params = np.stack([proxy * C, h['Ti'][()], h['Te'][()], h['vi'][()]], axis=1)
        err = np.stack([h[k][()] for k in ['dne', 'dTi', 'dTe', 'dvi']], axis=1)
        uerr = err.copy()
        err[:, 0] *= C
        digest = sha(path)
        fits.append(dict(
            id=hashlib.sha256((digest + ':/').encode()).hexdigest()[:24],
            start_unix_s=float(h['t0'][()]), end_unix_s=float(h['t1'][()]),
            range_m=h['rgs'][()] * 1000.0, parameters=params, standard_error=err,
            uncalibrated_density_proxy=proxy, uncalibrated_standard_error=uerr,
            transmitter_w=float(h['P_tx'][()]), system_temperature_k=float(h['T_sys'][()]),
            elevation_deg=float(h['el'][()]), azimuth_deg=float(h['az'][()]),
            space_object_count=h['space_object_count'][()] if 'space_object_count' in h else [],
            n_lpi_files=int(h['n_lpi_files'][()]) if 'n_lpi_files' in h else None))
        sources.append(dict(path=str(path), sha256=digest))
if not fits:
    raise SystemExit('no fits of %s in %s' % (a.date, folder))
fits.sort(key=lambda f: f['start_unix_s'])
ranges = fits[0]['range_m']
if any(not np.array_equal(f['range_m'], ranges) for f in fits):
    raise SystemExit('fit range grids differ')

raw, rawranges = [], None
for path in sorted(folder.glob('lpi-*.h5')):
    t = float(path.stem.split('-')[1])
    if not start <= t < end:
        continue
    with h5py.File(path, 'r') as h:
        r = h['rgs_km'][()]
        lags = h['lags'][()]
        lag = int(np.argmin(abs(lags - 1e-5)))
        if abs(lags[lag] - 1e-5) > 1e-9:
            continue
        acf = h['acfs_e'][:, lag].copy()
        cut = np.where(r > 300)[0][0]
        acf[:cut] = h['acfs_g'][:cut, lag] / 2.0
        power, alpha = float(h['P_tx'][()]), float(h['alpha'][()])
        values = 2 * C * acf.real / alpha * r * r / power if power > 0 and alpha > 0 else np.full(len(r), np.nan)
        time = float(h['t0'][()]) if 't0' in h else float(h['i0'][()])
    if rawranges is None:
        rawranges = r * 1000.0
    elif not np.array_equal(rawranges, r * 1000.0):
        raise SystemExit('raw range grids differ')
    raw.append((time, values))
raw.sort(key=lambda x: x[0])
rawtimes = np.array([x[0] for x in raw])
rawpixels = np.array([x[1] for x in raw]).T

base = json.loads(a.plasma_manifest.read_text())
if base['date'] != a.date or base['antenna'] != a.antenna:
    raise SystemExit('the plasma manifest is for another day or antenna')

times = np.array([(f['start_unix_s'] + f['end_unix_s']) / 2 for f in fits])
parameters = np.array([f['parameters'] for f in fits])
a.output.mkdir(parents=True)
np.asarray(parameters.transpose(2, 1, 0), dtype='<f4').tofile(a.output / 'parameters.f32')
np.asarray(rawpixels, dtype='<f4').tofile(a.output / 'raw-ne.f32')
manifest = dict(
    schema='ruisdap-millstone-quicklook-v1', station='Millstone Hill', antenna=a.antenna, date=a.date,
    radar_frequency_hz=440.2e6, fit_times=times, range_m=ranges, parameters_url='parameters.f32',
    parameter_names=['Ne / m⁻³', 'Ti / K', 'Te / K', 'Velocity / m s⁻¹'], fits=fits,
    raw_times=rawtimes, raw_range_m=rawranges, raw_url='raw-ne.f32',
    raw_definition='Raw Ne = 2 C Re ACF(10 us) / receiver gain * range_km^2 / transmitter_W from each '
                   '10-second isr_analysis inversion; Te/Ti = 1; ground-clutter-cancelled ACF below 300 km',
    raw_units='m^-3',
    plasma_times=base['plasma_times'], plasma_frequency_hz=base['plasma_frequency_hz'],
    plasma_url=a.plasma_url, plasma_records=base['plasma_records'],
    collapse_range_km=base['collapse_range_km'], collapse_definition=base['collapse_definition'],
    density_semantics='digisonde_calibrated_density',
    frequency_convention=base.get('frequency_convention', 'signed receiver offset; no sign reversal'),
    frequency_sign_verified=base.get('frequency_sign_verified', False),
    velocity_convention='isr_analysis fit_lpi.py convention; sign not independently verified',
    calibration=dict(gain=C, revision=1, magic_constant_units='W m^-3 K^-1 km^-2', points=[],
                     method=a.calibration_method),
    isr_analysis=dict(label=a.label, run=str(a.run), postfix=a.postfix, channel=a.channel,
                      commit=a.source_commit, fit_count=len(fits), raw_count=len(raw)))
(a.output / 'manifest.json').write_text(json.dumps(finite(manifest), allow_nan=False))
with h5py.File(a.output / 'quicklook.h5', 'w') as h:
    h.attrs['schema'] = manifest['schema']
    h.attrs['metadata_json'] = json.dumps(finite({k: v for k, v in manifest.items() if k not in ['fits', 'plasma_records']}), allow_nan=False)
    h.attrs['sources_json'] = json.dumps(sources)
    for k, v, unit in [('fit_time_unix_s', times, 'UTC Unix s'), ('range_m', ranges, 'm'),
                       ('parameters', parameters, 'm^-3,K,K,m s^-1'),
                       ('standard_error', np.array([f['standard_error'] for f in fits]), 'm^-3,K,K,m s^-1'),
                       ('raw_time_unix_s', rawtimes, 'UTC Unix s'), ('raw_range_m', rawranges, 'm'),
                       ('raw_ne', rawpixels, 'm^-3')]:
        h.create_dataset(k, data=v, compression='gzip')
        h[k].attrs['units'] = unit
for f in a.output.iterdir():
    f.chmod(0o644)
print('Export complete: %s %s, %d fits, %d raw profiles -> %s' % (a.date, a.antenna, len(fits), len(raw), a.output))
