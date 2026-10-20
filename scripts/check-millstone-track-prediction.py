#!/usr/bin/env python3
"""Held-pulse forecasts from other pulses' saved native amplitude diagnostics.

Metadata-only diagnostic. The queried pulse is excluded from association and
all parameter fits. Its saved fit is used afterwards to score the forecast,
never to select or modify it. This is not yet production subtraction.
"""
import argparse
import hashlib
import json
import pathlib

import h5py
import numpy as np

import millstone_track_prediction as track


def events_from_record(record, configuration):
    keys = record['pulse_sample'][:]
    qa = json.loads(bytes(record['debris/diagnostics_json'][:]))
    assert len(keys) == len(qa) == len(configuration['lpi']['science_pulse'])
    events = []
    for i, q in enumerate(qa):
        if not configuration['lpi']['science_pulse'][i] or q['reason'] in ['clipped_adc', 'invalid_background_noise', 'background_gate_contaminated']:
            continue
        for detection, amplitude, snr in zip(q['detections'], q['amplitudes'], q['amplitude_snr_db']):
            events.append(dict(sample=int(keys[i]), delay_samples=detection['delay_samples'],
                               doppler_hz=detection['doppler_hz'], amplitude=complex(*amplitude), snr_db=snr))
    return keys, qa, events


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--product', type=pathlib.Path, action='append', required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--queries', type=int, default=8)
    a = p.parse_args()
    assert not a.output.exists() and a.queries > 0
    with h5py.File(a.output, 'x') as report:
        report.attrs['complete'] = False
        report.attrs['scope'] = __doc__
        for name, path in [('script', pathlib.Path(__file__)), ('predictor', pathlib.Path(track.__file__))]:
            report.attrs[name+'_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        for path in a.product:
            with h5py.File(path) as source:
                for record_id, record in source['records'].items():
                    c = json.loads(bytes(record['configuration_json'][:]))
                    keys, qa, events = events_from_record(record, c)
                    selected = np.flatnonzero(c['lpi']['science_pulse'])
                    held = selected[np.linspace(0, len(selected)-1, min(a.queries, len(selected))).astype(int)]
                    out = report.create_group(path.stem+'/'+record_id)
                    out.attrs['source_path'] = str(path)
                    out.attrs['source_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
                    # Long pulses occur every 26.73 ms: a 0.5 s neighbourhood
                    # can supply sixteen other pulses even at a source edge.
                    radius = .5 if c['mode'] == 'long-pulse-480us' else .35
                    out.attrs['prediction_radius_s'] = radius
                    out['configuration_json'] = np.frombuffer(json.dumps(c).encode(), np.uint8)
                    for lo, hi in [(250., 500.), (525., 575.)]:
                        # The range bands are fixed in advance, not inferred
                        # from the held pulse's data or detection coordinates.
                        training_region = [v for v in events if lo <= v['delay_samples']*.149896229 <= hi]
                        group = out.create_group(f'range_{lo}_{hi}')
                        scores = []
                        for index in held:
                            sample = int(keys[index])
                            predictions = track.predict(training_region, sample, np.mean(c['lpi']['tx_gate']), radius_s=radius)
                            row = group.create_group(str(sample))
                            values = []
                            for prediction in predictions:
                                if prediction['supported']:
                                    assert sample not in prediction['training_samples']
                                    candidates = [e for e in training_region if e['sample'] == sample]
                                    # Assessment only, strictly after prediction.
                                    if candidates:
                                        measured = min(candidates, key=lambda v: ((v['delay_samples']-prediction['delay_samples'])/8)**2+((v['doppler_hz']-prediction['doppler_hz'])/400)**2)
                                        centre = measured['amplitude']*np.exp(2j*np.pi*measured['doppler_hz']*(np.mean(c['lpi']['tx_gate'])+measured['delay_samples'])/1e6)
                                        predicted = complex(*prediction['centre_amplitude'])
                                        prediction['assessment'] = dict(delay_error_samples=prediction['delay_samples']-measured['delay_samples'],
                                            doppler_error_hz=prediction['doppler_hz']-measured['doppler_hz'],
                                            centre_complex_relative_error=float(abs(predicted-centre)/abs(centre)) if abs(centre) else None,
                                            centre_phase_error_deg=float(np.angle(predicted/centre)*180/np.pi) if abs(centre) else None,
                                            measured_snr_db=measured['snr_db'])
                                        if abs(centre):
                                            scores.append((prediction['explained_fraction'], prediction['assessment']['centre_complex_relative_error']))
                                    else:
                                        prediction['assessment'] = dict(no_held_candidate=True)
                                values.append(prediction)
                            row['predictions_json'] = np.frombuffer(json.dumps(values).encode(), np.uint8)
                            row['held_diagnostic_json'] = np.frombuffer(json.dumps(qa[index]).encode(), np.uint8)
                            report.flush()
                        print(path.name, f'fixed range {lo}-{hi}', 'queries', len(held), 'supported/matched', len(scores),
                              'median explained/error', np.median(scores, axis=0).tolist() if scores else [], flush=True)
        report.attrs['complete'] = True


if __name__ == '__main__':
    main()
