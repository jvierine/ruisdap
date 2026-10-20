#!/usr/bin/env python3
"""Explore inter-pulse phase coherence from saved complex-amplitude diagnostics.

Metadata only; no RF reads, no change to production selection. Random-phase
controls condition on the selected pulse times and the same frequency/chirp
search. They do not establish a campaign false-alarm probability.
"""
import argparse
import hashlib
import json
import pathlib

import h5py
import numpy as np


def phase_search(times, amplitudes, frequencies, centres, null_trials, rng):
    # Fitted amplitudes refer to pulse-local zero. Refer them to the actual
    # echo centre to avoid phase jumps from the noisy per-pulse Doppler estimate.
    phases = amplitudes*np.exp(2j*np.pi*frequencies*centres/1e6)
    phases /= abs(phases)
    t = times+centres/1e6
    t -= t.mean()
    span = np.ptp(t)
    assert span > 0
    slope = float(np.polyfit(t, frequencies, 1)[0])
    f0 = float(np.median(frequencies))
    grid = np.arange(f0-1000, f0+1000+1/(4*span), 1/(4*span))
    # Fixed acceleration search is also applied to every random-phase control.
    chirps = np.linspace(-1000, 1000, 81)
    matrix = np.concatenate([np.exp(-2j*np.pi*(grid[:, None]*t + .5*chirp*t*t)) for chirp in chirps])
    matched = matrix@phases
    best = int(np.argmax(abs(matched)))
    score = abs(matched[best])**2/len(times)
    null = np.exp(2j*np.pi*rng.random((len(times), null_trials)))
    null_scores = np.max(abs(matrix@null)**2/len(times), axis=0)
    return dict(score=float(score), normalized_coherence=float(score/len(times)),
                best_frequency_hz=float(grid[best % len(grid)]), best_chirp_hz_s=float(chirps[best//len(grid)]),
                diagnostic_doppler_slope_hz_s=slope, random_phase_p=(1+np.count_nonzero(null_scores >= score))/(1+null_trials),
                search_cells=len(matrix), pulse_count=len(times)), null_scores


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--product', type=pathlib.Path, action='append', required=True)
    p.add_argument('--output', type=pathlib.Path, required=True)
    p.add_argument('--null-trials', type=int, default=64)
    p.add_argument('--anchors', type=int, default=4)
    a = p.parse_args()
    assert not a.output.exists() and a.null_trials > 0 and a.anchors > 0
    rng = np.random.default_rng(397222)
    with h5py.File(a.output, 'x') as report:
        report.attrs['complete'] = False
        report.attrs['scope'] = __doc__
        report.attrs['script_sha256'] = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
        for path in a.product:
            with h5py.File(path) as source:
                for record_id, record in source['records'].items():
                    c = json.loads(bytes(record['configuration_json'][:]))
                    qa = json.loads(bytes(record['debris/diagnostics_json'][:]))
                    keys = record['pulse_sample'][:]
                    events = []
                    for i, q in enumerate(qa):
                        if not c['lpi']['science_pulse'][i] or not q['amplitudes']:
                            continue
                        for j, (v, amplitude) in enumerate(zip(q['detections'], q['amplitudes'])):
                            if q['amplitude_snr_db'][j] < 14:
                                continue
                            events.append((int(keys[i]), v['delay_samples'], v['doppler_hz'], complex(*amplitude)))
                    out = report.create_group(path.stem+'/'+record_id)
                    out.attrs['source_path'] = str(path)
                    out.attrs['source_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
                    if not events:
                        continue
                    samples = np.array([r[0] for r in events], np.int64)
                    times = (samples-samples.min())/1e6
                    delays = np.array([r[1] for r in events])
                    frequencies = np.array([r[2] for r in events])
                    amplitudes = np.array([r[3] for r in events])
                    centres = delays+np.mean(c['lpi']['tx_gate'])
                    for lo, hi in [(250., 500.), (525., 575.)]:
                        eligible = np.flatnonzero((delays*299792.458/2e6 >= lo) & (delays*299792.458/2e6 <= hi))
                        if not len(eligible):
                            continue
                        group = out.create_group(f'range_{lo}_{hi}')
                        for number, index in enumerate(np.unique(eligible[np.linspace(0, len(eligible)-1, min(a.anchors, len(eligible))).astype(int)])):
                            dt = times-times[index]
                            shift = frequencies[index]/440.2e6*1e6*dt
                            distance = np.minimum(abs(delays-delays[index]-shift), abs(delays-delays[index]+shift))
                            use = np.flatnonzero((abs(dt) <= .25) & (distance <= 16) & (abs(frequencies-frequencies[index]) <= 2000))
                            # Use at most one component per pulse; prefer the
                            # closest range to avoid counting one pulse twice.
                            chosen = {}
                            for j in use:
                                key = int(samples[j])
                                if key not in chosen or distance[j] < distance[chosen[key]]:
                                    chosen[key] = j
                            use = np.array(sorted(chosen.values(), key=lambda j: abs(dt[j]))[:48])
                            if len(use) < 16:
                                continue
                            result, null = phase_search(times[use], amplitudes[use], frequencies[use], centres[use], a.null_trials, rng)
                            row = group.create_group(str(number))
                            row.attrs.update(result)
                            row['pulse_sample'] = samples[use]
                            row['null_max_scores'] = null
                            print(path.name, f'range={lo}-{hi}', f'anchor={samples[index]}', json.dumps(result), flush=True)
        report.attrs['complete'] = True


if __name__ == '__main__':
    main()
