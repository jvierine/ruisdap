#!/usr/bin/env python3
"""Independently verify eligible pulses using upstream DigitalRF/Metadata.

Reads original metadata and RF block indices only. Writes an HDF5 coverage
ledger, never RF, converted samples, or a combined recording directory.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import time

import digital_rf
import h5py
import numpy as np


def canonical(value):
    if isinstance(value, dict):
        return {k: canonical(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        return canonical(value.item()) if value.size == 1 else value.tolist()
    if isinstance(value, np.generic):
        return canonical(value.item())
    return value.decode() if isinstance(value, bytes) else value


def read_metadata(readers, start, end, columns, pad=False):
    records = {}
    for reader in readers:
        for key, value in reader.read(max(0, start), end, columns=columns,
                                      method='pad' if pad else None).items():
            value = canonical(value)
            if key in records and records[key] != value:
                raise ValueError(f'Conflicting original metadata at {key}')
            records[int(key)] = value
    return dict(sorted(records.items()))


def covered(keys, blocks, length=10000):
    """Whether each complete requested vector lies in one continuous block."""
    if not blocks:
        return np.zeros(len(keys), dtype=bool)
    starts = np.array(list(blocks), dtype=np.uint64)
    ends = starts + np.array(list(blocks.values()), dtype=np.uint64)
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]
    index = np.searchsorted(starts, keys, side='right') - 1
    return (index >= 0) & (keys + length <= ends[np.maximum(index, 0)])


def paired_science(keys, codes, begin, end):
    """Nearest adjacent same-code partner <40 ms; >=16 pulses per family/window."""
    if not len(keys):
        return np.empty(0, dtype=np.uint64)
    paired = np.zeros(len(keys), dtype=bool)
    close = (np.diff(keys) < 40000) & (codes[1:] == codes[:-1])
    paired[:-1] |= close
    paired[1:] |= close
    selected = keys[paired & (keys >= begin) & (keys < end)]
    return selected if len(selected) >= 16 else np.empty(0, dtype=np.uint64)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=pathlib.Path, action='append', required=True)
    p.add_argument('--campaign', type=pathlib.Path, required=True)
    p.add_argument('--date', required=True)
    p.add_argument('--seconds', type=int, default=86400)
    p.add_argument('--report', type=pathlib.Path, required=True)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    assert 0 < args.seconds <= 86400 and args.seconds % 10 == 0
    assert args.resume or not args.report.exists(), 'Choose a fresh report or explicitly resume it'
    start = int(dt.datetime.fromisoformat(args.date).replace(tzinfo=dt.timezone.utc).timestamp())
    day = args.campaign / args.date
    streams = {}
    for name, rate in [('id_metadata', 1000000), ('antenna_control_metadata', 1000000), ('powermeter', 1)]:
        readers = []
        for root in args.root:
            path = root / 'metadata' / name
            if not path.is_dir():
                continue
            with h5py.File(path / 'dmd_properties.h5') as props:
                assert int(props.attrs['sample_rate_numerator']) == rate
                assert int(props.attrs['sample_rate_denominator']) == 1
            readers.append(digital_rf.DigitalMetadataReader(str(path)))
        assert readers, 'No original metadata: ' + name
        streams[name] = readers
    rf = digital_rf.DigitalRFReader(
        [str(r / 'rf_data') for r in args.root], rdcc_nbytes=64 * 1024 * 1024,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    began = time.monotonic()
    total_candidate = total_science = 0
    with h5py.File(args.report, 'a' if args.resume else 'x') as report:
        if 'roots_json' in report:
            assert json.loads(bytes(report['roots_json'][:])) == list(map(str, args.root))
            assert report.attrs['date'] == args.date
            assert report.attrs['requested_seconds'] == args.seconds
        else:
            report['roots_json'] = np.frombuffer(json.dumps(list(map(str, args.root))).encode(), np.uint8)
        report.attrs['requested_seconds'] = args.seconds
        report.attrs['reference'] = 'digital_rf ' + digital_rf.__version__
        report.attrs['date'] = args.date
        report.attrs['complete_utc_day'] = args.seconds == 86400
        report.attrs['input_mode'] = 'Original DigitalMetadata and RF continuous-block indices; no RF samples read or written'
        report.attrs['selection'] = 'Powered >=400 kW monostatic ZENITH; codes 1..32 or 300; complete 10000-sample RX/TX vectors; real same-code adjacent partner <40 ms; >=16 paired pulses per family in each ten-second window'
        report.attrs['verified'] = False
        for lo in range(start, start + args.seconds, 300):
            hi = min(start + args.seconds, lo + 300)
            windows = list(range(lo, hi, 10))
            if all(str(w) in report and report[str(w)].attrs.get('verified', False) for w in windows):
                for window in windows:
                    g = report[str(window)]
                    assert g.attrs['audit_sha256'] == hashlib.sha256((day/'audit'/f'{window}.h5').read_bytes()).hexdigest(), 'Previously verified audit changed'
                    total_candidate += len(g['expected_candidate_sample'])
                    total_science += len(g['expected_science_sample'])
                continue
            ids = read_metadata(streams['id_metadata'], (lo-1)*1000000,
                                (hi+1)*1000000-1, 'sweepid')
            ants = read_metadata(streams['antenna_control_metadata'], (lo-1)*1000000,
                                 (hi+1)*1000000-1, ['rx_antenna', 'tx_antenna'], pad=True)
            meter = read_metadata(streams['powermeter'], lo-31, hi+41, 'zenith_power')
            keys = np.array(list(ids), dtype=np.uint64)
            codes = np.array(list(ids.values()), dtype=np.uint64)
            antkeys = np.array(list(ants), dtype=np.uint64)
            state = np.array([v['rx_antenna'] == v['tx_antenna'] == 'ZENITH' for v in ants.values()])
            index = np.searchsorted(antkeys, keys, side='right') - 1
            zenith = ((index >= 0) & state[np.maximum(index, 0)]) if len(state) else np.zeros(len(keys), bool)
            pt = np.array(list(meter), dtype=float)
            pv = np.array(list(meter.values()), dtype=float)
            seconds = keys.astype(float) / 1e6
            power = np.interp(seconds, pt, pv) if len(pt) else np.full(len(keys), np.nan)
            if len(pt):
                pos = np.searchsorted(pt, seconds)
                distance = np.minimum(abs(seconds-pt[np.maximum(pos-1, 0)]),
                                      abs(seconds-pt[np.minimum(pos, len(pt)-1)]))
                power[distance > 30] = np.nan
            eligible = zenith & np.isfinite(power) & (power >= 400000) & (((codes >= 1) & (codes <= 32)) | (codes == 300))
            keys, codes = keys[eligible], codes[eligible]
            available = np.ones(len(keys), bool)
            if len(keys):
                for channel in ['zenith-l', 'tx-h']:
                    blocks = rf.get_continuous_blocks(int(keys[0]), int(keys[-1])+9999, channel)
                    available &= covered(keys, blocks)
            for window in range(lo, hi, 10):
                begin, end = window*1000000, (window+10)*1000000
                expected_candidate = keys[(keys >= begin) & (keys < end)]
                expected_science = []
                for family in [codes <= 32, codes == 300]:
                    use = available & family
                    expected_science.extend(paired_science(keys[use], codes[use], begin, end))
                expected_science = np.sort(np.array(expected_science, np.uint64))
                path = day / 'audit' / f'{window}.h5'
                with h5py.File(path) as audit:
                    observed_candidate = np.sort(audit['candidate_science_pulse_sample'][:])
                    observed_science = np.sort(audit['science_pulse_sample'][:])
                if str(window) in report:
                    existing = report[str(window)]
                    if existing.attrs.get('verified', False):
                        assert existing.attrs['audit_sha256'] == hashlib.sha256(path.read_bytes()).hexdigest(), 'Previously verified audit changed'
                    else:
                        del report[str(window)]
                group = report.require_group(str(window))
                for name, expected in [('expected_candidate_sample', expected_candidate), ('expected_science_sample', expected_science)]:
                    if name in group:
                        np.testing.assert_array_equal(group[name][:], expected)
                    else:
                        group[name] = expected
                group.attrs['audit_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
                np.testing.assert_array_equal(observed_candidate, expected_candidate, err_msg=f'{window} candidate pulse coverage differs')
                np.testing.assert_array_equal(observed_science, expected_science, err_msg=f'{window} selected pulse coverage differs')
                group.attrs['verified'] = True
                total_candidate += len(expected_candidate)
                total_science += len(expected_science)
            report.flush()
            print(f'{args.date} verified through {hi-start}s: {total_science} science / {total_candidate} candidates; elapsed {time.monotonic()-began:.1f}s', flush=True)
        report.attrs['candidate_pulse_count'] = total_candidate
        report.attrs['science_pulse_count'] = total_science
        report.attrs['verified'] = True
    print('INDEPENDENT COVERAGE VERIFIED', args.report, flush=True)


if __name__ == '__main__':
    main()
