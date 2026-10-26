#!/usr/bin/env python3
"""Add (or replace, by id) one analysis in the Millstone viewer's selector.

The viewer lists, per date and antenna, the historical data and the
analyses in catalog/millstone-config.json. This adds
{id, label, data_root} to one dataset's "analyses", holding an exclusive lock
on millstone-config.lock and replacing the file atomically, so concurrent
publishers do not lose each other's entries. Nothing else in the file
changes. --remove takes an entry out again.

Usage:
    python3 scripts/register-millstone-analysis.py CONFIG --date 2024-04-08 --antenna ZENITH \
        --id isr-tracked-v1 --label "isr_analysis · tracked satellite columns" \
        --data-root data/millstone-2024-04-08-zenith-isr-tracked-v1/
    python3 scripts/register-millstone-analysis.py CONFIG --date ... --antenna ... --id ... --remove
"""
import argparse, fcntl, json, os
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument('config', type=Path, help='catalog/millstone-config.json')
p.add_argument('--date', required=True)
p.add_argument('--antenna', required=True)
p.add_argument('--id', required=True)
p.add_argument('--label')
p.add_argument('--data-root')
p.add_argument('--remove', action='store_true')
a = p.parse_args()
if not a.remove and not (a.label and a.data_root):
    raise SystemExit('--label and --data-root are required unless --remove')
if a.id == 'historical':
    raise SystemExit('"historical" is the dataset itself, not an analysis')

with open(a.config.with_suffix('.lock'), 'a') as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    c = json.loads(a.config.read_text())
    d = next((d for d in c['datasets'] if d['date'] == a.date and d['antenna'] == a.antenna), None)
    if d is None:
        raise SystemExit('no dataset %s %s in %s' % (a.date, a.antenna, a.config))
    choices = d.setdefault('analyses', [])
    before = len(choices)
    choices[:] = [x for x in choices if x['id'] != a.id]
    if not a.remove:
        choices.append(dict(id=a.id, label=a.label, data_root=a.data_root))
    tmp = a.config.with_suffix('.register.tmp')
    tmp.write_text(json.dumps(c))
    os.chmod(tmp, 0o644)
    tmp.replace(a.config)
print('%s %s %s: %d -> %d analyses' % (a.date, a.antenna, a.id, before, len(choices)))
