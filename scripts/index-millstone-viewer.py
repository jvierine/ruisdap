#!/usr/bin/env python
"""Index immutable native products for on-demand spectrum and ACF viewing."""
import argparse,json
from pathlib import Path
import h5py

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path)
p.add_argument('quicklook',type=Path)
p.add_argument('output',type=Path)
p.add_argument('--date',required=True)
a=p.parse_args()
index={'plasma':{},'datasets':{}}
pending={}
for antenna,channel in [('ZENITH','zenith-l'),('MISA','misa-l')]:
    manifest=a.quicklook/antenna.lower()/'manifest.json'
    if not manifest.is_file():continue
    m=json.loads(manifest.read_text())
    if m['date']!=a.date or m['antenna']!=antenna:raise ValueError('Quicklook identity mismatch')
    rows=[]
    for time in m['raw_times']:
        path=a.source/'lpi_30'/channel/f'lpi-{int(time)}.h5'
        if not path.is_file():raise ValueError(f'Missing native ACF: {path}')
        rows.append({'time':time,'path':str(path.resolve())})
    index['datasets'][a.date+'/'+antenna]={'raw':rows,'manifest':str(manifest.resolve())}
    for record in m['plasma_records']:pending[(antenna,record['source_key'])]=record
for directory in sorted((a.source/'metadata/integrated_plasma_line_metadata_hires').glob(a.date+'T*')):
    for path in sorted(directory.glob('*.h5')):
        with h5py.File(path)as h:
            for key,g in h.items():
                if not isinstance(g,h5py.Group)or'antenna'not in g:continue
                antenna=g['antenna'][()]
                if isinstance(antenna,bytes):antenna=antenna.decode()
                identity=(str(antenna),key)
                if identity not in pending:continue
                record=pending.pop(identity)
                index['plasma'][record['id']]={'path':str(path.resolve()),'key':key,'antenna':str(antenna),'date':a.date,'start':record['start'],'end':record['end']}
if pending:raise ValueError(f'{len(pending)} plasma records have no native source')
a.output.parent.mkdir(parents=True,exist_ok=True)
temporary=a.output.with_suffix('.json.new')
temporary.write_text(json.dumps(index));temporary.chmod(0o600);temporary.replace(a.output)
print(f"Indexed {len(index['plasma'])} plasma records and {sum(len(d['raw'])for d in index['datasets'].values())} ACF files")
