#!/usr/bin/env python3
"""Fit completed integrations during direct DigitalRF reduction; no RF staging.
All source windows must have atomic audits before fitting. Plans/checkpoints are
immutable. Public revisions are immutable; manifest replacement is the commit.
"""
import argparse, concurrent.futures, datetime as dt, fcntl, hashlib, json, os
import pathlib, subprocess, sys, tempfile, time, threading
import h5py
from isr_product import validate
PUBLIC = pathlib.Path('/srv/guisdap')
ANALYSIS = 'native-direct-live-20261019'
FIT_SERIES = 'legacy300-v3'
CPU_ASSIGNMENTS = {}
CPU_LOCK = threading.Lock()

def read_json(dataset):
    return json.loads(bytes(dataset[:]))

def atomic_json(path, value):
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as f:
        json.dump(value, f, allow_nan=False); tmp = pathlib.Path(f.name)
    tmp.chmod(0o644); tmp.replace(path)

def register(date, destination, analysis=ANALYSIS, label='Rust · live direct DigitalRF'):
    path = PUBLIC/'catalog/millstone-config.json'
    with path.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        config = json.loads(path.read_text())
        day = next(d for d in config['datasets'] if d['antenna']=='ZENITH' and d['date']==date)
        rows = day.setdefault('analyses', [])
        row = dict(id=analysis, label=label, data_root=destination.relative_to(PUBLIC).as_posix()+'/')
        rows[:] = [r for r in rows if r['id'] != analysis] + [row]
        atomic_json(path, config)

def partitions(date, day):
    """Legacy fit_lpi rule: anchor at first available ACF, keep starts < t0+300.
    A complete reduction watermark closes the integration, including empty and
    missing-RF audits. Unknown windows never masquerade as absent observations.
    The last integration is flushed when the day closes (legacy omitted it).
    """
    start = int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp())
    watermark = start
    available = []
    while watermark < start+86400:
        audit = day/'audit'/f'{watermark}.h5'
        if not audit.exists(): break
        with h5py.File(audit) as h:
            if read_json(h['outputs_json']): available.append(watermark)
        watermark += 10
    groups = []
    i = 0
    while i < len(available):
        anchor = available[i]
        end = min(anchor+300, start+86400)
        if watermark < end: break
        groups.append(list(range(anchor, end, 10)))
        i += 1
        while i < len(available) and available[i] < end: i += 1
    return groups

def plan_fit(day, slots, gain):
    audits = [day/'audit'/f'{t}.h5' for t in slots]
    if not all(p.exists() for p in audits): return None
    job = day/('stream-jobs-'+FIT_SERIES)/str(slots[0]); job.mkdir(parents=True, exist_ok=True)
    target = day/('stream-fits-'+FIT_SERIES)/f'fit-{slots[0]}.h5'
    if target.exists(): return None
    inputs = []
    for audit in audits:
        with h5py.File(audit) as h:
            for source in read_json(h['outputs_json']):
                path = pathlib.Path(source['path'])
                assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256'], 'Reduced source changed'
                inputs.append(str(path))
    if not inputs: return None
    inputs = sorted(inputs)
    m = dict(inputs=inputs,output=str(job/'fits'),integration_s=300., windows=[[slots[0],slots[-1]+10]],
        magic_constant_by_antenna={'ZENITH':gain},ions=[dict(mass_amu=16.,fraction=1.)],range_km=[90.,1075.],
        historical_settings=True,initial=dict(electron_density_m3=1e11,ion_temperature_k=1000.,electron_temperature_k=2000.,velocity_m_s=0.,ion_collision_s=0.),
        settings=dict(tie_electron_ion_temperature=False,fit_collision_rate=False,fixed_collision_s=0.,max_iterations=160,millstone_historical_bounds=True))
    path = job/'fit.json'
    if path.exists(): assert json.loads(path.read_text())==m, 'Completed integration changed'
    else: atomic_json(path,m)
    return path,target

def fit_one(binary, plan, target, cpus):
    # Bound worker threads to distinct existing physical cores. Concurrent live
    # and final controllers share a per-integration lease and immutable plans.
    with CPU_LOCK:
        cpu = CPU_ASSIGNMENTS.setdefault(threading.get_ident(), cpus[len(CPU_ASSIGNMENTS)%len(cpus)])
    with (plan.parent/'fit.lock').open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        if target.exists():return target
        env = {**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
        env.pop('DISPLAY',None)
        with (plan.parent/'fit.log').open('a') as log:
            subprocess.run(['taskset','-c',str(cpu),'nice','-n','5',str(binary),'fit',str(plan)],env=env,stdout=log,stderr=log,check=True)
        result = plan.parent/'fits/fit-000000.h5'
        with h5py.File(result) as h:
            validate(h)
            assert 'fit_checkpoint_sha256' in h
        os.link(result,target)
    return target

def _publish(root, date, context, binary_sha, complete=False):
    fits = root/date/('stream-fits-'+FIT_SERIES)
    files = sorted(fits.glob('fit-*.h5'))
    if not files: return
    destination = PUBLIC/'data'/f'millstone-{date}-zenith-native-direct-live-20261019-v1'
    destination.mkdir(exist_ok=True); destination.chmod(0o755)
    signature = hashlib.sha256(('\n'.join(p.name for p in files)+str(complete)+'debris-v2').encode()).hexdigest()
    old = json.loads((destination/'manifest.json').read_text()) if (destination/'manifest.json').exists() else {}
    if old.get('analysis_complete') and not complete:return
    if old.get('live_revision')==signature: return
    revisions = destination/'revisions'; revisions.mkdir(exist_ok=True); revisions.chmod(0o755)
    revision = revisions/signature
    if not revision.exists():
        staging = pathlib.Path(tempfile.mkdtemp(prefix='.export-',dir=revisions)); staging.rmdir()
        subprocess.run([sys.executable,'scripts/export-millstone-native-view.py','--historical',str(context),'--fits',str(fits),'--output',str(staging)],check=True)
        # Snapshot actual reductions/fit kernel provenance, never substitute a full-day claim.
        manifest = json.loads((staging/'manifest.json').read_text())
        start = int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp())
        manifest.update(analysis_id=ANALYSIS,analysis_label='Rust · live direct DigitalRF',live_revision=signature,
            analysis_live=not complete,analysis_complete=complete,analysis_default_window=[start,start+86400],
            kernel_sha256=binary_sha,input_mode='direct multi-directory DigitalRF',
            analysis_note='Live native fits of completed integration windows. Unprocessed windows are unavailable. Fits use finite-response amplitude-profiled starting densities to avoid the hot-ion local minimum. Integrations follow the legacy first-available-ACF 300 s maximum rule; antenna switching can shorten their exposure. Recorded transmit power and pulse response are normalized; initial calibration transfers the reviewed April 8 reference. Physical velocity sign remains unverified.')
        # Local asset URLs must resolve from the stable parent manifest.
        prefix='revisions/'+signature+'/'
        for key in ['parameters_url','debris_url']:
            manifest[key]=prefix+manifest[key]
        for f in manifest['fits']:
            for key in ['diagnostics_url','hdf5_url']: f[key]=prefix+f[key]
        atomic_json(staging/'manifest.json',manifest); staging.rename(revision)
    manifest=json.loads((revision/'manifest.json').read_text())
    atomic_json(destination/'manifest.json',manifest)
    register(date,destination)
    print(date,'LIVE',len(files),'complete fitted integrations',signature[:12],flush=True)

def publish(root, date, context, binary_sha, complete=False):
    destination=PUBLIC/'data'/f'millstone-{date}-zenith-native-direct-live-20261019-v1'
    destination.mkdir(exist_ok=True)
    with (destination/'.publish.lock').open('a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        return _publish(root,date,context,binary_sha,complete)

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--workers',type=int,default=8);p.add_argument('--once',action='store_true');p.add_argument('--drain',action='store_true');p.add_argument('--date');a=p.parse_args()
    dates=[a.date] if a.date else ['2024-04-08','2024-04-07','2024-04-09']
    assert 1<=a.workers<=128
    with h5py.File(a.root/'kernel.h5') as h:
        binary=a.root/str(h.attrs.get('fit_binary_name','guisdap-millstone'))
        binary_sha=h.attrs.get('fit_kernel_sha256',h.attrs['binary_sha256'])
    assert hashlib.sha256(binary.read_bytes()).hexdigest()==binary_sha
    # One hardware thread per physical core; all chosen cores already belong to reduction ranks.
    rows=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True).splitlines(); physical={}
    allowed=os.sched_getaffinity(0)
    for row in rows:
        if row.startswith('#'):continue
        cpu,core,socket=map(int,row.split(','))
        if cpu in allowed:physical.setdefault((socket,core),cpu)
    cpus=list(physical.values())[:128]; assert cpus
    pending={}; scheduled=set(); failures={}
    reference=PUBLIC/'data/millstone-2024-04-08-zenith-cal-d9ac4d6a1088'
    ref=json.loads((reference/'manifest.json').read_text());gain=ref['calibration']['gain']*480*(299792.458*30/(2*1e6))
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
        while True:
            for future in list(pending):
                if not future.done():continue
                key=pending.pop(future)
                try: print('FIT READY',future.result(),flush=True)
                except Exception as error:
                    failures[key]=str(error);print('FIT FAILED',key,str(error),flush=True)
            for date in dates:
                day=a.root/date
                if not day.exists():continue
                (day/('stream-fits-'+FIT_SERIES)).mkdir(exist_ok=True)
                context=reference if date=='2024-04-08' else PUBLIC/'data'/f'millstone-{date}-zenith-campaign-context-v1'
                for slots in partitions(date,day):
                    if len(pending)>=a.workers:break
                    key=(date,slots[0])
                    if key in scheduled or key in failures:continue
                    planned=plan_fit(day,slots,gain)
                    if planned:
                        pending[pool.submit(fit_one,binary,*planned,cpus)]=key;scheduled.add(key)
                publish(a.root,date,context,binary_sha)
            atomic_json(a.root/'live-stream.json',dict(pending=len(pending),failed={str(k):v for k,v in failures.items()},updated_unix_s=time.time()))
            state=a.root/'live/rank-0.json'
            if not pending and state.exists() and json.loads(state.read_text()).get('phase','').startswith('all available ZENITH eclipse days complete'):return
            if a.drain and not pending:
                assert not failures, failures
                return
            if a.once:
                for f in pending:f.result()
                for date in dates:
                    context=reference if date=='2024-04-08' else PUBLIC/'data'/f'millstone-{date}-zenith-campaign-context-v1'
                    if (a.root/date/('stream-fits-'+FIT_SERIES)).exists():publish(a.root,date,context,binary_sha)
                return
            time.sleep(15)

if __name__=='__main__':main()
