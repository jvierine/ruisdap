#!/usr/bin/env python3
"""128-core native campaign controller: original multi-root DigitalRF, no staging.
Reduce to atomic HDF5 ACFs, fit completed legacy 300 s windows, audit, publish.
The independent live scheduler can fit the same immutable jobs during reduction.
"""
import argparse, datetime as dt, fcntl, hashlib, importlib.util, json, os
import pathlib, shutil, subprocess, sys, tempfile, time
import h5py
import numpy as np
from isr_product import validate

PUBLIC=pathlib.Path('/srv/guisdap')
DATES=['2024-04-08','2024-04-07','2024-04-09']
REFERENCE=PUBLIC/'data/millstone-2024-04-08-zenith-cal-d9ac4d6a1088'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def atomic_json(path,value):
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False) as h:
        json.dump(value,h,allow_nan=False);tmp=pathlib.Path(h.name)
    tmp.chmod(0o644);tmp.replace(path)

def bytes_json(h,name,value):h[name]=np.frombuffer(json.dumps(value).encode(),np.uint8)

def audit_day(day,start,kernels):
    science=[];candidate=[];gaps=[];products=[];rows=[];rejected=[]
    original=kernels['binary_sha256']
    for begin in range(start,start+86400,10):
        with h5py.File(day/'audit'/f'{begin}.h5') as h:
            science.append(h['science_pulse_sample'][:]);candidate.append(h['candidate_science_pulse_sample'][:])
            gaps.extend(json.loads(bytes(h['gaps_json'][:])))
            outputs=json.loads(bytes(h['outputs_json'][:]))
        for row in outputs:
            path=pathlib.Path(row['path']);assert digest(path)==row['sha256'],'Reduced input changed'
            with h5py.File(path) as h:
                validate(h)
                for record in h['records'].values():
                    if 'debris/fully_rejected_pulse_sample' in record:
                        config=json.loads(bytes(record['configuration_json'][:]))
                        selected=record['pulse_sample'][:][np.array(config['lpi']['science_pulse'],dtype=bool)]
                        rejected.extend(np.intersect1d(selected,record['debris/fully_rejected_pulse_sample'][:]).tolist())
            kernel=row.get('kernel_sha256') or original
            assert kernel in set(kernels.values()),'Unknown reduction kernel'
            rows.append(dict(path=str(path),sha256=row['sha256'],kernel_sha256=kernel))
            products.append(str(path))
    science=np.concatenate(science);candidate=np.concatenate(candidate)
    assert len(np.unique(science))==len(science),'Duplicate science pulse'
    assert len(np.unique(candidate))==len(candidate),'Duplicate candidate pulse'
    assert np.isin(science,candidate).all(),'Science pulse missing from candidate audit'
    assert products and len(set(products))==len(products),'Missing or duplicated reduced products'
    path=day/'coverage.h5'
    if not path.exists():
        tmp=path.with_suffix('.tmp.h5')
        with h5py.File(tmp,'w') as h:
            h['science_pulse_sample']=science;h['candidate_science_pulse_sample']=candidate
            h['excluded_candidate_sample']=np.setdiff1d(candidate,science)
            h['fully_rejected_science_pulse_sample']=np.unique(np.asarray(rejected,dtype=np.uint64))
            bytes_json(h,'gaps_json',gaps);bytes_json(h,'products_json',rows)
            h.attrs['fit_kernel_sha256']=kernels.get('fit_kernel_sha256',original)
            h.attrs['input_mode']='original multi-directory DigitalRF; no RF intermediate'
            h.attrs['date']=day.name
        tmp.replace(path)
    else:
        with h5py.File(path) as h:
            assert np.array_equal(h['science_pulse_sample'][:],science),'Coverage changed'
            assert json.loads(bytes(h['products_json'][:]))==rows,'Reduced product ledger changed'
    return products,science,candidate,rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,action='append',required=True)
    p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    assert a.resume or not a.output.exists(),'Use a fresh version or explicitly resume'
    a.output.mkdir(parents=True,exist_ok=a.resume);(a.output/'live').mkdir(exist_ok=True)
    binary=a.output/'guisdap-millstone';kernel=a.output/'kernel.h5'
    if not a.resume:
        shutil.copy2('target/release/guisdap-millstone',binary);binary.chmod(0o700)
        with h5py.File(kernel,'w') as h:
            h.attrs['binary_sha256']=digest(binary);h.attrs['input_mode']='direct multi-directory DigitalRF; no staging'
            bytes_json(h,'roots_json',list(map(str,a.root)))
    with h5py.File(kernel) as h:
        kernels={k:str(v) for k,v in h.attrs.items() if k.endswith('_sha256')}
        assert digest(binary)==kernels['binary_sha256'],'Original frozen kernel changed'
        binary=a.output/str(h.attrs.get('fit_binary_name',binary.name))
        fit_sha=str(h.attrs.get('fit_kernel_sha256',kernels['binary_sha256']))
        assert digest(binary)==fit_sha,'Frozen fit kernel changed'
        assert json.loads(bytes(h['roots_json'][:]))==list(map(str,a.root)),'Archive roots changed'
        reducer=a.output/str(h.attrs.get('reduction_binary_name',binary.name))
        expected=str(h.attrs.get('reduction_kernel_sha256',kernels['binary_sha256']))
        assert digest(reducer)==expected,'Frozen reduction kernel changed'
    spec=importlib.util.spec_from_file_location('stream','scripts/stream-millstone-direct.py')
    stream=importlib.util.module_from_spec(spec);spec.loader.exec_module(stream)
    env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'};env.pop('DISPLAY',None)
    mpi=['/usr/bin/mpirun','-np','128','--map-by','core','--bind-to','core']
    def status(phase,completed=0,total=1):
        atomic_json(a.output/'live/rank-0.json',dict(host='avaruus',rank=0,ranks=1,phase=phase,completed=completed,total=total,updated_unix_s=time.time()))
        print(phase,flush=True)
    def run(command,log):
        with log.open('a') as h:subprocess.run(command,env=env,stdout=h,stderr=h,check=True)
    for date in DATES:
        day=a.output/date;day.mkdir(exist_ok=True);(day/'live').mkdir(exist_ok=True)
        destination=PUBLIC/'data'/f'millstone-{date}-zenith-native-direct-20261019-v1'
        if (destination/'manifest.json').exists():
            manifest=json.loads((destination/'manifest.json').read_text())
            assert manifest['kernel_sha256']==fit_sha and manifest['analysis_complete']
            stream.register(date,destination,'native-direct-20261019','Rust · direct DigitalRF · 45 lags')
            continue
        start=int(dt.datetime.fromisoformat(date).replace(tzinfo=dt.timezone.utc).timestamp())
        plan=day/'archive.json';m=dict(roots=list(map(str,a.root)),start_unix_s=start,end_unix_s=start+86400,output=str(day),live_directory=str(day/'live'))
        if plan.exists():assert json.loads(plan.read_text())==m,'Archive plan changed'
        else:atomic_json(plan,m)
        lock=pathlib.Path('/home/j/data/guisdap/2022-11-02/calibration/native-compute.lock');lock.parent.mkdir(parents=True,exist_ok=True)
        with lock.open('a') as lease:
            fcntl.flock(lease,fcntl.LOCK_EX)
            status(date+' direct DigitalRF reduction on 128 cores')
            run([*mpi,str(reducer),'process-archive',str(plan)],day/'reduce.log')
            status(date+' verifying pulse coverage, gaps and reduced products')
            products,science,candidate,rows=audit_day(day,start,kernels)
            status(date+' fitting completed legacy 300 s integration windows on 128 cores')
            run([sys.executable,'scripts/stream-millstone-direct.py','--root',str(a.output),'--date',date,'--workers','128','--drain'],day/'fit.log')
            fits=day/('stream-fits-'+stream.FIT_SERIES);fit_files=sorted(fits.glob('fit-*.h5'));used=[]
            assert fit_files,'No complete native fits'
            for path in fit_files:
                with h5py.File(path) as h:
                    validate(h);assert 'fit_checkpoint_sha256' in h
                    sources=json.loads(bytes(h['sources/combined/records_json'][:]))
                    for row in sources:
                        assert digest(pathlib.Path(row['path']))==row['sha256'],'Fitted input changed'
                        used.append(row['path'])
            assert len(used)==len(set(used)) and set(used)==set(products),'Fit coverage differs from reduced products'
            status(date+' publishing verified complete day')
            context=REFERENCE if date=='2024-04-08' else PUBLIC/'data'/f'millstone-{date}-zenith-campaign-context-v1'
            staging=destination.with_name('.'+destination.name)
            assert not staging.exists(),'Inspect interrupted export staging before resuming'
            run([sys.executable,'scripts/export-millstone-native-view.py','--historical',str(context),'--fits',str(fits),'--output',str(staging)],day/'export.log')
            manifest=json.loads((staging/'manifest.json').read_text())
            with h5py.File(day/'coverage.h5') as coverage: rejected_count=len(coverage['fully_rejected_science_pulse_sample'])
            manifest.update(analysis_id='native-direct-20261019',analysis_live=False,analysis_complete=True,
                analysis_label='Rust · direct DigitalRF · 45 lags',analysis_default_window=[start,start+86400],
                science_pulse_count=len(science),fully_rejected_science_pulse_count=rejected_count,candidate_science_pulse_count=len(candidate),excluded_candidate_count=len(candidate)-len(science),
                kernel_sha256=fit_sha,reduction_kernel_sha256=sorted(set(r['kernel_sha256'] for r in rows)),
                input_mode='direct multi-directory DigitalRF',integration_rule='Legacy first-available ACF anchor; source starts strictly below anchor + 300 s; trailing integration retained',
                analysis_note='Native DigitalRF analysis of all supported powered monostatic ZENITH coded and long pulses available on this day. 45 measured lags and historical range averaging; legacy 300 s maximum integration rule with shorter antenna-selected exposure. Missing RF is audited. Experimental CFAR-seeded amplitude subtraction; background-contaminated pulses are explicitly rejected. Finite-response, finite-Debye Rust fits; initial calibration transfers the reviewed April 8 reference. Physical velocity sign remains unverified.')
            atomic_json(staging/'manifest.json',manifest)
            shutil.copy2(day/'coverage.h5',staging/'coverage.h5');(staging/'coverage.h5').chmod(0o644)
            staging.rename(destination)
            stream.register(date,destination,'native-direct-20261019','Rust · direct DigitalRF · 45 lags')
            stream.publish(a.output,date,context,fit_sha,complete=True)
            print(date,'PUBLISHED',len(fit_files),'fits',len(science),'science pulses',flush=True)
    status('all available ZENITH eclipse days complete',3,3)

if __name__=='__main__':main()
