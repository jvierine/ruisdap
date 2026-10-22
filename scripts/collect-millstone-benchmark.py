#!/usr/bin/env python3
"""Freeze memo benchmark numerical evidence and exact source hashes into HDF5."""
import hashlib,json,pathlib,platform,re,subprocess
import h5py,numpy as np
root=pathlib.Path('/home/j/data/guisdap/20261019-lpi-benchmark-v1');repo=pathlib.Path('/home/j/src/guisdap26');output=root/'memo039-benchmark.h5';assert not output.exists()
with h5py.File(output,'w')as h:
 for name,file in [('stages','stages.h5'),('pipeline','pipeline.h5'),('reference_coded','reference-coded.h5'),('reference_long','reference-long.h5'),('cancellation','cancellation.h5')]:
  with h5py.File(root/file)as src:
   g=h.create_group(name)
   for key,value in src.attrs.items():g.attrs[key]=value
   for key in src:src.copy(key,g)
 prep=h.create_group('archive_preparation');log=(root/'archive-prepare.log').read_text();prep['wall_seconds']=47.25;prep['max_rss_kib']=int(re.search(r'Maximum resident set size \(kbytes\): (\d+)',log)[1]);prep.attrs['timing_scope']='one fresh process, three hours metadata discovery, ten seconds RF, concurrent campaign NFS readers; no cache flush'
 h['stages'].attrs['process_max_rss_kib']=888808
 for kind in ['coded','long']:
  f=root/f'archive-prepare/inputs/1712599210-ZENITH-{kind}.h5';p=pathlib.Path('/home/j/data/guisdap/20261019-native-campaign-pilot-v2/inputs')/f.name
  with h5py.File(f)as a,h5py.File(p)as b:
   for key in ['echo','transmit','pulse_sample','pair_index','science_pulse','detections']:np.testing.assert_array_equal(a[key][:],b[key][:])
   assert json.loads(bytes(a['configuration_json'][:]))==json.loads(bytes(b['configuration_json'][:]))
 fit=next((root/'fit').glob('*.h5'))
 with h5py.File(fit)as src:
  r=next(iter(src['records'].values()));rows=json.loads(bytes(r['ion_fit/results_json'][:]))
  counts={}
  for row in rows:
   result=row.get('result');status=result['status']if result else 'no_result';counts[status]=counts.get(status,0)+1
  h['pipeline'].attrs['fit_status_counts_json']=json.dumps(counts)
  h['pipeline'].attrs['fit_input_exposure_seconds']=10
  h['pipeline'].attrs['fit_configuration_json']=(root/'fit.json').read_text()
 h.attrs['host']=platform.node();h.attrs['cpu']='AMD EPYC 9754; 128 physical cores; SMT2; frequency boost enabled';h.attrs['rustc']=subprocess.check_output(['/home/j/.cargo/bin/rustc','--version'],text=True).strip();h.attrs['openmpi']='5.0.10-1 (Debian package)';h.attrs['benchmark_date_utc']='2026-10-19';h.attrs['thread_environment']='OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1'
 files=['crates/guisdap-millstone/src/'+n for n in ['lpi.rs','signal.rs','product.rs','debris.rs','fitting.rs']]+['crates/guisdap-millstone/examples/benchmark_lpi.rs','crates/guisdap-millstone/examples/benchmark_cancellation.rs','scripts/benchmark-millstone-lpi-reference.py','scripts/benchmark-millstone-pipeline.py']
 h.attrs['source_sha256_json']=json.dumps({f:hashlib.sha256((repo/f).read_bytes()).hexdigest()for f in files})
 h.attrs['source_git_base']=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip();h.attrs['source_state']='uncommitted working tree captured by source_sha256_json; frozen campaign kernel hash in pipeline metadata'
print(output)
print('fit status counts',counts)
