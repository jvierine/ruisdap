#!/usr/bin/env python3
"""One verified continuation after the original batch reports known guard errors.
Do not interrupt running workers. User stop or an unknown failure never resumes.
"""
import argparse,json,pathlib,subprocess,sys,time
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--root',action='append',required=True);p.add_argument('--unit',required=True);a=p.parse_args()
def properties():
 output=subprocess.check_output(['systemctl','--user','show',a.unit,'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True)
 return dict(row.split('=',1) for row in output.splitlines() if '=' in row)
while properties().get('ActiveState') in ['active','activating','deactivating']:time.sleep(15)
state=properties()
# Current old controller is expected to exit on its counted, already diagnosed
# per-window errors. Healthy MPI ranks must finish before this one continuation.
assert state.get('Result')=='exit-code' and state.get('ExecMainStatus')=='1',state
log=(a.output/'2024-04-08/reduce.log').read_text()
errors=[line for line in log.splitlines() if 'archive window ' in line]
assert errors and all('detection overlaps background window' in line or 'output exists; use a new product version' in line for line in errors),errors
print('Original batch finished; verifying completed hashes and resuming isolated window queues',flush=True)
subprocess.run([sys.executable,'scripts/repair-millstone-background-windows.py','--root',str(a.output)],check=True)
command=[sys.executable,'scripts/process-millstone-eclipse-direct.py','--output',str(a.output),'--resume']
for root in a.root:command+=['--root',root]
monitor=subprocess.Popen([sys.executable,'scripts/monitor-millstone-direct.py','--root',str(a.output),'--unit','ruisdap-millstone-direct-continuation-20261019.service'])
try:subprocess.run(command,check=True)
finally:monitor.terminate();monitor.wait()
