#!/usr/bin/env python3
"""Run one bounded, low-priority, single-threaded Magma job over SSH.

Example: python3 run_remote.py --host YOUR_HOST --case 5a --seconds 600
Remote files stay in ~/coates_sujatha_20261002. A lock forbids overlapping jobs
launched through this runner. No remote job survives its wall-time limit.
"""
import argparse
import datetime
import pathlib
import subprocess
import sys

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--host', required=True)
p.add_argument('--case', choices=['5a','5b','5c','5d','7a','7b'], required=True)
p.add_argument('--seconds', type=int, default=600)
p.add_argument('--effort', type=int, default=2)
p.add_argument('--grh', action='store_true')
p.add_argument('--diagnostic', action='store_true')
p.add_argument('--local-hints', action='store_true')
p.add_argument('--selmer-only', action='store_true')
p.add_argument('--twist', type=int, default=1)
args = p.parse_args()
assert 1 <= args.seconds <= 3600 and 1 <= args.effort <= 100
root = pathlib.Path(__file__).resolve().parent
ssh = ['ssh','-o','BatchMode=yes','-o','ClearAllForwardings=yes',
       '-o','ConnectTimeout=20',args.host]
remote = 'coates_sujatha_20261002'
source = (root/'magma/descent.m').read_bytes()
subprocess.run(ssh+[f'mkdir -p {remote} && cat > {remote}/descent.m'],
               input=source,check=True)
job = f'{args.case}_tw{args.twist}_grh{int(args.grh)}_e{args.effort}_{args.seconds}s'
if args.diagnostic:
    job += '_diagnostic'
if args.local_hints:
    job += '_hints'
if args.selmer_only:
    job += '_selmer'
cmd = (f'cd {remote} && '
       'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 '
       f'flock -n job.lock nice -n 10 timeout --signal=TERM --kill-after=10 {args.seconds} '
       f'magma -n -t 1 -l 4294967296 -b case_id:={args.case} '
       f'grh:={int(args.grh)} effort:={args.effort} twist:={args.twist} '
       f'diagnostic:={int(args.diagnostic)} local_hints:={int(args.local_hints)} '
       f'selmer_only:={int(args.selmer_only)} descent.m')
out = root/'data'/f'{job}.out'
with out.open('w') as f:
    f.write(f'# UTC {datetime.datetime.now(datetime.timezone.utc).isoformat()}\n')
    f.write('# One thread, nice 10, memory cap 4 GiB.\n')
    f.write(f'# Command: {cmd}\n')
    f.flush()
    proc = subprocess.Popen(ssh+[cmd],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    for line in proc.stdout:
        s = line.decode(errors='replace')
        sys.stdout.write(s); sys.stdout.flush()
        f.write(s); f.flush()
    rc = proc.wait()
    f.write(f'\n# EXIT {rc}\n')
print(f'LOG {out}\nEXIT {rc}',flush=True)
sys.exit(rc)
