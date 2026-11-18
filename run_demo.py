import argparse
import json
import os
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--mode',choices=['run','candidate_replace'],default='run')
parser.add_argument('--output',type=Path,default=ROOT/'results/normal')
args=parser.parse_args()
output=args.output.resolve()
output.mkdir(parents=True,exist_ok=False)
with np.load(ROOT/'data/resources.npz') as z:
    with (output/'inputs.f64').open('xb') as f:
        z['q'].tofile(f)
        z['bounds'].tofile(f)
command=[str(ROOT/'build/lobule_3000'),str(output/'inputs.f64'),str(output/'demo'),
         str(ROOT/'data/oxygen.txt'),str(ROOT/'data/hosts.txt'),str(ROOT/'data/host.xml'),
         str(ROOT/'data/fields.xml'),args.mode]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
with (output/'console.log').open('x') as log:
    r=subprocess.run(command,cwd=output,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
expected=0 if args.mode=='run' else 42
assert r.returncode==expected, f'Unexpected exit {r.returncode}; inspect {output}/console.log'
report=json.loads((output/'demo_coordinator.json').read_text())
assert report['canonical_parent_unchanged']
assert report['accepted_snapshot']==(args.mode=='run')
if args.mode=='run':
    child=json.loads((output/'demo_child.json').read_text())
    assert child['agents']==child['individual_models']==3000
    assert child['joint_LP_calls']==1884 and child['scheduler_entries']==1
    assert child['clock_min']==.0075/60 and child['commits']==1
    payload=(output/'demo_certified.f64').read_bytes()
    assert payload==(output/'demo_candidate.f64').read_bytes()==(output/'demo_committed.f64').read_bytes()
    a=np.frombuffer(payload,dtype=np.float64)
    assert a.shape==(589600,) and np.isfinite(a).all()
    q=a[:50000].reshape(2500,20)
    assert (q>=0).all()
else:
    assert report['rejection_reason']=='certified_payload_mismatch'
    assert (output/'demo_candidate.f64').read_bytes()==(output/'demo_committed.f64').read_bytes()
    assert (output/'demo_certified.f64').read_bytes()!=(output/'demo_committed.f64').read_bytes()
print(json.dumps(report,indent=2))
print('Output:',output)
