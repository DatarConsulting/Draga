import concurrent.futures, hashlib, json, os, platform, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
import torch
ROOT=Path(__file__).resolve().parents[1]
ARMS=['main_only','single','single_param','multi_path','multi_cert']
jobs=[dict(task=t,arm=a,seed=s,steps=1800 if t=='text' else 1500)
      for t in ['text','order'] for s in range(101,106) for a in ARMS]
def run(j):
    stem=f"{j['task']}_{j['arm']}_seed{j['seed']}"
    cmd=[sys.executable,str(ROOT/'code/benchmark.py')]
    for k,v in j.items(): cmd.extend(['--'+k,str(v)])
    start=time.perf_counter()
    with (ROOT/'results'/f'{stem}.log').open('w') as f:
        p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,
                         env={**os.environ,'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    out={**j,'exit_code':p.returncode,'wall_seconds':time.perf_counter()-start}
    print(json.dumps(out),flush=True)
    return out
if __name__=='__main__':
    for j in jobs:
        if (ROOT/'results'/f"{j['task']}_{j['arm']}_seed{j['seed']}.json").exists():
            raise RuntimeError('Existing result: refusing to overwrite frozen evidence')
    manifest={'started_utc':datetime.now(timezone.utc).isoformat(),'jobs':jobs,'workers':3,
              'hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in [ROOT/'PROTOCOL.md',ROOT/'code/benchmark.py',ROOT/'data/manifest.json']},
              'python':sys.version,'torch':torch.__version__,'platform':platform.platform()}
    (ROOT/'results/run_manifest.json').write_text(json.dumps(manifest,indent=2))
    start=time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex: outcomes=list(ex.map(run,jobs))
    (ROOT/'results/run_outcomes.json').write_text(json.dumps({'runs':outcomes,'wall_seconds':time.perf_counter()-start},indent=2))
    if any(x['exit_code'] for x in outcomes): raise RuntimeError('A training job failed; inspect logs')
