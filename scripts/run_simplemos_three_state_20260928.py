"""Run qualified fixed-state probes on frozen Codespaces binary; no solve."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

BASE=Path('/workspaces/simplemos-hfs-merged-20260926')
OUT=BASE/'ep_review_20260928'
EXE=BASE/'source-continuation-v2/build/vela_example_runner'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(devices):
    manifest=json.loads((OUT/'manifest.json').read_text())
    assert sha(EXE)==manifest['binary_expected_sha256']
    for name,digest in manifest['files'].items():assert sha(OUT/name)==digest,name
    expected=json.loads((OUT/'expected_currents.json').read_text())
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VELA_LINEAR_THREADS='1')
    results=[]
    for device in devices:
        for arm in ('a','b','c'):
            for contact in ('drain','source','substrate'):
                name=f'{device}_{arm}_{contact}'
                cfg=OUT/(name+'.json')
                output=OUT/(name+'.result.json')
                if output.exists():
                    r=json.loads(output.read_text())
                    assert r['config_sha256']==sha(cfg)
                else:
                    # Preserve the manually verified initial (a) replay.
                    stdout=OUT/(name+'.probe.stdout.txt');stderr=OUT/(name+'.probe.stderr.txt')
                    assert not stdout.exists() and not stderr.exists()
                    start=time.monotonic()
                    with stdout.open('w') as o,stderr.open('w') as e:
                        p=subprocess.run([str(EXE),'--config',str(cfg),'--log','off'],stdout=o,stderr=e,env=env)
                    r=json.loads(stdout.read_text().strip().splitlines()[-1])
                    r.update(exit_code=p.returncode,elapsed_seconds=time.monotonic()-start,config_sha256=sha(cfg),state_sha256=sha(OUT/f'{device}_{arm}.h5'))
                    output.write_text(json.dumps(r,indent=2)+'\n')
                assert r['exit_code']==0 and r['converged'] and r['read_only'],name
                if arm=='a':
                    target=expected[device][contact]
                    for field in ('current_A_per_um','contact_current_extractor_A_per_um'):
                        assert abs(r[field]-target)<=1e-12*abs(target),(name,field,r[field],target)
                results.append(dict(device=device,arm=arm,contact=contact,
                    operator=r['current_A_per_um'],total=r['contact_current_extractor_A_per_um'],
                    electron=r['contact_current_extractor_electron_A_per_um'],
                    hole_internal=r['contact_current_extractor_hole_A_per_um'],
                    hole_conventional=-r['contact_current_extractor_hole_A_per_um']))
                print(name,'passed',flush=True)
        keyed={(r['arm'],r['contact']):r for r in results if r['device']==device}
        differences=[]
        for contact in ('drain','source','substrate'):
            for lhs,rhs,meaning in (('b','a','projection_repacking'),('c','b','export_precision_state_response')):
                a,b=keyed[lhs,contact],keyed[rhs,contact]
                differences.append(dict(contact=contact,meaning=meaning,
                    **{f'{k}_delta_A_per_um':a[k]-b[k] for k in ('operator','total','electron','hole_conventional')},
                    total_delta_over_original_drain=(a['total']-b['total'])/abs(keyed['a','drain']['total'])))
        (OUT/f'{device}_three_state_summary.json').write_text(json.dumps(dict(binary_sha256=sha(EXE),
            values=[r for r in results if r['device']==device],differences=differences,
            all_split_restores_passed=True,baseline_reproduction_relative_limit=1e-12),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('devices',nargs='+',choices=['n23','n24']);a=p.parse_args();main(a.devices)
