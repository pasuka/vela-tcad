"""Four-run ablation: remove EP(128) from frozen decks, keep all other settings."""
from pathlib import Path
import difflib
import hashlib
import json
import shutil

ROOT=Path(__file__).resolve().parents[1]
EP=ROOT/'build/simplemos_ep_20260928'
OUT=ROOT/'build/simplemos_noep_20260929'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=False)
    upstream=json.loads((EP/'manifest.json').read_text())
    cases=[];diff=[]
    for c in upstream['cases']:
        if c['device'] not in ('n23','n24'):continue
        src=EP/'bundle'/c['device']/(c['name']+'_des.cmd')
        tdr=src.parent/'input_fps.tdr'
        assert sha(src)==c['deck_sha256'] and sha(tdr)==c['tdr_sha256']
        name=c['name'].replace('ep_','noep_',1)
        before=src.read_text();assert before.count('ExtendedPrecision(128) ')==1
        after=before.replace('ExtendedPrecision(128) ','').replace(c['name'],name)
        assert after.replace(name,c['name'])==before.replace('ExtendedPrecision(128) ','')
        folder=OUT/'bundle'/c['device'];folder.mkdir(parents=True,exist_ok=True)
        deck=folder/(name+'_des.cmd');deck.write_text(after,encoding='utf-8',newline='\n')
        shutil.copy2(tdr,folder/'input_fps.tdr')
        cases.append(dict(device=c['device'],algorithm=c['algorithm'],name=name,
            source_ep_name=c['name'],source_ep_deck_sha256=sha(src),deck_sha256=sha(deck),
            tdr_sha256=sha(tdr),only_ep_removed_except_output_names=True))
        diff.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile=c['name'],tofile=name))
    assert len(cases)==4
    runner='''#!/bin/bash
set -u
cd "$(dirname "$0")"
mkdir execution_lock || exit 2
sha256sum -c inputs.sha256 > input_verification.txt || exit 3
sdevice -h > sdevice_banner.txt 2>&1
date -Is > started.txt
run_device() {
  local device="$1"
  for algorithm in default direct; do
    local name="noep_${algorithm}_${device}_vd_0p05"
    date -Is > "${name}.started.txt"
    (cd "bundle/$device" && sdevice "${name}_des.cmd" > "${name}.console.log" 2>&1)
    local rc=$?
    printf '%s\\n' "$rc" > "${name}.exitcode"
    date -Is > "${name}.completed.txt"
  done
}
run_device n23 &
first=$!
run_device n24 &
second=$!
wait "$first"
wait "$second"
date -Is > completed.txt
sha256sum bundle/*/* > results.sha256
'''
    (OUT/'run.sh').write_text(runner,encoding='utf-8',newline='\n')
    manifest=dict(schema='vela.simplemos.noep-ablation.v1',cases=cases,
        intervention='Remove ExtendedPrecision(128) only; rename outputs for isolation',
        unchanged_math=dict(Digits=15,ErrRef_electron_cm3=100,ErrRef_hole_cm3=100,RhsMin=1e-15,Iterations=20),
        execution='two independent devices in parallel, default/direct sequential per device',
        acceptance_changed=False,prepare_script_sha256=sha(Path(__file__)),runner_sha256=sha(OUT/'run.sh'))
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'deck_changes.diff').write_text(''.join(diff))
    files=[OUT/'run.sh',OUT/'manifest.json',*sorted((OUT/'bundle').rglob('*'))]
    (OUT/'inputs.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n' for p in files if p.is_file()),newline='\n')
    print(json.dumps(dict(cases=len(cases),output=str(OUT),manifest_sha256=sha(OUT/'manifest.json'))))

if __name__=='__main__':main()
