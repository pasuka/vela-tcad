"""M60 default-observer field export at all 816 existing target biases.

Only output names and Plot sampling change. Existing M60 references are read-only.
"""
import hashlib
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_engineering_20260929/native_fields'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def without_plot(text):
    return '\n'.join(line.strip() for line in text.splitlines() if line.strip() and 'Plot(FilePrefix=' not in line)


def main():
    OUT.mkdir(exist_ok=False)
    old=ROOT/'build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst'
    manifest=json.loads((old/'sentaurus_manifest.json').read_text())
    cases=[]
    for c in manifest['cases']:
        if c['algorithm']!='default':continue
        src=ROOT/c['deck'];tdr=ROOT/c['input_tdr']
        assert sha(src)==c['deck_sha256'] and sha(tdr)==c['input_tdr_sha256']
        name='m60fields_'+c['base_case'];folder=OUT/'bundle'/c['device'];folder.mkdir(parents=True,exist_ok=True)
        original=src.read_text();text='\n'.join(line for line in original.splitlines() if 'Plot(FilePrefix=' not in line)
        text=text.replace(c['run_case'],name)
        anchor='CurrentPlot(Time=(Range=(0 1) Intervals=50))';assert text.count(anchor)==1
        text=text.replace(anchor,anchor+'\n    Plot(FilePrefix="'+name+'_state" NoOverWrite Time=(Range=(0 1) Intervals=50))')+'\n'
        assert without_plot(text.replace(name,c['run_case']))==without_plot(original)
        path=folder/(name+'_des.cmd');path.write_text(text,newline='\n')
        shutil.copy2(tdr,folder/'input_fps.tdr')
        cases.append(dict(name=name,case=c['base_case'],device=c['device'],vd=c['drain_voltage_V'],
            source_deck_sha256=sha(src),deck_sha256=sha(path),tdr_sha256=sha(tdr),expected_fields=51))
    assert len(cases)==16
    lines=['#!/bin/bash','set -u','cd "$(dirname "$0")"','mkdir execution_lock || exit 2',
        'sha256sum -c inputs.sha256 > input_verification.txt || exit 3','date -Is > started.txt',
        'run_case() {','  local dev="$1"','  local name="$2"','  date -Is > "$name.started.txt"',
        '  (cd "bundle/$dev" && sdevice "${name}_des.cmd" > "$name.console.log" 2>&1)',
        '  local rc=$?','  printf "%s\\n" "$rc" > "$name.exitcode"','  date -Is > "$name.completed.txt"','}',
        'lane_a() {']
    for c in cases[::2]:lines.append(f'  run_case {c["device"]} {c["name"]}')
    lines+=['}','lane_b() {']
    for c in cases[1::2]:lines.append(f'  run_case {c["device"]} {c["name"]}')
    lines+=['}','lane_a &','first=$!','lane_b &','second=$!','wait "$first"','wait "$second"',
        'date -Is > completed.txt','sha256sum bundle/*/* > results.sha256']
    (OUT/'run.sh').write_text('\n'.join(lines)+'\n',newline='\n')
    (OUT/'manifest.json').write_text(json.dumps(dict(cases=cases,numerical_settings_changed=False,
        original_reference_replaced=False,purpose='Independent initial states and full physical-field coverage'),indent=2)+'\n')
    files=[OUT/'run.sh',OUT/'manifest.json',*sorted((OUT/'bundle').rglob('*'))]
    (OUT/'inputs.sha256').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT).as_posix()}\n' for p in files if p.is_file()),newline='\n')
    print(json.dumps(dict(runs=16,target_states=816,path=str(OUT))))


if __name__=='__main__':main()
