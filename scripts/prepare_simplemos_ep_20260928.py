"""Prepare the agreed eight-run EP intervention without changing M60 evidence."""
from pathlib import Path
import hashlib
import json
import shutil
import difflib

ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / 'build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst/sentaurus_bundle'
OUT = ROOT / 'build/simplemos_ep_20260928'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    OUT.mkdir(exist_ok=True)
    manifest = dict(schema='vela.simplemos.ep-intervention.v1',
        intervention=dict(ExtendedPrecision=128, Digits=15, RhsMin=1e-15),
        diagnostic_expectations_only=dict(substrate_gap_relative=1e-4,
            kcl_relative=1e-4, vela_Id_relative=3e-4, direct_default_Id_relative=.01),
        acceptance_changed=False, cases=[])
    commands = ['#!/bin/bash', 'set -u', 'cd "$(dirname "$0")"',
        'sdevice -h > sdevice_banner.txt 2>&1',
        'date -Is > started.txt', ': > status.tsv']
    diffs = []
    for device in ('n23', 'n24', 'n17', 'n19'):
        target = OUT/'bundle'/device
        target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(OLD/device/'input_fps.tdr', target/'input_fps.tdr')
        for algorithm in ('default', 'direct'):
            oldname = f'm60_{algorithm}_{device}_vd_0p05'
            name = f'ep_{algorithm}_{device}_vd_0p05'
            original = (OLD/device/f'{oldname}_des.cmd').read_text()
            assert original.count('Digits=8') == 1 and 'ExtendedPrecision' not in original
            deck = original.replace('Digits=8', 'ExtendedPrecision(128) Digits=15 RhsMin=1e-15')
            deck = deck.replace(oldname, name)
            if 'Plot(FilePrefix=' not in deck:
                anchor = 'CurrentPlot(Time=(Range=(0 1) Intervals=50))'
                assert deck.count(anchor) == 1
                deck = deck.replace(anchor, anchor + '\n    Plot(FilePrefix="'+name+'_state" NoOverWrite Time=(0; 0.02; 0.04; 0.06; 0.32))')
            path = target/f'{name}_des.cmd'
            path.write_text(deck, encoding='utf-8', newline='\n')
            diffs.extend(difflib.unified_diff(original.splitlines(True), deck.splitlines(True), fromfile=oldname, tofile=name))
            manifest['cases'].append(dict(device=device, algorithm=algorithm, name=name,
                old_deck_sha256=sha(OLD/device/f'{oldname}_des.cmd'),
                deck_sha256=sha(path), tdr_sha256=sha(target/'input_fps.tdr')))
            commands += [f'(cd bundle/{device} && sdevice {name}_des.cmd > {name}.console.log 2>&1)',
                'rc=$?', f'printf "{name}\\t%s\\n" "$rc" >> status.tsv']
    commands += ['date -Is > completed.txt']
    (OUT/'run.sh').write_text('\n'.join(commands)+'\n', encoding='utf-8', newline='\n')
    manifest['runner_sha256'] = sha(OUT/'run.sh')
    manifest['preparation_script_sha256'] = sha(Path(__file__))
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (OUT/'deck_changes.diff').write_text(''.join(diffs))
    print(json.dumps(dict(output=str(OUT), cases=len(manifest['cases']), manifest_sha256=sha(OUT/'manifest.json'))))

if __name__ == '__main__':
    main()
