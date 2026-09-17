"""Frozen R9 seeds: paired thermal-HFS reuse and separate trace/timer diagnostics."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from analyze_templates_ldmos_contact_sweep import numerical
from analyze_templates_ldmos_newton_phases import analyze
from run_templates_ldmos_contact_sweep import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate


def comparable(result):
    result = numerical(result)
    for row in result['history']:
        row.pop('iteration_trace', None)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('reference', 'probe', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    runtime = out/'runtime'
    runtime.mkdir()
    probe = args.probe.resolve()
    exe = runtime/probe.name
    shutil.copy2(probe, exe)
    for path in probe.parent.glob('*.dll'):
        shutil.copy2(path, runtime/path.name)
    frozen = {p.name: sha(p) for p in runtime.iterdir()}
    ledger = read(args.reference/'results/ledger.json')
    report = dict(status='running', runtime_sha256=frozen, runs=[], diagnostics=[])
    def save():
        (out/'summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    save()
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    try:
        for target in (.375, 4., 16./3.):
            matches = [r for r in ledger['runs'] if abs(r['bias_V']-target)<1e-12]
            assert len(matches)==1
            source = Path(matches[0]['directory'])
            cfg, expected = read(source/'input.json'), comparable(read(source/'output.json'))
            for repeat, flags in enumerate(((False, True), (True, False), (True,))):
                diagnostic = repeat==2
                for enabled in flags:
                    name = 'v%g_r%d_%s'%(target, repeat, 'reuse' if enabled else 'base')
                    directory = out/name
                    directory.mkdir()
                    current = copy.deepcopy(cfg)
                    current.update(reuse_ialmob_thermal_high_field=enabled,
                        diagnostic_ialmob_kernel_timing=diagnostic, diagnostic_iteration_trace=diagnostic)
                    input_path = directory/'input.json'
                    input_path.write_text(json.dumps(current), encoding='utf-8')
                    start = time.perf_counter()
                    with (directory/'run.log').open('x', encoding='utf-8') as log:
                        subprocess.run([str(exe), str(input_path), str(directory/'output.json')],
                            cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
                    wall = time.perf_counter()-start
                    result = read(directory/'output.json')
                    assert comparable(result)==expected, (name, 'numerical trajectory changed')
                    gate = state_gate(result, target)
                    assert gate['pass_gate'], (name, gate)
                    row = dict(name=name, bias_V=target, repeat=repeat, reuse=enabled,
                        external_wall_seconds=wall, gate=gate, numerical_exact_to_frozen=True,
                        newton_updates=result['newton_updates'],
                        trials=sum(h['line_search_trials'] for h in result['history']),
                        performance=result['performance'], source_input_sha256=sha(source/'input.json'),
                        input_sha256=sha(input_path), output_sha256=sha(directory/'output.json'))
                    if diagnostic:
                        row['phase_analysis']=analyze(result)
                        candidates=[c for h in result['history'] for c in h['iteration_trace'].get('line_search_candidates',[])]
                        assert len(candidates)==row['trials']
                        row['candidate_cost']=dict(count=len(candidates),
                            rejected=sum(not c['accepted'] for c in candidates),
                            accepted_seconds=sum(c['trial_seconds'] for c in candidates if c['accepted']),
                            rejected_seconds=sum(c['trial_seconds'] for c in candidates if not c['accepted']))
                        report['diagnostics'].append(row)
                    else:
                        report['runs'].append(row)
                    save()
                    print(json.dumps({k:v for k,v in row.items() if k not in ('phase_analysis','performance')}),flush=True)
        for name,digest in frozen.items():
            assert sha(runtime/name)==digest
        report['status']='complete'
        save()
    except BaseException as error:
        report.update(status='failed',error=repr(error))
        save()
        raise


if __name__=='__main__':
    main()
