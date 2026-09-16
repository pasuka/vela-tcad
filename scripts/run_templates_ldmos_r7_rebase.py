"""Paired R7 high-voltage controls for one-time local QF representation."""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil
import statistics

from run_templates_ldmos_density_projection import read, sha
from run_templates_ldmos_poisson_initialization import run_point
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_predictor_study import state_delta

CASES = {(4, 19), (4, 30), (8, 17), (8, 28)}
FLAG = 'diagnostic_near_steady_qf_rebase'


def qualify(data, reference, bias):
    gate = state_gate(data, bias)
    delta = state_delta(data, reference)
    return dict(gate=gate, max_state_delta=delta,
                qualified=gate['pass_gate'] and all(d <= t for d, t in zip(delta, (1e-8, 1e-8, 1e-8, 1e-7))))


def analyze(root):
    s = read(root/'summary.json')
    if s['status'] != 'completed_diagnostic_matrix':
        raise ValueError('Incomplete matrix')
    for name, digest in s['runtime_sha256'].items():
        if sha(root/name) != digest: raise ValueError('Runtime changed')
    expected = {(g, i, r, v) for g, i in CASES for r in range(s['repeats']) for v in ('baseline', 'rebase')}
    if len(s['runs']) != len(expected) or {(r['vg'], r['index'], r['repeat'], r['variant']) for r in s['runs']} != expected:
        raise ValueError('Missing or duplicate runs')
    loaded = {}
    for row in s['runs']:
        directory = root/row['directory']
        for name in ('input', 'output'):
            if sha(directory/(name+'.json')) != row[name+'_sha256']: raise ValueError('Evidence changed')
        if sha(Path(row['reference_path'])) != row['reference_sha256']: raise ValueError('Reference changed')
        data = read(directory/'output.json')
        qualification = qualify(data, read(Path(row['reference_path'])), row['bias_V'])
        if any(row[k] != v for k, v in qualification.items()): raise ValueError('Qualification changed')
        counts = dict(newton_updates=data['newton_updates'], attempts=len(data['history']),
                      trials=sum(h['line_search_trials'] for h in data['history']),
                      assemblies=data['performance']['assembly_calls'], factorizations=data['performance']['factorizations'])
        if any(row[k] != v for k, v in counts.items()): raise ValueError('Cost mismatch')
        loaded[row['vg'], row['index'], row['repeat'], row['variant']] = (row, read(directory/'input.json'), data)
    findings = []
    for g, i in sorted(CASES):
        for variant in ('baseline', 'rebase'):
            first = loaded[g, i, 0, variant][2]
            for repeat in range(1, s['repeats']):
                data = loaded[g, i, repeat, variant][2]
                for key in ('history', 'state_interleaved', 'referenced_state_interleaved',
                            'electron_qf_reference_V', 'hole_qf_reference_V', 'residual',
                            'carrier_row_gate', 'electrical_block_gates', 'diagnostic_stop'):
                    if data[key] != first[key]: raise ValueError('Repeated trajectory changed: '+key)
        pairs = []
        for repeat in range(s['repeats']):
            b, bc, bd = loaded[g, i, repeat, 'baseline']
            c, cc, cd = loaded[g, i, repeat, 'rebase']
            expected_cfg = dict(bc, **{FLAG: True})
            if cc != expected_cfg: raise ValueError('Non-isolated configuration')
            event = cd['near_steady_qf_rebase']['events']
            if len(event) > 1: raise ValueError('Repeated rebase')
            n = event[0]['next_iteration']-1 if event else len(cd['history'])
            prefix = copy.deepcopy(cd['history'][:n])
            for h in prefix: h.pop('near_steady_qf_active', None)
            if prefix != bd['history'][:n]: raise ValueError('Pre-trigger history changed')
            if event:
                old = bd['history'][n]
                if event[0]['merit_before'] >= 1e-9 or not all(x['satisfied'] for x in old['iteration_trace']['before_gates']['blocks']):
                    raise ValueError('Trigger before original block gates')
                if old['iteration_trace']['before_gates']['row']['satisfied']: raise ValueError('Unneeded trigger')
            elif any(cd[k] != bd[k] for k in ('state_interleaved', 'referenced_state_interleaved', 'residual', 'carrier_row_gate')):
                raise ValueError('Inactive rebase changed result')
            if 'pseudo_transient' in cd or 'density_projection' in cd: raise ValueError('Unexpected algorithm')
            pairs.append(dict(repeat=repeat, qualified=b['qualified'] and c['qualified'],
                baseline={k:b[k] for k in ('newton_updates','attempts','trials','assemblies','factorizations','wall_seconds')},
                rebase={k:c[k] for k in ('newton_updates','attempts','trials','assemblies','factorizations','wall_seconds')},
                event=event, state_delta_to_paired_baseline=state_delta(cd,bd), prefix_exact=True,
                baseline_floor=sum(h['scaled_l2_before']<1e-9 for h in bd['history']),
                rebase_floor=sum(h['scaled_l2_before']<1e-9 for h in cd['history']),
                final_row=cd['carrier_row_gate']['max_ratio']))
        findings.append(dict(vg=g,index=i,pairs=pairs,repeated_trajectory_exact=True,
            median_baseline_seconds=statistics.median(p['baseline']['wall_seconds'] for p in pairs),
            median_rebase_seconds=statistics.median(p['rebase']['wall_seconds'] for p in pairs)))
    result=dict(summary_sha256=sha(root/'summary.json'), findings=findings,
        scope='Frozen original seeds, identical inputs except opt-in representation. All attempts and rebuilds charged. Child process wall includes trace/output; not full-curve production timing.')
    (root/'findings.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps(result,allow_nan=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path);p.add_argument('--probe',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--analyze',action='store_true');args=p.parse_args()
    if args.analyze: analyze(args.output);return
    if args.repeats<1 or args.manifest is None or args.probe is None: p.error('Manifest, probe and positive repeats required')
    cases=[c for c in read(args.manifest) if c['variant']=='linear_baseline' and (c['gate_V'],c['baseline_index']) in CASES]
    if len(cases)!=4: raise ValueError('Expected four frozen high-voltage seeds')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=out/args.probe.name;shutil.copy2(args.probe,binary)
    for path in args.probe.parent.glob('*.dll'): shutil.copy2(path,out/path.name)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
    s=dict(schema='vela.r7_local_rebase.v1',status='running',repeats=args.repeats,
           manifest_sha256=sha(args.manifest),runtime_sha256={f.name:sha(f) for f in [binary,*out.glob('*.dll')]},runs=[])
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(s,indent=2,allow_nan=False),encoding='utf-8');os.replace(tmp,out/'summary.json')
    save()
    for repeat in range(args.repeats):
        for case in cases:
            original=read(Path(case['input']));reference=read(Path(case['expected']))
            for variant in (('baseline','rebase') if repeat%2==0 else ('rebase','baseline')):
                cfg=copy.deepcopy(original)
                if variant=='rebase': cfg[FLAG]=True
                name='vg%d_p%d_r%d_%s'%(case['gate_V'],case['baseline_index'],repeat,variant)
                row,data=run_point(binary,out/name,cfg)
                row.update(vg=case['gate_V'],index=case['baseline_index'],bias_V=case['bias_V'],repeat=repeat,variant=variant,
                    source_input_sha256=sha(Path(case['input'])),reference_path=case['expected'],reference_sha256=sha(Path(case['expected'])))
                if data is not None: row.update(qualify(data,reference,case['bias_V']))
                s['runs'].append(row);save()
                if data is None or row['exit_code']: raise RuntimeError('Probe execution failed; saved all evidence')
    s['status']='completed_diagnostic_matrix';save();analyze(out)


if __name__=='__main__': main()
