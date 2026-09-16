"""Bounded frozen-state controls for mass, density inversion and QF representation."""
import argparse
import copy
from decimal import Decimal, localcontext
import json
import os
from pathlib import Path
import shutil

from analyze_templates_ldmos_predictor_study import state_delta
from run_templates_ldmos_density_projection import read, sha
from run_templates_ldmos_electrothermal_curve import state_gate
from run_templates_ldmos_poisson_initialization import restored_config, run_point


def rebase_state(state, reference):
    """Change QF representation only; never substitute the reference solution."""
    state_delta(state, reference)
    result=copy.deepcopy(state)
    with localcontext() as context:
        context.prec=80
        for k,key in ((1,'electron_qf_reference_V'),(2,'hole_qf_reference_V')):
            for i,new_reference in enumerate(reference[key]):
                physical=Decimal(state[key][i])+Decimal(state['referenced_state_interleaved'][4*i+k])
                result[key][i]=new_reference
                result['referenced_state_interleaved'][4*i+k]=float(physical-Decimal(new_reference))
                result['state_interleaved'][4*i+k]=float(physical)
    if any(d>1e-14 for d in state_delta(state,result)):
        raise ValueError('Rebase cannot preserve physical state at diagnostic precision')
    return result


def zero_small_qf_references(state):
    """Experimental state-only rule, matching the existing 1 mV recenter window."""
    reference=copy.deepcopy(state)
    with localcontext() as context:
        context.prec=80
        for k,key in ((1,'electron_qf_reference_V'),(2,'hole_qf_reference_V')):
            for i,old in enumerate(state[key]):
                physical=Decimal(old)+Decimal(state['referenced_state_interleaved'][4*i+k])
                if abs(physical)<=Decimal('0.001'):reference[key][i]=0.
    return rebase_state(state,reference)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-attempts',type=int,default=8)
    p.add_argument('--reference-followup',action='store_true',help='State-only small-QF rebase controls, no reference-solution representation')
    args=p.parse_args()
    if args.max_attempts<1:raise ValueError('Positive attempt budget required')
    source=args.source.resolve();summary=read(source/'summary.json')
    if summary['status']!='completed_diagnostic_matrix':raise ValueError('Source incomplete')
    if sha(args.manifest)!=summary['manifest_sha256']:raise ValueError('Manifest changed')
    if sha(source/'electrothermal_probe.exe')!=summary['binary_sha256']:raise ValueError('Binary changed')
    for name,digest in summary['runtime_sha256'].items():
        if sha(source/name)!=digest:raise ValueError('Runtime changed')
    cases=[c for c in read(args.manifest) if c['variant']=='linear_baseline'
           and (c['gate_V'],c['baseline_index']) in ((4,11),(8,9))]
    if len(cases)!=2:raise ValueError('Expected both frozen points')
    out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=out/'electrothermal_probe.exe';shutil.copy2(source/binary.name,binary)
    for name in summary['runtime_sha256']:shutil.copy2(source/name,out/name)
    if os.name=='nt':
        import ctypes
        ctypes.windll.kernel32.SetErrorMode(0x0001|0x0002|0x8000)
    report=dict(schema='vela.near_steady_isolation.v1',status='running',
        source_directory=str(source),source_summary_sha256=sha(source/'summary.json'),
        binary_sha256=sha(binary),runtime_sha256=summary['runtime_sha256'],
        manifest_sha256=sha(args.manifest),max_attempts=args.max_attempts,
        reference_followup=args.reference_followup,sources=[],runs=[])
    def save():
        tmp=out/'summary.tmp';tmp.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
        os.replace(tmp,out/'summary.json')
    save()
    for case in cases:
        vg=case['gate_V'];label='vg%d_p%d'%(vg,case['baseline_index'])
        record=next(r for r in summary['runs'] if r['vg']==vg and r['variant']=='ptc_defect_model_s1')
        directory=source/record['directory']
        if sha(directory/'output.json')!=record['output_sha256']:raise ValueError('Terminal changed')
        terminal=read(directory/'output.json');reference=read(Path(case['expected']))
        original=read(Path(case['input']));base=restored_config(original,terminal)
        rebased=restored_config(original,rebase_state(terminal,reference))
        report['sources'].append(dict(vg=vg,terminal_sha256=sha(directory/'output.json'),
            reference_path=case['expected'],reference_sha256=sha(Path(case['expected'])),
            original_input_path=case['input'],original_input_sha256=sha(Path(case['input'])),
            terminal_tau_s=terminal['pseudo_transient']['next_tau_s'],
            rebase_delta=state_delta(base,rebased),prior_charged_cost=record['charged_cost']))
        def run(variant,cfg):
            row,data=run_point(binary,out/(variant+'_'+label),cfg)
            row.update(variant=variant,vg=vg,qualified=False)
            if data is not None:
                gate=state_gate(data,case['bias_V']);delta=state_delta(data,reference)
                row.update(gate=gate,max_state_delta=delta,state_change=state_delta(terminal,data),
                    row_gate=data['carrier_row_gate'],blocks=data['electrical_block_gates'],
                    qualified=row['exit_code']==0 and gate['pass_gate']
                    and all(d<=t for d,t in zip(delta,(1e-8,1e-8,1e-8,1e-7))))
            report['runs'].append(row);save()
            if row['exit_code'] or data is None:raise RuntimeError('Probe failed: '+variant)
            return data
        if args.reference_followup:
            local=zero_small_qf_references(base)
            run('audit_local',dict(local,diagnostic_newton_max_iterations=0))
            run('r7_local',dict(local,diagnostic_newton_max_iterations=args.max_attempts))
            run('density_local',dict(local,diagnostic_newton_max_iterations=args.max_attempts,diagnostic_density_projection='v1'))
            continue
        for name,cfg in (('audit_terminal',base),('audit_rebased',rebased),
                         ('audit_reference',restored_config(original,reference))):
            run(name,dict(cfg,diagnostic_newton_max_iterations=0))
        for name,cfg in (('r7',base),('r7_rebased',rebased),('density',dict(base,diagnostic_density_projection='v1'))):
            run(name,dict(cfg,diagnostic_newton_max_iterations=args.max_attempts))
        # Calibrate the near-state median(Mii/Jii), then use the original terminal tau.
        # Both mass controls start with this same operator; subsequent SER is unchanged.
        probe=run('mass_calibration',dict(base,diagnostic_pseudo_transient=True,
                  diagnostic_newton_max_iterations=1,diagnostic_pseudo_direction_audit=True))
        actual=probe['pseudo_transient'].get('initial_tau_s')
        if not actual or actual<=0:raise ValueError('No mass calibration at unqualified terminal')
        scale=terminal['pseudo_transient']['next_tau_s']/actual
        for name,projection in (('mass_qf','off'),('mass_density','v1')):
            data=run(name,dict(base,diagnostic_pseudo_transient=True,diagnostic_pseudo_time_scale=scale,
                diagnostic_pseudo_acceptance='defect_model',diagnostic_density_projection=projection,
                diagnostic_newton_max_iterations=args.max_attempts))
            got=data['pseudo_transient']['initial_tau_s'];wanted=terminal['pseudo_transient']['next_tau_s']
            if abs(got/wanted-1)>1e-12:raise ValueError('Initial tau isolation failed')
    report['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':main()
