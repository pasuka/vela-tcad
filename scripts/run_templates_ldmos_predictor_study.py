"""Serial diagnostic predictor study, with frozen R7 inputs and no acceptance relaxation.

Phase A replays the accepted R7 target sequence from its qualified zero-bias
state; it excludes initialization and is NOT an end-to-end qualification.
Keep all attempted solves, including failures, and report screening overhead.
"""
import argparse
import hashlib
import json
import os
from electrothermal_state import read_bound_record, write as write_state_record
from state_archive import mesh_identity
from pathlib import Path, PurePosixPath
import subprocess
import time

from run_templates_ldmos_electrothermal_curve import state_gate


def read(p):
    return read_bound_record(p)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root',type=Path,required=True)
    parser.add_argument('--profile',type=Path,required=True)
    parser.add_argument('--runner',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--maximum-bias',type=float,default=28/3)
    parser.add_argument('--phase',choices=('A','B','C','D'),default='A')
    parser.add_argument('--variants',nargs='+',choices=('baseline','residual','nearest','actual_step','tangent','localfit'),default=['baseline','residual','nearest'])
    args=parser.parse_args()
    root=args.evidence_root.resolve();out=args.output.resolve();runner=args.runner.resolve()
    profile=read(args.profile)
    for name,digest in profile['files_sha256'].items():
        if sha(root/name)!=digest:raise ValueError('Frozen input changed: '+name)
    recorded=PurePosixPath(profile['dependencies'][0]['recorded_path']).parents[2]
    def mapped(name):
        relative=PurePosixPath(name).relative_to(recorded)
        if '..' in relative.parts:raise ValueError('Invalid evidence path')
        return root.joinpath(*relative.parts)
    out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    allowed={'A':('baseline','residual','nearest'),'B':('baseline','actual_step'),'C':('baseline','tangent'),'D':('baseline','nearest','localfit')}[args.phase]
    if any(v not in allowed for v in args.variants):raise ValueError('Variant is not part of this phase')
    report=dict(schema='vela.ldmos.predictor_study.v1',status='running',phase=args.phase,
        scope=({'A':'Accepted R7 fixed target sequence','B':'Adaptive reference versus actual successful step growth','C':'Accepted R7 fixed target sequence','D':'R8 density updates with adaptive history controls'}[args.phase])+' from qualified zero state; no initialization; no new full-curve or native-time qualification',
        runner_sha256=sha(runner),profile_sha256=sha(args.profile),maximum_bias_V=args.maximum_bias,runs=[])
    def save():
        p=out/'summary.tmp';p.write_text(json.dumps(report,indent=2),encoding='utf-8');os.replace(str(p),str(out/'summary.json'))
    save()
    for gate in (4,8):
        base=read(root/'results'/('r7_full_vg%d'%gate)/'ledger.json')
        accepted=[r for r in base['runs'] if r['gate']['pass_gate'] and r['bias_V']<=args.maximum_bias+1e-12]
        targets=[r['bias_V'] for r in accepted]
        if args.phase in ('B','D'):
            targets=[v for v in read(root/'cases_r7'/('vela_vg%d.json'%gate))['sweep']['bias_points_V'] if v<=args.maximum_bias+1e-12]
        if not targets or targets[0]!=0.:raise ValueError('Missing zero target')
        initial_path=mapped(accepted[0]['directory'])/'output.json';initial=read(initial_path)
        if not state_gate(initial,0.)['pass_gate']:raise ValueError('Unqualified initial state')
        for variant in args.variants:
            case=out/('%s_vg%d'%(variant,gate));case.mkdir()
            cfg=read(root/'cases_r7/data'/('input_vg%d.json'%gate))
            cfg['mesh_file']=str(root/'cases/data/mesh.json')
            if 'ialmob' in cfg['mobility_SI']:
                cfg['mobility_SI']['ialmob']['geometry_file']=str(root/'cases/data/ialmob_geometry.json')
            for key in ('state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V'):
                cfg[key]=initial[key]
            cfg['initialization']='provided_state'
            cfg['diagnostic_density_update_iterations']=60 if args.phase=='D' else 0
            inp=case/'input.json'
            mesh=read(Path(cfg['mesh_file']))
            write_state_record(inp,cfg,dict(mode='electrothermal',mesh_sha256=mesh_identity(mesh,cfg.get('coordinate_to_metres',1.)),potential_origin_V=cfg.get('potential_origin_V',0.)))
            deck=read(root/'cases_r7'/('vela_vg%d.json'%gate))
            deck.update(input_file=str(inp),output_directory=str(case/'results'),initialization=dict(mode='provided_state'))
            deck['sweep'].update(step_policy='fixed_targets' if args.phase in ('A','C') else 'actual_step' if variant=='actual_step' else 'adaptive',bias_points_V=targets,
                predictor='tangent_guarded' if variant=='tangent' else 'local_guarded' if variant=='localfit' else 'linear',
                predictor_guard='residual' if variant in ('residual','nearest') else 'none',
                predictor_history='nearest' if variant=='nearest' else 'legacy')
            if args.phase=='D':deck['sweep'].update(density_update_requires_prediction=True,density_update_maximum_bias_V=1.)
            deckpath=case/'deck.json';deckpath.write_text(json.dumps(deck,indent=2),encoding='utf-8')
            row=dict(gate_V=gate,variant=variant,initial_state_sha256=sha(initial_path),input_sha256=sha(inp),deck_sha256=sha(deckpath),targets_V=targets)
            before=None
            if os.name!='nt':
                import resource
                before=resource.getrusage(resource.RUSAGE_CHILDREN)
            if sha(runner)!=report['runner_sha256']:raise ValueError('Runner changed during experiment')
            with (case/'run.log').open('x') as log:
                started=time.perf_counter()
                process=subprocess.Popen([str(runner),'--config',str(deckpath)],stdout=log,stderr=subprocess.STDOUT,env=env)
                while True:
                    try:code=process.wait(timeout=30);break
                    except subprocess.TimeoutExpired:
                        progress=dict(gate_V=gate,variant=variant,elapsed_s=time.perf_counter()-started)
                        ledgerpath=case/'results/ledger.json'
                        if ledgerpath.exists():
                            ledger=read(ledgerpath);progress['accepted_bias_V']=ledger['accepted_bias_V']
                        print(json.dumps(progress),flush=True)
                row.update(exit_code=code,wall_seconds=time.perf_counter()-started)
            if before is not None:
                after=resource.getrusage(resource.RUSAGE_CHILDREN)
                row['child_cpu_seconds']=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime
            ledger=read(case/'results/ledger.json');row['status']=ledger['status'];row['attempts']=len(ledger['runs'])
            row['failed_attempts']=sum(not r['gate']['pass_gate'] for r in ledger['runs'])
            details=[]
            for run in ledger['runs']:
                path=Path(run['directory'])/'output.json'
                if not path.exists():
                    details.append(dict(bias_V=run['bias_V'],solver_error=run['gate']));continue
                result=read(path);check=state_gate(result,run['bias_V'])
                if check['pass_gate']!=run['gate']['pass_gate']:raise ValueError('Independent gate disagrees')
                perf=result.get('performance',{})
                details.append(dict(bias_V=run['bias_V'],gate=check,newton_updates=result['newton_updates'],
                    line_search_trials=sum(h.get('line_search_trials',0) for h in result['history']),
                    prediction=run['prediction'],performance=perf))
            row['details']=details
            for name in ('newton_updates','line_search_trials'):row[name]=sum(d.get(name,0) for d in details)
            for name in ('assembly_calls','factorizations','residual_only_calls','assembly_seconds','factorization_seconds'):
                row[name]=sum(d.get('performance',{}).get(name,0) for d in details)
            row['screening_seconds']=sum(d.get('prediction',{}).get('screening',{}).get('wall_seconds',0) for d in details)
            row['tangent_preparation_seconds']=sum(d.get('prediction',{}).get('tangent_preparation',{}).get('wall_seconds',0) for d in details)
            row['local_fit_preparation_seconds']=sum(d.get('prediction',{}).get('local_fit',{}).get('wall_seconds',0) for d in details)
            report['runs'].append(row);save()
            print(json.dumps({k:v for k,v in row.items() if k not in ('details','targets_V')}),flush=True)
    report['status']='completed_diagnostic_matrix';save()


if __name__=='__main__':
    main()
