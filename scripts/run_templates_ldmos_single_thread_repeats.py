"""STRUMPACK 1+1 qualification, then five reuse backends: fresh D5 three-round batch."""
import argparse
import os
from pathlib import Path
import sys

from run_templates_ldmos_linked_d5 import read, write, digest, preflight
from run_templates_ldmos_analysis_reuse import compare
from run_templates_ldmos_full_curve_timing import (
    ROOT, CONFIGS, freeze, verify_frozen, completed, run_case, check_ready, require_joint_pass)
from analyze_templates_ldmos_full_curve_timing import summarize

IDS = ('U1','T1','L1','M1','S1')
CONFIGURATIONS = {key:CONFIGS[key] for key in IDS}
STAGES = ('S0','R0','R1','R2')


def schedule():
    rows=[]
    def add(stage, repeat, gate, config):
        rows.append(dict(stage=stage,profile='D5',points=31,round=repeat,gate=gate,
                         config=config,key=f'd5_p31_r{repeat}_vg{gate}_{config}'))
    # Qualified STRUMPACK curves count as its predeclared first repeat.
    for gate in (4,8):add('S0',0,gate,'T1')
    for gate in (4,8):
        order=['U1','L1','M1','S1']
        for config in (order if gate==4 else order[::-1]):add('R0',0,gate,config)
    for repeat in (1,2):
        order=list(IDS[repeat:]+IDS[:repeat])
        for gate in (4,8):
            for config in (order if gate==4 else order[::-1]):add(f'R{repeat}',repeat,gate,config)
    return rows


def validate_resume(report):
    if report['status']!='paused' or report.get('active'):
        raise ValueError('Resume requires inactive paused batch')
    if report.get('configurations')!=CONFIGURATIONS or report.get('schedule')!=schedule() or report.get('control_config')!='U1':
        raise ValueError('Single-thread repeat contract changed')


def freeze_qualification_source(out, report, old_root):
    source=read(old_root/'summary.json')
    if source['runner_sha256']!=report['runner_sha256']:
        raise ValueError('Qualification control executable differs')
    previous=read(old_root/'binary/U0.json')['runtime_sha256']
    previous_dlls={Path(p).name:sha for p,sha in previous.items() if Path(p).suffix.lower()=='.dll'}
    current_dlls={p.name:digest(p) for p in (out/'binary').glob('*.dll')}
    if previous_dlls!=current_dlls:raise ValueError('Qualification runtime DLLs differ')
    controls={}
    for gate in (4,8):
        case=next(c for c in source['cases'] if c['key']==f'd5_p31_r0_vg{gate}_U0')
        if case['status']!='pass':raise ValueError('Qualification control is not passed')
        directory=old_root/case['name'];plan=read(directory/'plan.json');ledger=read(directory/'fixed/ledger.json')
        if ledger['status']!='completed' or len(ledger['exact_points'])!=31:
            raise ValueError('Incomplete qualification control')
        for path,expected in plan['frozen_files'].items():
            if digest(Path(path))!=expected:raise ValueError('Control input changed: '+path)
            report['hashes'][path]=expected
        for path in (old_root/'summary.json',directory/'plan.json',directory/'fixed/ledger.json'):
            report['hashes'][str(path)]=digest(path)
        for point in ledger['exact_points']:
            if digest(Path(point['state']))!=point['sha256']:raise ValueError('Control state changed')
            report['hashes'][point['state']]=point['sha256']
        controls[str(gate)]=str(directory)
    require_joint_pass(source['joint_qualifications']['d5_r0_U0'])
    report['qualification_controls']=controls
    report['qualification_comparisons']={}


def check_qualification(out, report):
    done=completed(report)
    for gate in (4,8):
        key=f'd5_p31_r0_vg{gate}_T1'
        if key in done and str(gate) not in report['qualification_comparisons']:
            result=compare(Path(report['qualification_controls'][str(gate)]),out/done[key]['name'])
            report['qualification_comparisons'][str(gate)]=result
            if not result['equivalent']:raise ValueError('STRUMPACK single-thread state qualification failed')
    if all(f'd5_p31_r0_vg{g}_T1' in done for g in (4,8)):
        require_joint_pass(report['joint_qualifications']['d5_r0_T1'])
        if not all(report['qualification_comparisons'][str(g)]['equivalent'] for g in (4,8)):
            raise ValueError('STRUMPACK qualification incomplete')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--qualification-source',type=Path,default=ROOT/'reference_staging/templates_ldmos_full_curve_timing_20260920')
    parser.add_argument('--through',choices=STAGES,default='R2')
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    if not __debug__ or os.name!='nt':raise RuntimeError('Requires Windows and assertions')
    out=args.output.resolve()
    if args.resume:
        report=read(out/'summary.json');validate_resume(report);verify_frozen(report)
        if (out/'PAUSE').exists():raise ValueError('Remove PAUSE before explicit resume')
    else:
        out.mkdir(parents=True,exist_ok=False);report=freeze(out)
        report.update(configurations=CONFIGURATIONS,schedule=schedule(),control_config='U1',
            completed_stages=[],scope='D5 five backends, all solver/BLAS 1+1, reuse enabled, three fresh rounds')
        freeze_qualification_source(out,report,args.qualification_source.resolve())
    lock=out/'batch.lock'
    with lock.open('x') as handle:handle.write(str(os.getpid()))
    def save():write(out/'summary.json',report)
    try:
        verify_frozen(report)
        bundle=ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json'
        for gate in (4,8):preflight(bundle,ROOT,out/'binary/U1.json',gate)
        report.update(status='running',through=args.through);save()
        check_ready(out,report);check_qualification(out,report);save()
        for item in report['schedule']:
            if STAGES.index(item['stage'])>STAGES.index(args.through):break
            if item['key'] in completed(report):continue
            if (out/'PAUSE').exists():raise KeyboardInterrupt('PAUSE requested')
            report['current_stage']=item['stage']
            run_case(out,report,item,save)
            check_ready(out,report);check_qualification(out,report);save()
            stage=item['stage']
            if all(row['key'] in completed(report) for row in report['schedule'] if row['stage']==stage):
                if stage not in report['completed_stages']:report['completed_stages'].append(stage)
                save();write(out/'analysis.json',summarize(out))
                print('SINGLE_THREAD_STAGE_COMPLETE',stage,len(completed(report)),flush=True)
        verify_frozen(report)
        report['status']='pass' if len(completed(report))==len(report['schedule']) else 'paused'
        report['stop_reason']='complete' if report['status']=='pass' else 'stage_boundary'
        save();write(out/'analysis.json',summarize(out))
    except KeyboardInterrupt as error:
        report.update(status='paused',stop_reason='user_interrupt',error=repr(error))
    except BaseException as error:
        report.update(status='failed',error=repr(error));raise
    finally:
        save();lock.unlink()
    print('SINGLE_THREAD_BATCH_END',report['status'],len(completed(report)),flush=True)


if __name__=='__main__':main()
