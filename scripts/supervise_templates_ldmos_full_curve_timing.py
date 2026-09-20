"""Advance an authorized frozen batch only across qualified stage boundaries."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

from run_templates_ldmos_linked_d5 import read, write, digest
from run_templates_ldmos_full_curve_timing import CONFIGS, ROOT, STAGES, schedule
from analyze_templates_ldmos_full_curve_timing import summarize


def review_reasons(report, through):
    reasons=[]
    if report['status'] not in ('paused','pass') or report.get('stop_reason') not in ('stage_boundary','complete'):
        return ['Batch is not at an ordinary completed stage boundary']
    done={c['key']:c for c in report['cases'] if c['status']=='pass'}
    expected=[x for x in schedule() if STAGES.index(x['stage'])<=STAGES.index(through)]
    if any(x['key'] not in done for x in expected): reasons.append('Stage is incomplete')
    if reasons: return reasons
    for key,result in report['joint_qualifications'].items():
        if any(result.get(k,{}).get('status')!='pass' for k in ('engineering','final')):
            reasons.append('Joint qualification failed: '+key)
    if through=='P1':
        for config in ('L1','M1','S1'):
            for gate in (4,8):
                base=done[f'd5_p31_r0_vg{gate}_U1'];candidate=done[f'd5_p31_r0_vg{gate}_{config}']
                if candidate['wall_seconds']<base['wall_seconds']:
                    reasons.append(f'{config} faster than U1 at Vg{gate}; inspect repeat shortlist')
                if candidate['audit']['rollbacks']<base['audit']['rollbacks']:
                    reasons.append(f'{config} has fewer outer rollbacks at Vg{gate}')
    return reasons


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();root=args.root.resolve()
    lock=root/'supervisor.lock'
    with lock.open('x') as handle: handle.write(str(os.getpid()))
    state=dict(status='running',pid=os.getpid(),script_sha256=digest(Path(__file__)),
               scope='Advance existing authorized batch, stop on failure, user pause or shortlist review',events=[])
    def save():write(root/'supervision.json',state)
    save()
    try:
        while True:
            report=read(root/'summary.json')
            if report['status']=='running':
                # A missing owner lock never licenses an automatic duplicate run.
                if not (root/'batch.lock').exists():
                    state.update(status='review_required',reasons=['Running batch has no owner lock']);break
                time.sleep(20);continue
            if (root/'PAUSE').exists() or report.get('stop_reason')=='user_interrupt':
                state.update(status='paused',reasons=['User pause preserved']);break
            stage=report.get('through','P1');reasons=review_reasons(report,stage)
            if reasons:
                state.update(status='review_required',reasons=reasons);break
            # Wait for the current driver's final lock removal before resuming.
            if (root/'batch.lock').exists():time.sleep(1);continue
            result=summarize(root);write(root/'analysis.json',result)
            state['events'].append(dict(stage=stage,completed=result['completed_cases'],status=report['status']))
            if report['status']=='pass':state['status']='complete';break
            next_stage=STAGES[STAGES.index(stage)+1]
            command=[sys.executable,ROOT/'scripts/run_templates_ldmos_full_curve_timing.py',
                     '--output',root,'--resume','--through',next_stage]
            state['next_stage']=next_stage;save()
            print('ADVANCING',next_stage,flush=True)
            env=dict(os.environ);env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env.get('PATH','')
            with (root/f'supervised_{next_stage}.log').open('x') as log:
                proc=subprocess.Popen(list(map(str,command)),cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
                state['driver_pid']=proc.pid;save()
                while proc.poll() is None:
                    try:proc.wait(timeout=20)
                    except subprocess.TimeoutExpired:pass
                state['events'].append(dict(stage=next_stage,returncode=proc.returncode));save()
                if proc.returncode:
                    state.update(status='review_required',reasons=[f'{next_stage} driver failed']);break
    except BaseException as error:
        state.update(status='review_required',error=repr(error));raise
    finally:
        save();lock.unlink()
    print('SUPERVISION_STOPPED',state['status'],flush=True)


if __name__=='__main__':main()
