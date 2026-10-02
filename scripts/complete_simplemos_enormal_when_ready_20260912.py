"""Current-run postprocessing dependency chain. Does not launch HFS or alter solvers."""
import json,os,subprocess,time
from pathlib import Path
from datetime import datetime

ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'build-release/enormal_curves_20260912'
OUT=ROOT/'reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912'
STOP=LOCAL/'stop_postprocessing'


def record(**values):
    (LOCAL/'postprocessing_status.json').write_text(json.dumps(dict(time=datetime.now().astimezone().isoformat(),coordinator_pid=os.getpid(),**values),indent=2),encoding='utf-8')


def main():
    record(state='waiting_for_two_curve_arms')
    while not all((OUT/('vela_'+arm+'_evidence.json')).exists() for arm in ('continuation','native')):
        if STOP.exists():record(state='stopped_before_postprocessing');return
        time.sleep(10)
    if STOP.exists():record(state='stopped_before_postprocessing');return
    import simplemos_enormal_curves_20260912 as e
    e.a.verify(OUT/'supplemental_source_evidence.json')
    for arm in ('continuation','native'):e.a.verify(OUT/('vela_'+arm+'_evidence.json'))
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    actions=[('comparison','qualify_simplemos_enormal_curves_20260912.py',['analyze'],'comparison_evidence.json'),
        ('fields','summarize_simplemos_enormal_curve_fields_20260912.py',[],'fields/evidence.json'),
        ('archive','finish_simplemos_enormal_curves_20260912.py',[],'completion_evidence.json')]
    for name,script,args,evidence in actions:
        if STOP.exists():record(state='stopped_between_steps',next_step=name);return
        if (OUT/evidence).exists():e.a.verify(OUT/evidence);continue
        cmd=['D:/msys64/ucrt64/bin/python.exe','-u','-X','utf8',str(ROOT/'scripts'/script),*args]
        with (LOCAL/('post_'+name+'.log')).open('x',encoding='utf-8') as log:
            child=subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
            record(state='running',step=name,child_pid=child.pid,command=cmd)
            while child.poll() is None:
                if STOP.exists():
                    child.terminate();child.wait();record(state='stopped_during_step',step=name,child_pid=child.pid);return
                time.sleep(1)
        if child.returncode:
            record(state='failed',step=name,exit_code=child.returncode);raise RuntimeError(f'{name} failed; see preserved log')
        e.a.verify(OUT/evidence)
    assert e.a.read(OUT/'final_summary.json')['all_qualified']
    record(state='completed_enormal_only',high_field_simulation_started=False)
    print('Enormal postprocessing completed; HFS remains separately gated.',flush=True)

if __name__=='__main__':main()
