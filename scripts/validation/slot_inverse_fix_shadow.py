"""Bounded read-only P2 stage analysis for the independent inverse-fix queue."""
from pathlib import Path
import datetime, hashlib, importlib.util, json, time

R=Path.cwd();B=R.parent
END=datetime.datetime.fromisoformat('2026-10-22T19:16:00+08:00').timestamp()
module=B/'s22p1_complete_r2/slot/stage_report.py'
spec=importlib.util.spec_from_file_location('stage_report',module)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def put(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_suffix('.new');t.write_text(json.dumps(x,indent=2,allow_nan=False));t.replace(p)
def get(p):return json.loads(Path(p).read_text())
def main():
    (R/'shadow.lock').mkdir()
    done={};out=R/'shadow';out.mkdir()
    while time.time()<END-45:
        status=get(R/'status.json') if (R/'status.json').exists() else {}
        for rec in status.get('records',[]):
            name=rec['name']
            if name in done:continue
            reports=[]
            for f in sorted((R/name/'out').glob('p1_stage_step_*')):
                if time.time()>=END-45:break
                try:
                    q=m.summarize_step(f);put(out/name/(f.name+'.json'),q)
                    reports.append(dict(step=f.name,capture_status=q.get('capture_status'),report=str(out/name/(f.name+'.json'))))
                except Exception as e:reports.append(dict(step=f.name,status='not_evaluated',error=repr(e)))
            done[name]=dict(reports=reports,mechanism_only=rec['mechanism_only'],
                prospective=name.startswith(('prospective_','adjacent_','direct_')),
                legacy_stop=rec.get('stop_reason'),legacy_converged=rec.get('converged'),
                promoted=rec.get('promoted'),sample_role='mechanism replay' if rec['mechanism_only'] else 'independent finite target')
        put(R/'p2_shadow_progress.json',dict(completed=bool(status.get('completed')),deadline_epoch=END,
            analyzer_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),samples=done,
            thresholds='not_calibrated',acceptance_changed=False,
            limitations=['No new threshold selected','No second mesh',
                'Failure states and single-step mechanism states are not production parents',
                'Stage metrics use production fixed DOS and documented inverse eta; missing inputs remain not_evaluated']))
        if status.get('completed'):return
        time.sleep(min(30,max(0,END-time.time()-45)))
    put(R/'shadow_deadline.json',dict(deadline_epoch=END,completed=False,reason='deadline reserve',samples=list(done)))
if __name__=='__main__':main()
