"""Lightweight, read-only-to-results status; no solver/helper imports or launch."""
import json,math
from pathlib import Path
from datetime import datetime

ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'build-release/enormal_curves_20260912'


def main():
    rows=[];pending=[]
    for path in (LOCAL/'vela').rglob('result.json'):
        try:row=json.loads(path.read_text(encoding='utf-8'))
        except (FileNotFoundError,json.JSONDecodeError):pending.append(str(path));continue
        rows.append(row)
    good={}
    for row in rows:
        if row['qualified']:good.setdefault((row['case'],int(row['index']),row['arm']),row)
    pairs=[]
    for (case,index,arm),row in good.items():
        if arm!='continuation' or (case,index,'native') not in good:continue
        other=good[(case,index,'native')]
        denominator=other['current_A_per_um']
        pairs.append(abs(row['current_A_per_um']/denominator-1) if denominator else math.inf)
    result=dict(recorded_at=datetime.now().astimezone().isoformat(),target_states=408,
        finished_attempts=len(rows),qualified_states=len(good),failed_attempts=sum(not r['qualified'] for r in rows),files_being_written=len(pending),
        paired_currents=len(pairs),max_paired_Id_relative=max(pairs,default=None),
        interpretation='Live status only; paired current comparisons do not certify dual field qualification.',
        groups=[dict(device=device,vd=vd,arm=arm,qualified=sum(r['device']==device and r['vd']==vd and r['arm']==arm for r in good.values())) for device in ('n19','n23') for vd in (.05,1.) for arm in ('continuation','native')],
        failures=[r for r in rows if not r['qualified']])
    (LOCAL/'progress.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return result

if __name__=='__main__':main()
