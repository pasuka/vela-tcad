"""Audit both isolated first-eight-point experiments and report paired costs."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from analyze_templates_ldmos_d5_cost_controls import analyze


def summarize(root):
    out={}
    for stage,control,candidate in [('files8','baseline','files'),('efficiency8','files','efficiency')]:
        directory=root/stage;raw=(directory/'summary.json').read_bytes();batch=json.loads(raw)
        result=analyze(batch)
        if batch['points']!=8 or len(batch['cases'])!=4:raise ValueError('Unexpected screening scope')
        identities=[(c['mode'],c['gate'],c['repeat']) for c in batch['cases']]
        if set(identities)!={(m,g,0) for m in (control,candidate) for g in (4,8)}:
            raise ValueError('Missing or duplicate paired curve')
        for row,c in zip(result['rows'],batch['cases']):
            if not c.get('cold_final_audit'):raise ValueError('Missing cold final audit')
            if c['mode']=='files' and c['max_potential_difference_V']!=0.:
                raise ValueError('File reuse changed a numerical state')
            case=directory/Path(c['directory'].replace('\\','/')).name
            ledger_raw=(case/'fixed/ledger.json').read_bytes()
            if hashlib.sha256(ledger_raw).hexdigest()!=c['performance']['ledger_sha256']:
                raise ValueError('Ledger changed after performance summary')
            ledger=json.loads(ledger_raw);feedback=ledger.get('cost_feedback',[])
            row['feedback_reasons']=dict(Counter(f['reason'] for f in feedback))
            row['max_actual_step_V']=max((t['accepted_step_V'] for t in ledger['transfers']),default=0.)
            row['stage_seconds']=c['performance']['stage_seconds']
        paired=[]
        for g in (4,8):
            a=next(r for r in result['rows'] if r['gate']==g and r['mode']==control)
            b=next(r for r in result['rows'] if r['gate']==g and r['mode']==candidate)
            paired.append(dict(gate=g,wall_ratio=b['wall_seconds']/a['wall_seconds'],
                update_difference=b['updates']-a['updates'],service_difference=b['services']-a['services'],
                parent_seconds_saved=a['parent_exclusive_seconds']-b['parent_exclusive_seconds']))
        result.update(pairs=paired,summary_sha256=hashlib.sha256(raw).hexdigest(),
                      qualification='single paired screening round, not full curves or stable speedup')
        out[stage]=result
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=summarize(a.root);a.output.write_text(json.dumps(result,indent=2),encoding='utf8')
    for stage,data in result.items():print(stage,json.dumps(data['pairs']))
