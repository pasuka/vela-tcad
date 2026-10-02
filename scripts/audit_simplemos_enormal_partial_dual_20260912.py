"""Immutable snapshot of completed dual pairs; not a full-curve acceptance."""
import argparse,json,re
from pathlib import Path
import simplemos_enormal_curves_20260912 as e
a,d=e.a,e.d


def main(tag):
    assert re.fullmatch('[a-z][a-z0-9_-]*',tag)
    e.configure();a.verify(e.O/'supplemental_source_evidence.json');a.verify(e.O/'control_evidence.json')
    out=e.O/'partial_dual'/tag;assert not (out/'evidence.json').exists()
    good={}
    for path in (e.L/'vela').rglob('result.json'):
        try:row=a.read(path)
        except (FileNotFoundError,json.JSONDecodeError):continue
        if row['qualified']:good.setdefault((row['case'],int(row['index']),row['arm']),row)
    cases={c['case']:c for c in a.read(e.O/'vela_contract.json')['cases']}
    native=a.rows(e.O/'native_points.csv');results=[]
    files=[Path(__file__).resolve(),e.O/'supplemental_source_evidence.json',e.O/'control_evidence.json',e.O/'native_evidence.json']
    for (case,index,arm),first in sorted(good.items()):
        if arm!='continuation' or (case,index,'native') not in good:continue
        second=good[(case,index,'native')];geo,mask=e.v.V.m.previous.prior.support(cases[case])
        roots=[Path(r['dest']) for r in (first,second)]
        for root in roots:a.verify(root/'input_freeze.json')
        states=[d.ordered(root/'state.csv',geo.count) for root in roots]
        delta=e.v.V.m.previous.prior.delta_states(*states,mask)
        ratio=abs(first['current_A_per_um']/second['current_A_per_um']-1)
        ref=next(r for r in native if r['case']==case and int(r['index'])==index)
        qualified=e.v.dual_qualified(True,True,delta,ratio) and ref['native_qualified']=='True'
        results.append(dict(case=case,device=first['device'],vd=first['vd'],vg=first['vg'],index=index,**delta,dual_Id_relative=ratio,qualified=qualified))
        files += [root/name for root in roots for name in ('input_freeze.json','result.json','state.csv','config.status.json','independent_acceptance.json','functional.status.json','all_row.csv')]
    assert results,'No completed overlapping targets yet; no snapshot written.'
    a.write_csv(out/'comparison.csv',results)
    summary=dict(pairs=len(results),qualified=sum(r['qualified'] for r in results),all_available_pairs_qualified=all(r['qualified'] for r in results),
        max_dual_Id_relative=max(r['dual_Id_relative'] for r in results),max_psi_V=max(r['psi_max_V'] for r in results),
        max_phin_V=max(r['phin_max_V'] for r in results),max_phip_V=max(r['phip_max_V'] for r in results),max_density_relative=max(r['density_max_relative'] for r in results),
        gates=e.v.GATES,scope='Snapshot of already completed qualified pairs only. This does not certify the remaining curve targets or replace the final 204-point comparison.')
    a.write(out/'summary.json',summary);d.matrix.freeze(out/'evidence.json',files+[out/'comparison.csv',out/'summary.json'])
    print(summary,flush=True);assert summary['all_available_pairs_qualified']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--tag',required=True);args=ap.parse_args();main(args.tag)
