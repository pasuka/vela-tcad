"""Report the explicit Enormal candidate without combining unqualified seeds."""
import math
from pathlib import Path
import validate_simplemos_enormal_candidate_20260912 as q

a,d,O=q.a,q.d,q.O

def main():
    a.verify(O/'dc_evidence.json')
    attempts=a.rows(O/'attempts.csv');refs={r['name']:r for r in a.rows(q.p.O/'native_points.csv')}
    comparison=[];selected=[]
    for cc in a.read(O/'inputs.json')['cases']:
        group=[r for r in attempts if r['case']==cc['case'] and int(r['index'])==cc['index']]
        chosen={}
        for arm in ('phumob_seed','native_seed'):
            same=[r for r in group if r['arm']==arm];assert same
            chosen[arm]=next((r for r in same if r['qualified']=='True'),same[-1])
            selected.append(dict(key=cc['key'],**chosen[arm]))
        first,second=chosen.values();geo,mask=q.V.m.previous.prior.support(cc)
        delta=dict(psi_max_V=math.inf,phin_max_V=math.inf,phip_max_V=math.inf,density_max_relative=math.inf)
        if all((Path(r['dest'])/'state.csv').exists() for r in (first,second)):
            delta=q.V.m.previous.prior.delta_states(*[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in (first,second)],mask)
        currents=[float(r.get('current_A_per_um') or 'nan') for r in (first,second)]
        dualI=abs(currents[0]/currents[1]-1) if currents[1]!=0 else math.inf
        dual=q.v.dual_qualified(first['qualified']=='True',second['qualified']=='True',delta,dualI)
        ref=refs[cc['key']+'_baseline'];nativeI=float(ref['Id_A_per_um'])
        comparison.append(dict(key=cc['key'],device=cc['device'],vd=cc['vd'],vg=cc['vg'],native_Id_A_per_um=nativeI,
            phumob_seed_Id_A_per_um=currents[0],native_seed_Id_A_per_um=currents[1],
            error_percent=100*(currents[0]/nativeI-1),native_seed_error_percent=100*(currents[1]/nativeI-1),
            phumob_seed_qualified=first['qualified']=='True',native_seed_qualified=second['qualified']=='True',
            dual_Id_relative=dualI,**delta,dual_qualified=dual,comparison_qualified=dual and ref['qualified']=='True'))
    summary=dict(attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),
        states=len(selected),qualified_states=sum(r['qualified']=='True' for r in selected),
        points=len(comparison),dual_qualified=sum(r['dual_qualified'] for r in comparison),
        comparison_qualified=sum(r['comparison_qualified'] for r in comparison),
        max_row_ratio=max(float(r['max_row_ratio']) for r in selected),
        max_port_relative=max(float(r['port_relative']) for r in selected),
        max_dual_Id_relative=max(r['dual_Id_relative'] for r in comparison),
        max_dual_psi_V=max(r['psi_max_V'] for r in comparison),max_dual_phin_V=max(r['phin_max_V'] for r in comparison),max_dual_phip_V=max(r['phip_max_V'] for r in comparison),
        max_dual_density_relative=max(r['density_max_relative'] for r in comparison),
        acceptance_changed=False,current_error_gate='Report the original relative metric; do not invent a new absolute-current tolerance.',
        all_qualified=all(r['comparison_qualified'] for r in comparison))
    q.v.csv_union(O/'selected_states.csv',selected);a.write_csv(O/'comparison.csv',comparison);a.write(O/'summary.json',summary)
    d.matrix.freeze(O/'comparison_evidence.json',[Path(__file__).resolve(),O/'dc_evidence.json',q.p.O/'native_points.csv',O/'selected_states.csv',O/'comparison.csv',O/'summary.json'])
    print(summary,comparison,flush=True)

if __name__=='__main__':main()
