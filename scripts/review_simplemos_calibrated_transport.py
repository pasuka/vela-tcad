"""Independent edge conservation, exact error budget and paired-budget review."""
import math
from pathlib import Path
import numpy as np
import decompose_simplemos_calibrated_transport as d


def main():
    d.a.verify(d.FREEZE)
    cases=d.a.rows(d.OUT/'case.csv')
    components=d.a.rows(d.OUT/'component.csv')
    pairs=d.a.rows(d.OUT/'pair.csv')
    jv=d.a.rows(d.OUT/'jvp.csv')
    assert len(cases)==4 and len(jv)==12
    assert len({(r['case'],r['step']) for r in jv})==12
    assert all(r['passed']=='True' for r in jv)
    checks=[]
    for c in d.a.read(d.CONTRACT)['cases']:
        root=d.LOCAL/c['case']
        count=c['nodes']
        nodes=d.a.rows(root/'secant.csv')
        edges=d.a.rows(root/'secant_edges.csv')
        adj=d.ordered(Path(c['adjoint']),count)
        weight=np.array([float(r['lambda_electron']) for r in adj])
        differences=[]
        for mask in range(8):
            group=[r for r in edges if int(r['mask'])==mask]
            assert len({r['edge_id'] for r in group})==len(group)
            acc=np.zeros(count,dtype=np.longdouble)
            scale=np.zeros(count,dtype=np.longdouble)
            for r in group:
                i,j=int(r['node0']),int(r['node1'])
                val=np.longdouble(r['electron_flux'])
                if r['node0_constrained']=='0':acc[i]+=val;scale[i]+=abs(val)
                if r['node1_constrained']=='0':acc[j]-=val;scale[j]+=abs(val)
            nrows=sorted([r for r in nodes if int(r['mask'])==mask],key=lambda r:int(r['node_id']))
            assert len(nrows)==count
            expected=np.array([np.longdouble(r['electron_flux']) for r in nrows])
            norm=float(np.linalg.norm(acc-expected)/max(np.linalg.norm(scale),1e-300))
            assert norm<=1e-8
            differences.append(float(np.dot(weight,acc-expected)))
        record=next(r for r in cases if r['case']==c['case'])
        budget=math.fsum(float(r['current_A_per_um']) for r in components if r['case']==c['case'])
        exact=budget+float(record['nonlinear_remainder_A_per_um'])+float(record['target_gap_A_per_um'])
        closure=abs(exact/float(record['actual_gap_A_per_um'])-1)
        assert closure<=1e-10
        assert abs(budget/float(record['predicted_state_gap_A_per_um'])-1)<=1e-10
        checks.append({'case':c['case'],'error_budget_relative':closure,
            'max_secant_edge_node_weighted_difference_A_per_um':max(map(abs,differences))})
    pair_checks=[]
    for vd in (.05,1.):
        lo,hi=[next(r for r in cases if r['device']==device and float(r['vd'])==vd) for device in ('n19','n23')]
        target=float(hi['actual_gap_A_per_um'])/float(hi['native_Id_A_per_um'])-float(lo['actual_gap_A_per_um'])/float(lo['native_Id_A_per_um'])
        budget=math.fsum(float(r['high_minus_low_relative_error']) for r in pairs if float(r['vd'])==vd)
        assert abs(budget/target-1)<=1e-10
        pair_checks.append({'vd':vd,'exact_high_minus_low_relative_error':target,'budget':budget})
    d.a.write(d.OUT/'review.json',{'checks':checks,'paired_checks':pair_checks,'passed':True,
        'input_hashes':{d.a.rel(Path(__file__).resolve()):d.a.sha(Path(__file__).resolve())}})
    print('All eight secant masks conserve edge/node flux; four exact budgets and two paired budgets close',flush=True)


if __name__=='__main__':main()
