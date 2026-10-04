"""Fixed-state algebraic local row effect; not a self-consistent prediction."""
from pathlib import Path
import math
import prepare_simplemos_production_consistency_20260907 as run

a=run.a; OUT=run.OUT


def main():
    a.verify(OUT/'analysis_evidence.json'); nodes=a.rows(OUT/'local_node_ledger.csv'); edges=a.rows(OUT/'local_edge_ledger.csv'); rows=[]
    for n in nodes:
        i=int(n['node_id'])
        for car in ('electron','hole'):
            local=[e for e in edges if e['key']==n['key'] and e['carrier']==car and i in (int(e['node0']),int(e['node1']))]
            delta=math.fsum((1 if int(e['node0'])==i else -1)*float(e['fixed_native_state_mobility_change_per_m_s']) for e in local)
            convention=math.fsum((1 if int(e['node0'])==i else -1)*float(e['fixed_native_state_remaining_convention_per_m_s']) for e in local)
            before=float(n[car+'_mapped_residual_per_m_s']); after=before+delta
            rows.append(dict(key=n['key'],device=n['device'],vg=n['vg'],vd=n['vd'],node_id=i,carrier=car,incident_native_positive_edges=len(local),
                native_state_joint_residual_per_m_s=before,native_mobility_only_delta_per_m_s=delta,
                fixed_state_counterfactual_residual_per_m_s=after,residual_absolute_ratio=abs(after)/max(abs(before),1e-300),
                remaining_edge_state_convention_delta_per_m_s=convention,self_consistent=False,causal_Id_prediction=False))
    a.write_csv(OUT/'local_mobility_counterfactual.csv',rows)
    run.d.matrix.freeze(OUT/'local_counterfactual_evidence.json',[Path(__file__).resolve(),OUT/'analysis_evidence.json',OUT/'local_mobility_counterfactual.csv'])
    print([x for x in rows if x['device']=='n23' and x['vd']=='1.0' and x['node_id'] in (1000,1009)],flush=True)


if __name__=='__main__':main()
