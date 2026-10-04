"""Retain input-state replay failures; compare the operator's exported live state."""
from pathlib import Path
import math
import numpy as np
import analyze_simplemos_production_consistency_20260907 as audit

a=audit.a; d=audit.d; p=audit.p; OUT=audit.OUT; LOCAL=audit.LOCAL


def main():
    a.verify(OUT/'analysis_evidence.json'); files=[Path(__file__).resolve(),OUT/'analysis_evidence.json']; result=[]; worst=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask=audit.prior.support(c); root=LOCAL/c['key']/'native_referenced_joint'
        source=d.ordered(root/'state.csv',geo.count); residual=d.ordered(root/'residual.csv',geo.count)
        doping=d.ordered(p.prior.LOCAL/'vela'/c['key']/'strict/residual.csv',geo.count)
        edges=a.rows(root/'edges.csv'); effective={}
        for e in edges:
            for side in ('0','1'):
                node=int(e['node'+side]); row=tuple(float(e[k+side+unit]) for k,unit in [('psi','_V'),('electron_density','_m3'),('hole_density','_m3')])
                if node in effective:assert effective[node]==row
                effective[node]=row
        assert len(effective)==geo.count
        nativervol=geo.volumes['all_cell']*np.array([float(x.split()[1]) for x in (Path(c['ratios'])/'nodes.txt').read_text().splitlines()])
        initial=np.array([[float(x[k]) for x in source] for k in ('psi','electrons_m3','holes_m3')])
        live=np.array([effective[i] for i in range(geo.count)]).T; net=np.array([float(x['net_doping_m3']) for x in doping])
        actual=np.array([float(x['psi_residual']) for x in residual])/factor
        def replay(state):return geo.matrices['legacy']@state[0]+d.fixed.Q*(state[1]-state[2]-net)*nativervol
        old=replay(initial); new=replay(live); den=np.linalg.norm(actual[mask]); before=float(np.linalg.norm((actual-old)[mask])/den); after=float(np.linalg.norm((actual-new)[mask])/den)
        # This records a different diagnostic gate, never rewrites original rows.
        result.append(dict(key=c['key'],input_state_relative=before,live_state_relative=after,live_state_qualified=after<=1e-8,
            actual_residual_norm_C_per_m=float(den),input_remainder_norm_C_per_m=float(np.linalg.norm((actual-old)[mask])),live_remainder_norm_C_per_m=float(np.linalg.norm((actual-new)[mask])),
            max_packed_psi_change_V=float(max(abs(live[0]-initial[0]))),max_density_relative=float(max(abs(live[1:,mask]/initial[1:,mask]-1).flat))))
        order=sorted(np.flatnonzero(mask),key=lambda i:abs(actual[i]-old[i]),reverse=True)[:5]
        for i in order:worst.append(dict(key=c['key'],node_id=int(i),x_um=geo.coords[i][0],y_um=geo.coords[i][1],input_remainder_C_per_m=actual[i]-old[i],live_remainder_C_per_m=actual[i]-new[i],
            input_electron_m3=initial[1,i],live_electron_m3=live[1,i],input_hole_m3=initial[2,i],live_hole_m3=live[2,i]))
    a.write_csv(OUT/'packed_Poisson_replay.csv',result); a.write_csv(OUT/'packed_Poisson_replay_worst.csv',worst)
    files += [OUT/'packed_Poisson_replay.csv',OUT/'packed_Poisson_replay_worst.csv']
    d.matrix.freeze(OUT/'packed_Poisson_evidence.json',files); print(result,flush=True)


if __name__=='__main__':main()
