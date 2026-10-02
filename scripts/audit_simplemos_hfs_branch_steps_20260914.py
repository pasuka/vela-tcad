"""Read-only geometric branch audit for the frozen real-state Jv directions.

Use independent native triangle geometry and raw physical potentials. This
classifies finite perturbations; it does not replace production Jv checks.
"""
from pathlib import Path
import numpy as np
import validate_simplemos_hfs_candidate_20260914 as h
import calibrate_simplemos_hfs_formation_v2_20260914 as f

def main():
    rows=[]
    for cc in h.a.read(h.O/'inputs.json')['cases']:
        g,contacts,boundary=f.geometry(cc['device']);state=h.a.rows(Path(cc['native_seed']))
        fields={name:{int(r['node_id']):float(r[name]) for r in state} for name in ('psi','phin','phip')}
        for label,selection in [('single',{1000}),('hotspots',{320,792,1000,1009,1057,1089,1091})]:
            ids={cid for n in selection for cid in g['nodecells'].get(n,[])}
            for mode in ('psi','phin','phip'):
                for cid in ids:
                    ns=g['cells'][cid]['nodes'];grad=g['gradient'][cid]*1e4
                    dg=np.array([float(n in selection)-float(ns[0] in selection) for n in ns[1:]])@grad
                    gradient={name:np.array([fields[name][n]-fields[name][ns[0]] for n in ns[1:]])@grad for name in fields}
                    contact=any(n in contacts for n in ns);proj=np.eye(2)
                    if boundary[cid]:
                        proj=boundary[cid][0][1]
                        if any(np.linalg.norm(p-proj)>1e-10 for _,p in boundary[cid][1:]):proj=np.zeros((2,2))
                    en0=gradient['psi']@g['gd'][cid]
                    for car,qf in [('e','phin'),('h','phip')]:
                        raw0=gradient['psi'] if contact else proj@gradient[qf]
                        change=dg if contact and mode=='psi' else proj@dg if not contact and mode==qf else np.zeros(2)
                        for step in (1e-6,1e-10,1e-12):
                            values=[float(np.linalg.norm(raw0+sgn*step*change)) for sgn in (0,1,-1)]
                            en=[float(en0+sgn*step*(dg@g['gd'][cid] if mode=='psi' else 0.)) for sgn in (0,1,-1)]
                            rows.append(dict(key=cc['key'],direction=label+'_'+mode,cell=cid,nodes=' '.join(map(str,ns)),carrier=car,step_V=step,contact=contact,
                                F0_V_cm=values[0],Fplus_V_cm=values[1],Fminus_V_cm=values[2],hfs_cross=len({x>=1 for x in values})>1,
                                Enormal0_signed_V_cm=en[0],Enormal_plus_signed_V_cm=en[1],Enormal_minus_signed_V_cm=en[2],enormal_sign_cross=min(en)<0<max(en)))
    h.a.write_csv(h.O/'branch_steps.csv',rows)
    summary=[dict(step_V=step,samples=sum(r['step_V']==step for r in rows),hfs_cross=sum(r['step_V']==step and r['hfs_cross'] for r in rows),enormal_sign_cross=sum(r['step_V']==step and r['enormal_sign_cross'] for r in rows)) for step in (1e-6,1e-10,1e-12)]
    h.a.write(h.O/'branch_steps_summary.json',summary)
    h.d.matrix.freeze(h.O/'branch_steps_evidence.json',[Path(__file__).resolve(),Path(f.__file__),h.O/'input_evidence.json',h.O/'branch_steps.csv',h.O/'branch_steps_summary.json'])
    print(summary,flush=True)

if __name__=='__main__':main()
