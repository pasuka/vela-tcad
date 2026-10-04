"""Select a convergent component reference independently of the tested Jacobian."""
from pathlib import Path
import numpy as np
import validate_simplemos_weak_columns_20260906 as w

a=w.a


def richardson(coarse,half,quarter):
    r1=(4*half-coarse)/3;r2=(4*quarter-half)/3
    norm=float(np.linalg.norm(r2))
    return r2,float(np.linalg.norm(r2-r1))/max(norm,1e-300)


def main():
    a.verify(w.OUT/'freeze.json');summary=[];trials=[];files=[Path(__file__).resolve(),w.OUT/'freeze.json']
    for c in a.read(w.OUT/'contract.json')['cases']:
        path=Path(c['config']);cfg=a.read(path);raw=a.rows(path.parent/'components.csv');status=a.read(path.with_suffix('.status.json'))
        assert status['exit_code']==0 and status['component_reconstruction_scaled_error']<=1e-10
        files += [path.parent/'components.csv',path.with_suffix('.status.json')]
        for t in cfg['targets']:
            rows=[r for r in raw if (int(r['column_node']),int(r['column_block']),int(r['row_block']))==(t['node'],t['column_block'],t['row_block'])]
            ids=sorted({int(r['row_node']) for r in rows});indices={n:i for i,n in enumerate(ids)};steps={}
            for r in rows:
                step=float(r['step_V']);values=steps.setdefault(step,{k:np.zeros(len(ids)) for k in ('component_fd','full_residual_fd','jacobian_1e6','jacobian_1e7','jacobian_1e8')})
                for k,arr in values.items():arr[indices[int(r['row_node'])]]=float(r[k])
            candidates=[]
            for step in (1e-3,1e-4,1e-5,1e-6):
                arrays=[steps[s] for s in (step,step/2,step/4)]
                ref,stability=richardson(*(r['component_fd'] for r in arrays))
                direct,direct_stability=richardson(*(r['full_residual_fd'] for r in arrays))
                norm=float(np.linalg.norm(ref));passed=norm>0 and stability<=1e-4
                common=dict(case=c['case'],column_node=t['node'],column_block=t['column_block'],row_block=t['row_block'])
                row=dict(**common,coarse_step_V=step,reference_norm=norm,reference_step_relative=stability,reference_qualified=passed,
                    full_fd_step_relative=direct_stability,full_fd_vs_component=float(np.linalg.norm(direct-ref))/max(norm,1e-300))
                for label in ('1e6','1e7','1e8'):
                    jac=arrays[0]['jacobian_'+label]
                    row['jacobian_'+label+'_relative']=float(np.linalg.norm(jac-ref))/max(norm,1e-300)
                    row['jacobian_'+label+'_passed']=passed and row['jacobian_'+label+'_relative']<=1e-3
                trials.append(row);candidates.append(row)
            valid=[r for r in candidates if r['reference_qualified']]
            selected=min(valid or candidates,key=lambda r:r['reference_step_relative'])
            summary.append(dict(**selected,reference_selection='Minimum component Richardson step disagreement; no Jacobian information used.'))
    assert len(summary)==33
    a.write_csv(w.OUT/'step_trials.csv',trials);a.write_csv(w.OUT/'summary.csv',summary)
    a.write(w.OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},targets=33,
        reference_qualified=sum(r['reference_qualified'] for r in summary),
        legacy_step_passed=sum(r['jacobian_1e6_passed'] for r in summary),
        fine_step_passed=sum(r['jacobian_1e7_passed'] for r in summary),
        finer_step_passed=sum(r['jacobian_1e8_passed'] for r in summary),production_promoted=False))
    print([(r['case'],r['column_node'],r['column_block'],r['row_block'],r['reference_step_relative'],r['jacobian_1e6_relative'],r['jacobian_1e7_relative']) for r in summary],flush=True)


if __name__=='__main__':main()
