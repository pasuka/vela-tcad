"""The isolated selector's disabled branch must preserve existing trajectories."""
from pathlib import Path
import validate_simplemos_e1089_stable_merit_20260909 as m
s=m.s;OUT=m.OUT/'disabled_identity';LOCAL=m.LOCAL/'disabled_identity'

def run():
    s.a.verify(m.OUT/'build_evidence.json');jobs=[];files=[Path(__file__).resolve(),m.OUT/'build_evidence.json']
    for label,base in [('e1089',m.audit.LOCAL/'failed_e1089/plain'),('lowvd_source',m.w.LOCAL/'settled_source/dc/m65_n19_vd_0p050000_endpoint/plus_full')]:
        dest=LOCAL/label;cfg=s.a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv')
        if 'local_update_diagnostics' in cfg['solver']:cfg['solver']['local_update_diagnostics']['csv_file']=str(dest/'updates.csv')
        s.a.write(dest/'config.json',cfg);jobs.append(dict(label=label,base=str(base),dest=str(dest)));files += [base/'state.csv',base/'config.status.json',dest/'config.json',Path(cfg['state_file'])]
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Replay the already completed e1089 single-step and n19 low-Vd positive full-source trajectory with the new binary and selector disabled; require identical saved state, exit, iteration count, currents and final residual.'))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);rows=[]
    for j in jobs:
        dest,base=Path(j['dest']),Path(j['base']);env=s.V.environment();env.pop('VELA_STABLE_MERIT',None)
        if j['label']=='lowvd_source':
            c=next(c for c in s.a.read(m.w.OUT/'settled_source/contract.json')['cases'] if c['device']=='n19');env=s.env(c,.001)
        result=s.V.execute(dest/'config.json',m.RUNNER,env);original=s.a.read(base/'config.status.json')
        keys=('iterations','exit_code','converged','failure_reason','final_residual','contact_currents_A_per_um');same=s.a.sha(base/'state.csv')==s.a.sha(dest/'state.csv') and all(result[k]==original[k] for k in keys)
        rows.append(dict(label=j['label'],state_and_status_identity=same));assert same,j
    s.a.write_csv(OUT/'checks.csv',rows);s.d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'checks.csv']+[p for p in LOCAL.rglob('*') if p.is_file()]);print(rows,flush=True)

if __name__=='__main__':run()
