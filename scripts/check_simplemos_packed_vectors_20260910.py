"""Compare the loaded full residual to the original in-Newton audit vector."""
from pathlib import Path
import validate_simplemos_packed_restart_v2_20260910 as v

def main():
    a,d=v.a,v.d;a.verify(v.OUT/'roundtrip_v2_evidence.json')
    audit=v.OLDOUT.parent/'packed_roundtrip/replay_evidence.json';a.verify(audit)
    original=a.rows(v.REPO/'build-release/pr_audit/replay/state_0.csv')
    base=v.LOCAL/'roundtrip_v2/failed';seed=base/'state.csv';rows=a.rows(seed);n=len(rows)
    fields=('packed_psi','packed_electron_qf_increment','packed_hole_qf_increment')
    mismatches=sum(float(r['x'])!=float(rows[int(r['node'])][fields[int(r['block'])]]) for r in original)
    assert len(original)==3*n and mismatches==0
    cfg=a.read(base/'config.json');cfg.update(simulation_type='terminal_current_functional_probe',state_file=str(seed),contact='drain',
        residual_output_csv=str(base/'vector.csv'),contact_edge_output_csv=str(base/'edges.csv'))
    cfg.pop('output_state_file',None);cfg['solver']['stable_merit_comparison']=False
    a.write(base/'vector.json',cfg);s=v.run.V.execute(base/'vector.json',v.RUNNER,v.run.V.environment());assert s['exit_code']==0
    rr=a.rows(base/'vector.csv');names=('psi_residual','phin_residual','phip_residual')
    diffs=[float(r['R'])-float(rr[int(r['node'])][names[int(r['block'])]]) for r in original]
    result=dict(coordinates=3*n,coordinate_mismatches=mismatches,residual_mismatches=sum(z!=0 for z in diffs),max_residual_difference=max(abs(z) for z in diffs))
    a.write(v.OUT/'vector_summary.json',result)
    d.matrix.freeze(v.OUT/'vector_evidence.json',[Path(__file__).resolve(),audit,v.OUT/'roundtrip_v2_evidence.json',v.OUT/'vector_summary.json',seed]+[p for p in base.glob('vector*') if p.is_file()]+[base/'edges.csv'])
    print(result,flush=True);assert result['residual_mismatches']==0

if __name__=='__main__':main()
