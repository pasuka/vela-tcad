"""Fresh tangent and signed responses after independently qualified zero plateau."""
from pathlib import Path
import copy
import validate_simplemos_restart_source_followup_20260909 as w
import complete_simplemos_phumob_local_source_20260909 as j
s=w.s;LOCAL=w.LOCAL/'settled_source';OUT=w.OUT/'settled_source'

def run():
    s.a.verify(w.OUT/'zero_plateau/evidence.json')
    assert all(r['stationary']=='True' and r['qualified']=='True' for r in s.a.rows(w.OUT/'zero_plateau/results.csv'))
    cases=[c for c in s.a.read(w.OUT/'contract.json')['cases'] if c['vd']==.05];files=[Path(__file__).resolve(),w.OUT/'zero_plateau/evidence.json',w.OUT/'freeze.json',w.RUNNER]
    for c in cases:
        base=w.LOCAL/'zero_plateau'/c['key']/'2';c['baseline']=str(base);cfg=s.a.read(base/'config.json')
        files += [base/'state.csv',base/'all_row.csv',base/'config.json',Path(c['source_file'])]
        for job in c['jobs']:
            dest=LOCAL/'dc'/c['key']/job['label'];deck=copy.deepcopy(cfg);deck.update(state_file=str(base/'state.csv'),output_state_file=str(dest/'state.csv'));s.a.write(dest/'config.json',deck);s.V.post_config(deck,dest);job['dest']=str(dest);files+=list(dest.glob('*.json'))
        for label in ('zero','unit'):
            dest=LOCAL/'fixed'/c['key']/label;files+=s.V.probes(cfg,dest,base/'state.csv');deck=s.a.read(dest/'functional.json');deck.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'));s.a.write(dest/'jacobian.json',deck);files.append(dest/'jacobian.json')
    s.a.write(OUT/'contract.json',dict(cases=cases,gates=s.a.read(w.OUT/'contract.json')['gates'],source='Identical previously frozen physical particle source; do not recompute from the polished state. Two signed amplitudes unchanged.',baseline='Second of two zero-drift restarts of a qualified zero-source state. Original historical baseline experiment remains 6/8; this is an independently frozen calibration.',fresh_full_J_tangent=True,finite_volume_replacement=False))
    scope=[r for r in s.a.rows(w.OLDOUT/'scope.csv') if r['key'] in {c['key'] for c in cases}];s.a.write_csv(OUT/'scope.csv',scope)
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',OUT/'scope.csv'])
    j.RUNNER=w.OLD/'jacobian_adapter/runner.exe';s.d.matrix.freeze(OUT/'jacobian_build_evidence.json',[w.OLDOUT/'jacobian_build_evidence.json',j.RUNNER])
    s.LOCAL=LOCAL;s.OUT=OUT;s.RUNNER=w.RUNNER
    j.preflight();assert all(r['qualified']=='True' for r in s.a.rows(OUT/'preflight.csv'))
    s.run();s.analyze()

if __name__=='__main__':run()
