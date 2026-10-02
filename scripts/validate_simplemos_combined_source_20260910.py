"""Four-control source calibration with packed Poisson and exact contact rows."""
import argparse
from pathlib import Path
import validate_simplemos_exact_contact_step_20260910 as contact
import calibrate_simplemos_packed_poisson_response_20260910 as old
s=contact.s;w=old.w
LOCAL=contact.LOCAL/'full_source';OUT=contact.OUT/'full_source'
ORIGINAL_NATIVE=s.OUT/'native/comparison'

def bind():
    def execute(path,case,alpha):
        env=s.env(case,alpha);env.update(VELA_STABLE_MERIT='1',VELA_POISSON_PRECISION='packed',VELA_EXACT_CONTACT_STEP='1')
        return s.V.execute(path,old.RUNNER if path.stem=='jacobian' else contact.RUNNER,env)
    s.LOCAL=LOCAL;s.OUT=OUT;s.RUNNER=contact.RUNNER;s.execute=execute

def prepare():
    s.a.verify(contact.OUT/'dc_evidence.json');s.a.verify(old.OUT/'preflight_evidence.json')
    cases=s.a.read(w.OUT/'contract.json')['cases'];settled={c['key']:c for c in s.a.read(w.OUT/'settled_source/contract.json')['cases']}
    files=[Path(__file__).resolve(),contact.OUT/'dc_evidence.json',old.OUT/'preflight_evidence.json',w.OUT/'freeze.json',w.OUT/'settled_source/freeze.json',old.RUNNER,contact.RUNNER]
    for case in cases:
        if case['key'] in settled:case.update(settled[case['key']])
        base=Path(case['baseline']);files += [base/'state.csv',base/'all_row.csv',Path(case['source_file'])]
        for job in case['jobs']:
            original=Path(job['dest']);dest=LOCAL/'dc'/case['key']/job['label'];cfg=s.a.read(original/'config.json');cfg['output_state_file']=str(dest/'state.csv')
            assert Path(cfg['state_file'])==base/'state.csv'
            s.a.write(dest/'config.json',cfg);s.V.post_config(cfg,dest);job['dest']=str(dest);files+=list(dest.glob('*.json'))
        cfg=s.a.read(Path(case['jobs'][0]['dest'])/'config.json')
        for label in ('zero','unit'):
            dest=LOCAL/'fixed'/case['key']/label;files+=s.V.probes(cfg,dest,base/'state.csv');deck=s.a.read(dest/'functional.json');deck.update(simulation_type='parameter_jacobian',output_csv=str(dest/'jacobian.csv'));s.a.write(dest/'jacobian.json',deck);files.append(dest/'jacobian.json')
    s.a.write(OUT/'contract.json',dict(cases=cases,gates=s.a.read(w.OUT/'contract.json')['gates'],scope='Combined exact contact step, stable merit, restart-coordinate repair and packed Poisson. Four Vg=0 controls, original sources and +/-0.001/0.0005 amplitudes, fresh full-J preflight, unchanged 200-iteration and physical gates. Diagnostic implementation only.'))
    s.a.write_csv(OUT/'scope.csv',s.a.rows(w.OLDOUT/'scope.csv'));s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json',OUT/'scope.csv'])
    bind();s.preflight()

def run():bind();s.run()

def analyze():
    bind();s.analyze();rows=s.a.rows(OUT/'ports.csv');native=s.a.rows(ORIGINAL_NATIVE/'terminal_derivatives.csv');dc=s.a.rows(OUT/'dc.csv');comp=[]
    for row in rows:
        if row['port'] not in ('drain','substrate'):continue
        target=next(r for r in native if all(r[k]==row[k] for k in ('key','amplitude','port')));actual=float(row['derivative_A_per_um']);expected=float(target['derivative_A_per_um']);error=abs(actual/expected-1)
        qualified=all(r['qualified']=='True' for r in dc if r['key']==row['key']) and error<=1e-3
        comp.append(dict(key=row['key'],amplitude=row['amplitude'],port=row['port'],vela_A_per_um=actual,native_A_per_um=expected,relative=error,qualified=qualified))
    s.a.write_csv(OUT/'native_ports.csv',comp);summary=s.a.read(OUT/'summary.json');summary.update(native_port_responses=len(comp),qualified_native_ports=sum(r['qualified'] for r in comp),max_native_relative=max(r['relative'] for r in comp),numerical_combination='packed Poisson + exact contact step + stable merit + corrected restart coordinates',production_modified=False)
    s.a.write(OUT/'combined_summary.json',summary);s.d.matrix.freeze(OUT/'combined_evidence.json',[OUT/'evidence.json',ORIGINAL_NATIVE/'evidence.json',OUT/'native_ports.csv',OUT/'combined_summary.json']);print(summary,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run','analyze'));globals()[p.parse_args().action]()
