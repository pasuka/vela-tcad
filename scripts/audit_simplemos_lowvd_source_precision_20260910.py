"""Source-aware, output-only replay of low-Vd final rejected Newton trials."""
import argparse
import inspect
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import audit_simplemos_e1089_linesearch_20260909 as previous
import validate_simplemos_e1089_stable_merit_20260909 as m

s=m.s
LOCAL=s.REPO/'build-release/phumob_source_precision_20260910'
OUT=s.REPO/'reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910'
RUNNER=LOCAL/'runner.exe'

def build():
    s.a.verify(m.OUT/'build_evidence.json')
    # Reuse the frozen output-only instrumentation; change only its input
    # implementations to the already calibrated source and stable-merit pair.
    body=inspect.getsource(previous.build)
    replacements={
      "source=(s.REPO/'src/equation/CoupledDDAssembler.cpp').read_text(encoding='utf-8')":
      "source=(s.LOCAL/'scaled_source/CoupledDDAssembler.cpp').read_text(encoding='utf-8')",
      "source=(s.LOCAL/'restart_coordinates/NewtonSolver.cpp').read_text(encoding='utf-8')":
      "source=(m.LOCAL/'NewtonSolver.cpp').read_text(encoding='utf-8')",
      "('NewtonSolver',s.LOCAL/'restart_coordinates/compile_command.json')":
      "('NewtonSolver',m.LOCAL/'compile_command.json')",
      "s.REPO/'src/equation/CoupledDDAssembler.cpp',s.OUT/'restart_coordinates/build_evidence.json'":
      "s.LOCAL/'scaled_source/CoupledDDAssembler.cpp',m.OUT/'build_evidence.json',Path(previous.__file__)"}
    for before,after in replacements.items():
        assert body.count(before)==1,before
        body=body.replace(before,after)
    scope=dict(vars(previous),LOCAL=LOCAL,OUT=OUT,RUNNER=RUNNER,m=m,previous=previous,__file__=__file__)
    exec(compile(body,__file__,'exec'),scope)
    scope['build']()

def run():
    s.a.verify(OUT/'build_evidence.json');s.a.verify(m.OUT/'dc_evidence.json')
    cases=s.a.read(m.OUT/'contract.json')['jobs']
    results={(r['key'],r['label']):r for r in s.a.rows(m.OUT/'dc.csv') if r['kind']=='source'}
    jobs=[];files=[OUT/'build_evidence.json',m.OUT/'dc_evidence.json']
    for case in cases:
        if case['kind']!='source' or case['label'] not in ('plus_full','minus_full'):continue
        base=Path(case['dest']);dest=LOCAL/case['key']/case['label']
        cfg=s.a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv')
        s.a.write(dest/'config.json',cfg)
        jobs.append(dict(case,base=str(base),dest=str(dest),trace_iteration=int(results[case['key'],case['label']]['iterations'])+1))
        files += [dest/'config.json',Path(cfg['state_file']),Path(case['source_file'])]
    s.a.write(OUT/'contract.json',dict(jobs=jobs,scope='Replay the four +/-full low-Vd source failures from identical original inputs, with stable merit and corrected restart coordinates. Trace only the final rejected iteration. Preserve 200-iteration budget, 1e-6 row gate, 1e-8 KCL gate, damping, source and Jacobian.',production_changed=False))
    s.d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])
    def one(j):
        dest=Path(j['dest']);env=s.env(j,j['alpha'])
        env.update(VELA_STABLE_MERIT='1',VELA_LS_AUDIT_DIR=str(dest),VELA_LS_AUDIT_ITERATION=str(j['trace_iteration']))
        status=s.V.execute(dest/'config.json',RUNNER,env);old=s.a.read(Path(j['base'])/'config.status.json')
        same=s.a.sha(dest/'state.csv')==s.a.sha(Path(j['base'])/'state.csv') and all(status[k]==old[k] for k in ('iterations','exit_code','failure_reason','converged','final_residual'))
        dumps=list(dest.glob('*_direction.csv'))
        repeat=all(all(r['R']==r['repeated_R'] for r in s.a.rows(f)) for f in dumps)
        row=dict(key=j['key'],label=j['label'],identity=same,repeat_identity=repeat,snapshots=len(dumps),iterations=status['iterations'])
        print(row,flush=True);assert same and repeat and len(dumps)==14,row
        return row
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,jobs))
    s.a.write_csv(OUT/'identity.csv',rows)
    s.d.matrix.freeze(OUT/'replay_evidence.json',[OUT/'freeze.json',OUT/'identity.csv']+[p for j in jobs for p in Path(j['dest']).rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('build','run'));globals()[p.parse_args().action]()
