"""Exact checkpoint roundtrip plus unchanged sixteen-point dual-seed cohort."""
import argparse
import copy
import csv
from pathlib import Path
import validate_simplemos_production_weighted_20260910 as previous
a,d,run=previous.a,previous.d,previous.run
REPO=previous.REPO;OLDOUT=previous.OUT;OLDLOCAL=previous.LOCAL
LOCAL=REPO/'build-release/pr_fix';OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910'
RUNNER=previous.RUNNER
SOURCES=('src/solver/NewtonSolver.cpp','include/vela/solver/GummelSolver.h','src/io/DDSolutionCsv.cpp',
         'tests/test_dc_sweep.cpp','tests/test_newton_solver.cpp')

def configure():
    previous.LOCAL=LOCAL/'cohort';previous.OUT=OUT/'cohort';previous.configure()

def prepare():
    a.verify(OUT/'before_evidence.json')
    contract=a.read(OLDOUT/'validation_contract.json');files=[Path(__file__).resolve(),OUT/'before_evidence.json',RUNNER,REPO/'build-release/libvela_core.a']
    for j in contract['jobs']:files.extend([Path(j['original_config']),Path(j['seed'])])
    contract.update(change='Only preserve the original packed coordinates in binary128-Poisson checkpoint states and reuse validated matching hints. Original physical models, controls, seed files and one-reload protocol unchanged.')
    a.write(OUT/'cohort/validation_contract.json',contract)
    d.matrix.freeze(OUT/'cohort/validation_freeze.json',files+[REPO/p for p in SOURCES]+[OUT/'cohort/validation_contract.json'])

def roundtrip():
    a.verify(OUT/'cohort/validation_freeze.json');files=[OUT/'cohort/validation_freeze.json'];rows=[]
    controls=[('failed',OLDLOCAL/'dc/phumob/m65_n23_vd_1p000000_endpoint/native/vg_000/attempt_1'),
              ('qualified',OLDLOCAL/'dc/phumob/m65_n23_vd_1p000000_endpoint/vela/vg_000/attempt_0')]
    for name,base in controls:
        dest=LOCAL/'roundtrip'/name;cfg=a.read(base/'config.json');cfg['output_state_file']=str(dest/'state.csv');a.write(dest/'config.json',cfg)
        status=run.V.execute(dest/'config.json',RUNNER,run.V.environment());old=a.read(base/'config.status.json')
        original=list(csv.DictReader((base/'state.csv').open()));current=list(csv.DictReader((dest/'state.csv').open()))
        same=all(all(r[k]==t[k] for k in r) for r,t in zip(original,current)) and len(original)==len(current)
        assert same and all(status[k]==old[k] for k in ('iterations','final_residual','exit_code','failure_reason'))
        files += [base/'config.json',base/'state.csv']
        packed_keys=[k for k in current[0] if k.startswith('packed_')];assert len(packed_keys)==4
        legacy=dest/'legacy.csv'
        with legacy.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=[k for k in current[0] if k not in packed_keys]);writer.writeheader();writer.writerows([{k:r[k] for k in writer.fieldnames} for r in current])
        for label,seed in (('packed',dest/'state.csv'),('legacy',legacy)):
            deck=copy.deepcopy(cfg);p=dest/label;deck.update(state_file=str(seed),output_state_file=str(p/'state.csv'))
            deck['solver']['max_iter']=0;deck['solver']['carrier_row_convergence']['mode']='report'
            a.write(p/'config.json',deck);s=run.V.execute(p/'config.json',RUNNER,run.V.environment());assert s['iterations']==0
            row=dict(control=name,encoding=label,original_norm=status['final_residual'],reload_norm=s['initial_residual'],
                     same_norm=s['initial_residual']==status['final_residual'],baseline_physical_status_identical=same)
            if label=='packed':assert row['same_norm'],row
            rows.append(row)
        files += [p for p in dest.rglob('*') if p.is_file()]
    a.write_csv(OUT/'roundtrip.csv',rows);d.matrix.freeze(OUT/'roundtrip_evidence.json',files+[OUT/'roundtrip.csv']);print(rows,flush=True)

def solve():configure();previous.solve()
def analyze():configure();previous.analyze()
def jvp():configure();previous.jvp()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','roundtrip','solve','analyze','jvp'));globals()[p.parse_args().action]()
