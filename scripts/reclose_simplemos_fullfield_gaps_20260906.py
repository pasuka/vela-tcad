"""Stricter Newton stopping control for the two remaining KCL qualification gaps."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import run_simplemos_fullfield_vela_20260906 as b

a=b.a;d=b.d;OUT=b.OUT;LOCAL=b.LOCAL


def prepare():
    a.verify(OUT/'recovery_freeze.json')
    gaps=[r for r in a.rows(OUT/'vela_recovery.csv') if r['qualified']=='False']
    assert len(gaps)==2
    files=[Path(__file__).resolve(),OUT/'recovery_freeze.json',OUT/'vela_recovery.csv',b.v.RUNNER];jobs=[]
    for r in gaps:
        index=round(float(r['vg'])/.05);old=LOCAL/'native_seed'/r['case']/f'vg_{index:03d}'
        dest=LOCAL/'reclosure'/r['case']/f'vg_{index:03d}'
        cfg=a.read(old/'config.json')
        cfg['solver']['reltol']=1e-10;cfg['solver']['abstol']=1e-14
        cfg['output_state_file']=str(dest/'state.csv');a.write(dest/'config.json',cfg)
        files += [old/'config.json',Path(cfg['state_file']),dest/'config.json']
        jobs.append(str(dest/'config.json'))
    a.write(OUT/'reclosure_contract.json',dict(status='frozen_before_execution',configs=jobs,
        axis='Newton reltol 1e-7 -> 1e-10 and abstol 1e-12 -> 1e-14; same native-coherent initial state, physics and strict carrier/global/KCL acceptance.',
        previous_failures_retained=True,physical_or_production_changes=False))
    files.append(OUT/'reclosure_contract.json');d.matrix.freeze(OUT/'reclosure_freeze.json',files)


def one(path):
    row=b.qualify(Path(path).parent);row['initialization']='native_coherent_tighter_newton';print(row,flush=True);return row


def run():
    a.verify(OUT/'reclosure_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'reclosure_contract.json')['configs']))
    a.write_csv(OUT/'vela_reclosure.csv',rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));globals()[p.parse_args().action]()
