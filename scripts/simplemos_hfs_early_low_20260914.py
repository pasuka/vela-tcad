"""Run independent HFS low controls while native references are being computed.
Cache only completed attempts with input hashes and output provenance verified.
"""
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import simplemos_hfs_curves_20260914 as e
a,d,v=e.a,e.d,e.v
ROOT=e.L/'early_low'

def prepare():
    e.gate();a.verify(e.ENORMAL/'completion_evidence.json')
    assert not (e.O/'early_low_freeze.json').exists()
    source=a.rows(e.ENORMAL/'continuation_attempts.csv');cases=[];files=[Path(__file__).resolve(),e.q.O/'completion_evidence.json',e.ENORMAL/'completion_evidence.json']
    for cc in a.read(e.q.O/'inputs.json')['cases']:
        if cc['index']!=40:continue
        for i in (0,10):
            seed=next(r for r in source if r['case']==cc['case'] and int(r['index'])==i and r['qualified']=='True')
            row=dict(case=cc['case'],device=cc['device'],vd=cc['vd'],index=i,template=cc['template'],seed=str(Path(seed['dest'])/'state.csv'))
            cases.append(row);files += [Path(row['template']),Path(row['seed']),Path(seed['dest'])/'independent_acceptance.json']
    a.write(e.O/'early_low_contract.json',dict(cases=cases,gates=v.GATES,runner=str(e.q.RUNNER),runner_sha256=a.sha(e.q.RUNNER),scope='Only independent Enormal-start low controls; native qualification and full dual gate required later.',acceptance_changed=False))
    d.matrix.freeze(e.O/'early_low_freeze.json',files+[e.O/'early_low_contract.json'])

def run():
    e.configure();a.verify(e.O/'early_low_freeze.json');e.gate()
    v.LOCAL=ROOT;v.OUT=e.O/'early_low'
    cases=a.read(e.O/'early_low_contract.json')['cases']
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows=[r for group in pool.map(lambda cc:v.solve_target(cc,cc['index'],'enormal_seed',Path(cc['seed'])),cases) for r in group]
    v.csv_union(e.O/'early_low_attempts.csv',rows)
    d.matrix.freeze(e.O/'early_low_evidence.json',[e.O/'early_low_freeze.json',e.O/'early_low_attempts.csv']+[p for p in ROOT.rglob('*') if p.is_file()])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run'));globals()[p.parse_args().action]()
