"""Overlap completed native curves with identical frozen local solve jobs.

Staged raw files must match the final whole-run archive before final reporting.
No incomplete curve, changed model, or alternative acceptance is permitted.
"""
import argparse
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import simplemos_masetti_curve_vela_20260908 as v


def main(case):
    a, d, n = v.a, v.d, v.n
    a.verify(v.OUT/'vela_freeze.json')
    a.verify(v.OUT/'native_freeze.json')
    c = next(c for c in a.read(v.OUT/'vela_contract.json')['cases'] if c['case']==case)
    raw = v.LOCAL/'stage'/case
    assert int((raw/'exit_code.txt').read_text())==0
    log = (raw/'console.log').read_text(errors='replace')
    assert 'T-2022.03-SP2' in log and 'Good Bye' in log
    for p in (v.LOCAL/'bundle'/case).iterdir():
        assert a.sha(p)==a.sha(raw/p.name)
    jobs, refs = [], []
    for index, vg in enumerate(n.GRID):
        rows = n.exporter.pltrows(raw/f'point_{index:03d}_native_des.plt')
        assert len(rows)==1
        r = rows[0]
        assert abs(r['gate OuterVoltage']-vg)<=1e-10 and abs(r['drain OuterVoltage']-c['vd'])<=1e-10
        cc = [r[x+' TotalCurrent'] for x in ('drain', 'source', 'gate', 'substrate')]
        assert abs(math.fsum(cc))/abs(cc[0])<=1e-8
        refs.append(dict(index=index, vg=vg, Id=cc[0]))
        jobs.append(dict(case=case, index=index, tdr=str(raw/f'vg_{index:03d}_des.tdr'),
                         export=str(v.LOCAL/'native_exports'/case/f'vg_{index:03d}')))
    out = v.OUT/'staged'/case
    a.write(out/'contract.json', dict(case=case, target_states=51, refs=refs,
        method='Run identical frozen native-initialization target jobs from a completed curve snapshot while other native curves run. Require raw byte identity to final archive before final report.'))
    d.matrix.freeze(out/'freeze.json', [Path(__file__).resolve(), v.OUT/'vela_freeze.json', out/'contract.json']+
                    [x for x in raw.iterdir() if x.is_file()])
    with ThreadPoolExecutor(max_workers=2) as pool:
        list(pool.map(n.exporter.export_one, jobs))
    geo, _ = v.V.m.previous.prior.support(c)
    for job in jobs:
        assert v.mapping.original.m73.coordinate_error(Path(job['export']), geo)<=1e-12
    def one(index):
        return v.solve_target(c, index, 'native', v.native_seed(c, index))
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = [r for group in pool.map(one, range(51)) for r in group]
    v.csv_union(out/'attempts.csv', results)
    d.matrix.freeze(out/'evidence.json', [out/'freeze.json', out/'attempts.csv']+
                    [x for x in (v.LOCAL/'vela'/case/'native').rglob('*') if x.is_file()])
    print('Staged curve complete:', case, sum(r['qualified'] for r in results), '/', len(results), flush=True)


if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case')
    main(p.parse_args().case)
