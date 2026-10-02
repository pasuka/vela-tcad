"""R1b: match the independently exported OldSlotboom intrinsic-density convention.

Only a copied Si.ni changes. No current data enter parameter selection.
"""
import argparse
import copy
import math
import statistics
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import simplemos_bgn_restore_vela_20260908 as w

a, d, V, v = w.a, w.d, w.V, w.v
REPO, RUNNER = w.REPO, w.RUNNER
LOCAL = w.LOCAL/'intrinsic_control'
OUT = w.OUT/'intrinsic_control'


def native_intrinsic(semantics):
    on = [r for r in semantics if r['model']=='old_slotboom']
    assert len(on)==16
    assert all(float(r['source_identity_max_V'])<=1e-12 and float(r['bgn_max_abs_error_eV'])<=1e-12 for r in on)
    bounds = [float(r[k]) for r in on for k in ('inferred_base_ni_min_cm3','inferred_base_ni_max_cm3')]
    chosen = statistics.median(bounds)
    assert max(abs(x/chosen-1) for x in bounds) <= 1e-12
    return chosen


def prepare():
    a.verify(w.OUT/'semantics_evidence.json')
    cases = a.read(w.OUT/'vela_contract.json')['cases']
    chosen = native_intrinsic(a.rows(w.OUT/'semantics.csv'))
    source = Path(a.read(Path(cases[0]['template']))['materials_file'])
    materials = a.read(source)
    si = next(m for m in materials['materials'] if m['name']=='Si')
    previous = si['ni']
    si['ni'] = chosen
    a.write(LOCAL/'materials.json', materials)
    a.write(OUT/'contract.json', dict(cases=cases, indices=w.n.INDICES, target_states=16, paths=2,
        source_ni_cm3=previous, native_derived_ni_cm3=chosen, gates=w.GATES,
        single_axis='Only copied material Si.ni; compare against R1 OldSlotboom with original no-BGN ni. Same native OldSlotboom states, mesh, mobility, BGN, SRH, constants, geometry and solver.',
        selection='Median of independently inferred native base-ni bounds from 16 states; require spatial/cross-bias constancy 1e-12. Derived from exported ni_eff, deltaEg and source Boltzmann identity; no current fitting.',
        initialization='Qualified R1 OldSlotboom Vela-start target versus native target potentials with densities recomputed at the matched base ni.',
        recovery='At most one same-bias reload; preserve all failures; no propagation between biases.',
        model_restored=False, defaults_changed=False, acceptance_changed=False))
    files = [Path(__file__).resolve(), Path(w.__file__), w.OUT/'semantics_evidence.json', w.OUT/'semantics.csv',
             w.OUT/'vela_freeze.json', LOCAL/'materials.json', source, OUT/'contract.json', RUNNER]
    d.matrix.freeze(OUT/'freeze.json', files)
    print('Frozen R1b: Si.ni', previous, '->', chosen, 'cm^-3', flush=True)


def reference_seed(c, index, arm):
    if arm == 'vela':
        root = w.LOCAL/'vela/old_slotboom'/c['case']/'vela'/f'vg_{index:03d}'
        eligible = [p for p in sorted(root.glob('attempt_*')) if a.read(p/'result.json')['qualified']]
        assert eligible, root
        return eligible[0]/'state.csv'
    dest = LOCAL/'initial'/c['case']/f'vg_{index:03d}'
    if not (dest/'freeze.json').exists():
        src = w.native_seed(c, 'old_slotboom', index)
        rows = a.rows(src)
        contract = a.read(OUT/'contract.json')
        ratio = contract['native_derived_ni_cm3']/contract['source_ni_cm3']
        for r in rows:
            r['electrons_m3'] = float(r['electrons_m3'])*ratio
            r['holes_m3'] = float(r['holes_m3'])*ratio
        a.write_csv(dest/'state.csv', rows)
        d.matrix.freeze(dest/'freeze.json', [src, dest/'state.csv', OUT/'contract.json'])
    a.verify(dest/'freeze.json')
    return dest/'state.csv'


def attempt(c, index, arm, seed, number):
    dest = LOCAL/'vela'/c['case']/arm/f'vg_{index:03d}'/f'attempt_{number}'
    cfg = a.read(w.LOCAL/'vela/old_slotboom'/c['case']/'template.json')
    for contact in cfg['contacts']:
        if contact['name']=='gate':
            contact['bias'] = w.n.prior.GRID[index]
    cfg.update(materials_file=str(LOCAL/'materials.json'), state_file=str(seed), output_state_file=str(dest/'state.csv'))
    if not (dest/'freeze.json').exists():
        a.write(dest/'config.json', cfg)
        V.post_config(cfg, dest)
        fun = copy.deepcopy(a.read(dest/'acceptance_edges.json'))
        fun.update(simulation_type='terminal_current_functional_probe', contact='drain',
                   residual_output_csv=str(dest/'port_residual.csv'), contact_edge_output_csv=str(dest/'port_edges.csv'))
        fun.pop('output_csv', None)
        a.write(dest/'functional.json', fun)
        d.matrix.freeze(dest/'freeze.json', [Path(seed), dest/'config.json', dest/'all_row.json', dest/'acceptance_edges.json', dest/'functional.json', OUT/'freeze.json'])
    else:
        assert a.read(dest/'config.json')==cfg
    a.verify(dest/'freeze.json')
    status = V.execute(dest/'config.json', RUNNER, V.environment())
    row = dict(case=c['case'], device=c['device'], vd=c['vd'], vg=w.n.prior.GRID[index], index=index,
               arm=arm, attempt=number, dest=str(dest), qualified=False, exit_code=status['exit_code'],
               failure=status.get('failure_reason',''), elapsed_seconds=status['elapsed_seconds'])
    if (dest/'state.csv').exists():
        for name in ('all_row','acceptance_edges','functional'):
            V.execute(dest/(name+'.json'), RUNNER, V.environment())
        row.update(w.old.prior.old.qualify(c, dest))
        port = a.read(dest/'functional.status.json')
        error = abs(port['current_A_per_um']/port['contact_current_extractor_A_per_um']-1)
        row.update(port_relative=error, qualified=row['qualified'] and port['exit_code']==0 and error<=w.GATES['port_relative'])
    v.write_same(dest/'result.json', row)
    print('R1b', c['device'],c['vd'],row['vg'],arm,number,row['qualified'],'row',row.get('max_row_ratio'),'Id',row.get('current_A_per_um'),flush=True)
    return row


def group(job):
    c, arm = job
    rows = []
    for index in w.n.INDICES:
        result = [attempt(c,index,arm,reference_seed(c,index,arm),0)]
        if not result[-1]['qualified'] and (Path(result[-1]['dest'])/'state.csv').exists():
            result.append(attempt(c,index,arm,Path(result[-1]['dest'])/'state.csv',1))
        rows += result
    return rows


def run():
    a.verify(OUT/'freeze.json')
    a.verify(w.OUT/'vela_vela_evidence.json')
    a.verify(w.OUT/'export_evidence.json')
    jobs = [(c,arm) for c in a.read(OUT/'contract.json')['cases'] for arm in ('vela','native')]
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [r for block in pool.map(group,jobs) for r in block]
    v.csv_union(OUT/'attempts.csv',rows)
    d.matrix.freeze(OUT/'run_evidence.json', [OUT/'freeze.json',OUT/'attempts.csv']+[p for p in LOCAL.rglob('*') if p.is_file()])


def analyze():
    a.verify(OUT/'run_evidence.json')
    attempts = a.rows(OUT/'attempts.csv')
    native = a.rows(w.OUT/'native_points.csv')
    rows = []
    for c in a.read(OUT/'contract.json')['cases']:
        geo,mask = V.m.previous.prior.support(c)
        for index in w.n.INDICES:
            ref = next(r for r in native if r['case']==c['case'] and int(r['index'])==index and r['model']=='old_slotboom')
            row = dict(case=c['case'],device=c['device'],vd=c['vd'],vg=w.n.prior.GRID[index],index=index,
                       sentaurus_Id_A_per_um=float(ref['Id_A_per_um']),sentaurus_qualified=ref['native_qualified']=='True')
            chosen = {}
            for arm in ('vela','native'):
                choices=[r for r in attempts if r['case']==c['case'] and int(r['index'])==index and r['arm']==arm]
                chosen[arm]=next((r for r in choices if r['qualified']=='True'),choices[-1])
                row[arm+'_qualified']=chosen[arm]['qualified']=='True'
                row[arm+'_Id_A_per_um']=float(chosen[arm].get('current_A_per_um') or 'nan')
                row[arm+'_error_percent']=100*(row[arm+'_Id_A_per_um']/row['sentaurus_Id_A_per_um']-1)
            delta=dict(psi_max_V=math.inf,phin_max_V=math.inf,phip_max_V=math.inf,density_max_relative=math.inf)
            if all((Path(r['dest'])/'state.csv').exists() for r in chosen.values()):
                states=[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen.values()]
                delta=V.m.previous.prior.delta_states(*states,mask)
            error=abs(row['vela_Id_A_per_um']/row['native_Id_A_per_um']-1)
            row.update(**delta,dual_Id_relative=error)
            row['dual_qualified']=v.dual_qualified(row['vela_qualified'],row['native_qualified'],delta,error)
            row['comparison_qualified']=row['dual_qualified'] and row['sentaurus_qualified']
            rows.append(row)
    v.csv_union(OUT/'comparison.csv',rows)
    summary=dict(points=len(rows),comparison_qualified=sum(r['comparison_qualified'] for r in rows),
                 attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),
                 max_abs_error_percent=max(abs(r['vela_error_percent']) for r in rows if r['comparison_qualified']))
    a.write(OUT/'summary.json',summary)
    d.matrix.freeze(OUT/'comparison_evidence.json',[OUT/'run_evidence.json',w.OUT/'native_evidence.json',OUT/'comparison.csv',OUT/'summary.json'])
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('action',choices=('prepare','run','analyze'))
    globals()[parser.parse_args().action]()
