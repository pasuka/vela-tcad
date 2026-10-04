"""Extend the qualified production combination without changing acceptance.

Two paths: continuation from a qualified Vela .8 V anchor and independent
native-coherent seeds. A failed attempt may be reloaded once at the same bias;
its original failure is retained. Only qualified outputs may seed a new bias.
"""
import argparse
import copy
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import simplemos_masetti_curve_native_20260908 as n
import validate_simplemos_generated_box_mobility_20260907 as old
import run_simplemos_m80b_state_semantics as mapping

a, d = n.a, n.d
LOCAL, OUT, REPO = n.LOCAL, n.OUT, n.REPO
V = old.v
RUNNER = REPO/'build-release/vela_example_runner.exe'
PREVIOUS = old.LOCAL/'current_extraction_fix'
GATES = dict(all_row=1e-6, kcl_over_Id=1e-8, global_tolerance=1e-6, source_floor=1e-10,
             initialization_phi_max_V=1e-6, initialization_density_relative=1e-4,
             initialization_Id_relative=1e-6, port_relative=1e-8)


def write_same(path, obj):
    if path.exists():
        assert a.read(path) == obj, path
    else:
        a.write(path, obj)


def csv_union(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    a.write_csv(path, [{k:r.get(k, '') for k in keys} for r in rows])


def prepare():
    cases = []
    files = [Path(__file__).resolve(), Path(n.__file__), RUNNER, OUT/'native_contract.json',
             REPO/'build-release/libvela_core.a']
    for c in a.read(OUT/'native_contract.json')['cases']:
        root = PREVIOUS/(c['case']+'_vg016')/'generated'
        cfg = a.read(root/'config.json')
        assert cfg['solver']['mobility'] == dict(model='masetti', doping_concentration_basis='total_impurity', edge_averaging='element_box')
        assert cfg['solver']['carrier_row_convergence']['eps_row'] == GATES['all_row']
        assert cfg['solver']['linear_refinement_iterations'] == 4
        for field in ('state_file', 'output_state_file'):
            cfg.pop(field)
        dest = LOCAL/'vela'/c['case']
        a.write(dest/'template.json', cfg)
        cases.append(dict(**c, template=str(dest/'template.json'), anchor=str(root/'state.csv')))
        files += [dest/'template.json', root/'config.json', root/'state.csv', root/'independent_acceptance.json']
        files += [Path(cfg[k]) for k in ('mesh_file', 'materials_file', 'node_doping_file')]
        assert a.read(root/'independent_acceptance.json')['result']['qualified']
    # Freeze the transitive in-repository Python helpers and production sources.
    files += list((REPO/'scripts').glob('*simplemos*.py'))
    files += [x for directory in ('src', 'include') for x in (REPO/directory).rglob('*') if x.suffix in ('.cpp', '.h')]
    a.write(OUT/'vela_contract.json', dict(cases=cases, gates=GATES, target_states=204, paths=2,
        current_metric='100*(Id_Vela/Id_Sentaurus-1); no new absolute-current acceptance threshold is invented.',
        qualification='Native qualification plus both Vela paths strict, all 1814 free Si carrier rows, independent global closure, KCL, port and dual state/current gates.',
        continuation='At .8 V reclose accepted Vela state; descend .02 V to zero and independently ascend .02 V to one from the .8 anchor. A failed target does not update the last qualified seed; later targets are still attempted.',
        independent='Each native target psi/phin/phip supplies an independent six-column coherent Vela initial state; density recalculated using unchanged Vela ni and thermal voltage.',
        recovery='At most one saved-state reload per failed attempt, same target and identical solver settings. Preserve both attempts; no unqualified seed propagates across bias.',
        coverage='Four 51-point curves. Full PhuMob/Lombardi/HFS and 16-condition expansion are separate scopes.',
        gates_changed=False, defaults_changed=False, source_changes=False))
    d.matrix.freeze(OUT/'vela_freeze.json', files+[OUT/'vela_contract.json'])
    print('Frozen 204 targets x 2 paths; at most one same-bias reload per failed solve.', flush=True)


def run_attempt(c, index, arm, seed, attempt):
    dest = LOCAL/'vela'/c['case']/arm/f'vg_{index:03d}'/f'attempt_{attempt}'
    cfg = copy.deepcopy(a.read(Path(c['template'])))
    for contact in cfg['contacts']:
        if contact['name'] == 'gate':
            contact['bias'] = n.GRID[index]
    cfg.update(state_file=str(seed), output_state_file=str(dest/'state.csv'))
    if not (dest/'input_freeze.json').exists():
        a.write(dest/'config.json', cfg)
        V.post_config(cfg, dest)
        fun = copy.deepcopy(a.read(dest/'acceptance_edges.json'))
        fun.update(simulation_type='terminal_current_functional_probe', contact='drain',
                   residual_output_csv=str(dest/'port_residual.csv'), contact_edge_output_csv=str(dest/'port_edges.csv'))
        fun.pop('output_csv', None)
        a.write(dest/'functional.json', fun)
        d.matrix.freeze(dest/'input_freeze.json', [Path(seed), dest/'config.json', dest/'all_row.json', dest/'acceptance_edges.json', dest/'functional.json'])
    else:
        assert a.read(dest/'config.json') == cfg
    a.verify(dest/'input_freeze.json')
    status = V.execute(dest/'config.json', RUNNER, V.environment())
    row = dict(case=c['case'], device=c['device'], vd=c['vd'], vg=n.GRID[index], index=index,
               arm=arm, attempt=attempt, dest=str(dest), qualified=False,
               exit_code=status['exit_code'], reason=status.get('convergence_reason', ''),
               failure=status.get('failure_reason', ''), elapsed_seconds=status['elapsed_seconds'])
    if (dest/'state.csv').exists():
        for name in ('all_row', 'acceptance_edges', 'functional'):
            V.execute(dest/(name+'.json'), RUNNER, V.environment())
        row.update(old.prior.old.qualify(c, dest))
        port = a.read(dest/'functional.status.json')
        error = abs(port['current_A_per_um']/port['contact_current_extractor_A_per_um']-1)
        row.update(port_relative=error, qualified=row['qualified'] and port['exit_code']==0 and error<=GATES['port_relative'])
    write_same(dest/'result.json', row)
    print(c['device'], c['vd'], f'{n.GRID[index]:.2f}', arm, attempt, row['qualified'],
          'row', row.get('max_row_ratio'), 'Id', row.get('current_A_per_um'), flush=True)
    return row


def solve_target(c, index, arm, seed):
    rows = [run_attempt(c, index, arm, seed, 0)]
    if not rows[-1]['qualified'] and (Path(rows[-1]['dest'])/'state.csv').exists():
        rows.append(run_attempt(c, index, arm, Path(rows[-1]['dest'])/'state.csv', 1))
    return rows


def continuation(c):
    rows = []
    seed = Path(c['anchor'])
    middle = solve_target(c, 40, 'continuation', seed)
    rows += middle
    if middle[-1]['qualified']:
        seed = Path(middle[-1]['dest'])/'state.csv'
    for indices in (range(39, -1, -1), range(41, 51)):
        current = seed
        for index in indices:
            result = solve_target(c, index, 'continuation', current)
            rows += result
            if result[-1]['qualified']:
                current = Path(result[-1]['dest'])/'state.csv'
    return rows


def native_seed(c, index):
    dest = LOCAL/'initial'/c['case']/f'vg_{index:03d}'
    path = dest/'state.csv'
    if path.exists():
        a.verify(dest/'freeze.json')
        return path
    src = LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
    geo, mask = V.m.previous.prior.support(c)
    psi, en, hp, spread = mapping.coherent_state(src, geo)
    assert spread < 1e-12
    fn = mapping.original.m73.scalar(src/'fields/eQuasiFermiPotential_region0.csv')
    fp = mapping.original.m73.scalar(src/'fields/hQuasiFermiPotential_region0.csv')
    rows = [dict(node_id=i, psi=psi[i], phin=fn.get(i, 0.), phip=fp.get(i, 0.),
                 electrons_m3=en[i], holes_m3=hp[i]) for i in range(geo.count)]
    a.write_csv(path, rows)
    d.matrix.freeze(dest/'freeze.json', [path] + list((src/'fields').glob('*.csv')))
    return path


def native_path(c):
    rows = []
    native = a.rows(OUT/'native_points.csv')
    for index in range(51):
        ref = next(r for r in native if r['case']==c['case'] and int(r['index'])==index)
        if ref['native_qualified'] != 'True':
            rows.append(dict(case=c['case'], index=index, arm='native', qualified=False, failure='native target unqualified'))
            continue
        rows += solve_target(c, index, 'native', native_seed(c, index))
    return rows


def run(arm):
    a.verify(OUT/'vela_freeze.json')
    if arm == 'native':
        a.verify(OUT/'export_evidence.json')
    cases = a.read(OUT/'vela_contract.json')['cases']
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [r for group in pool.map(continuation if arm=='continuation' else native_path, cases) for r in group]
    csv_union(OUT/(arm+'_attempts.csv'), rows)
    d.matrix.freeze(OUT/('vela_'+arm+'_evidence.json'), [OUT/'vela_freeze.json', OUT/(arm+'_attempts.csv')] +
                    [x for c in cases for x in (LOCAL/'vela'/c['case']/arm).rglob('*') if x.is_file()])


def dual_qualified(first, second, delta, current_error):
    return (first and second and math.isfinite(current_error) and current_error <= GATES['initialization_Id_relative']
            and all(math.isfinite(delta[k]) and delta[k] <= GATES['initialization_phi_max_V'] for k in ('psi_max_V', 'phin_max_V', 'phip_max_V'))
            and math.isfinite(delta['density_max_relative']) and delta['density_max_relative'] <= GATES['initialization_density_relative'])


def analyze():
    for arm in ('continuation', 'native'):
        a.verify(OUT/('vela_'+arm+'_evidence.json'))
    native = a.rows(OUT/'native_points.csv')
    attempts = {arm:a.rows(OUT/(arm+'_attempts.csv')) for arm in ('continuation', 'native')}
    comparison = []
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo, mask = V.m.previous.prior.support(c)
        for index in range(51):
            ref = next(x for x in native if x['case']==c['case'] and int(x['index'])==index)
            row = dict(case=c['case'], device=c['device'], vg=n.GRID[index], vd=c['vd'], index=index,
                       sentaurus_Id_A_per_um=float(ref['Id_A_per_um']), sentaurus_qualified=ref['native_qualified']=='True')
            chosen = {}
            for arm in ('continuation', 'native'):
                group = [x for x in attempts[arm] if x['case']==c['case'] and int(x['index'])==index]
                chosen[arm] = next((x for x in group if x['qualified']=='True'), group[-1])
                r = chosen[arm]
                row[arm+'_qualified'] = r['qualified']=='True'
                row[arm+'_Id_A_per_um'] = float(r.get('current_A_per_um') or 'nan')
                row[arm+'_error_percent'] = 100*(row[arm+'_Id_A_per_um']/row['sentaurus_Id_A_per_um']-1)
            first, second = chosen['continuation'], chosen['native']
            delta = dict(psi_max_V=math.inf, phin_max_V=math.inf, phip_max_V=math.inf, density_max_relative=math.inf)
            if all((Path(r['dest'])/'state.csv').exists() for r in (first, second)):
                states = [d.ordered(Path(r['dest'])/'state.csv', geo.count) for r in (first, second)]
                delta = V.m.previous.prior.delta_states(*states, mask)
            # Compare the two Vela currents, separately from either reference error.
            error = abs(row['continuation_Id_A_per_um']/float(second.get('current_A_per_um') or 'nan')-1)
            row.update(**delta, dual_Id_relative=error)
            row['dual_qualified'] = dual_qualified(row['continuation_qualified'], row['native_qualified'], delta, error)
            row['comparison_qualified'] = row['dual_qualified'] and ref['native_qualified']=='True'
            comparison.append(row)
    csv_union(OUT/'comparison.csv', comparison)
    summary = dict(points=len(comparison), native_qualified=sum(r['sentaurus_qualified'] for r in comparison),
                   continuation_qualified=sum(r['continuation_qualified'] for r in comparison),
                   independent_qualified=sum(r['native_qualified'] for r in comparison),
                   dual_qualified=sum(r['dual_qualified'] for r in comparison),
                   comparison_qualified=sum(r['comparison_qualified'] for r in comparison), curves=[])
    for device in ('n19', 'n23'):
        for vd in (.05, 1.):
            rows = [r for r in comparison if r['device']==device and r['vd']==vd and r['comparison_qualified']]
            worst = max(rows, key=lambda r:abs(r['continuation_error_percent'])) if rows else None
            summary['curves'].append(dict(device=device, vd=vd, qualified=len(rows), total=51,
                error_min_percent=min((r['continuation_error_percent'] for r in rows), default=None),
                error_max_percent=max((r['continuation_error_percent'] for r in rows), default=None),
                worst_vg=worst['vg'] if worst else None))
    a.write(OUT/'summary.json', summary)
    d.matrix.freeze(OUT/'comparison_evidence.json', [OUT/'vela_continuation_evidence.json', OUT/'vela_native_evidence.json', OUT/'native_evidence.json', OUT/'comparison.csv', OUT/'summary.json'])
    print(summary, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('prepare', 'continuation', 'native', 'analyze'))
    action = p.parse_args().action
    if action in ('continuation', 'native'):
        run(action)
    else:
        globals()[action]()
