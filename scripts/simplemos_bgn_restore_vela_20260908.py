"""Paired OldSlotboom restoration with unchanged production geometry and gates.

This is a bounded model restoration experiment, not promotion of a new default.
"""
import argparse
import copy
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import simplemos_bgn_restore_native_20260908 as n
import simplemos_masetti_curve_vela_20260908 as v

a, d, V, old, mapping = v.a, v.d, v.V, v.old, v.mapping
LOCAL, OUT, REPO = n.LOCAL, n.OUT, n.REPO
GATES = dict(v.GATES)
RUNNER = v.RUNNER
VT = 1.380649e-23 * 300 / 1.602176634e-19


def delta_eg(total_cm3, model):
    if model == 'no_bgn' or total_cm3 <= 0:
        return 0.
    assert model == 'old_slotboom'
    x = math.log(total_cm3 / 1e17)
    return .009 * (x + math.sqrt(x*x + .5))


def model_config(template, model):
    assert model in n.MODELS
    cfg = copy.deepcopy(template)
    assert cfg['solver']['mobility'] == dict(model='masetti', doping_concentration_basis='total_impurity', edge_averaging='element_box')
    assert cfg['solver']['bandgap_narrowing']['model'] == 'none'
    cfg['solver']['bandgap_narrowing']['model'] = 'old_slotboom' if model == 'old_slotboom' else 'none'
    return cfg


def prepare():
    assert not (OUT/'vela_freeze.json').exists()
    cases = a.read(OUT/'native_contract.json')['cases']
    files = [RUNNER, REPO/'build-release/libvela_core.a', OUT/'native_contract.json']
    for c in cases:
        template = v.LOCAL/'vela'/c['case']/'template.json'
        cfg = a.read(template)
        files += [template] + [Path(cfg[k]) for k in ('mesh_file', 'node_doping_file', 'materials_file')]
        c['template'] = str(template)
        for i in n.INDICES:
            seed = v.LOCAL/'vela'/c['case']/'continuation'/f'vg_{i:03d}'/'attempt_0'
            assert a.read(seed/'result.json')['qualified']
            files += [seed/'result.json', seed/'state.csv', seed/'independent_acceptance.json']
        for model in n.MODELS:
            path = LOCAL/'vela'/model/c['case']/'template.json'
            a.write(path, model_config(cfg, model))
            files.append(path)
    files += list((REPO/'scripts').glob('*simplemos*.py'))
    files += [x for directory in ('src', 'include') for x in (REPO/directory).rglob('*') if x.suffix in ('.h', '.cpp')]
    a.write(OUT/'vela_contract.json', dict(cases=cases, target_states=32, paths=2, gates=GATES,
        single_axis='Only solver.bandgap_narrowing.model; all material parameters, constants, geometry, mobility, SRH volume and solver settings are retained.',
        paths_description='Each target independently starts from its qualified no-BGN Vela state, or corresponding native target potentials with Vela-consistent BGN densities.',
        recovery='At most one same-bias saved-state reload on failure; preserve original attempt. No failed state propagates to another target.',
        current_metric='100*(Id_Vela/Id_Sentaurus-1); report physical BGN increment in each solver separately. No absolute-current acceptance threshold invented.',
        semantics=dict(delta_eg_max_abs_eV=1e-12, native_identity_max_V=1e-12,
                       base_ni_relative=1e-5, note='Model parity gate only; never retune ni or constants.'),
        jvp=dict(states='Both models, four n19/n23 x Vd cases at Vg=.8; first qualified Vela-start state.',
                 directions='Every third free Si node, separately psi/phin/phip, three amplitudes 1e-4/3e-5/1e-5 V.',
                 metric='True per-block norm error divided by max(analytic, finite-difference), undo built-in max(1,FD) floor.',
                 strong_block_relative=1e-4, qualification='Two smaller amplitudes must pass strong nonzero blocks. Weak cross-carrier source blocks reported separately; this is not a full Jacobian audit.'),
        defaults_changed=False, gates_changed=False, source_changes=False))
    d.matrix.freeze(OUT/'vela_freeze.json', files+[OUT/'vela_contract.json'])
    print('Frozen 32 configurations x 2 initializations.', flush=True)


def run_attempt(c, model, index, arm, seed, attempt):
    dest = LOCAL/'vela'/model/c['case']/arm/f'vg_{index:03d}'/f'attempt_{attempt}'
    cfg = a.read(LOCAL/'vela'/model/c['case']/'template.json')
    for contact in cfg['contacts']:
        if contact['name'] == 'gate':
            contact['bias'] = n.prior.GRID[index]
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
    row = dict(case=c['case'], device=c['device'], vd=c['vd'], model=model, vg=n.prior.GRID[index],
               index=index, arm=arm, attempt=attempt, dest=str(dest), qualified=False,
               exit_code=status['exit_code'], reason=status.get('convergence_reason', ''),
               failure=status.get('failure_reason', ''), elapsed_seconds=status['elapsed_seconds'])
    if (dest/'state.csv').exists():
        for name in ('all_row', 'acceptance_edges', 'functional'):
            V.execute(dest/(name+'.json'), RUNNER, V.environment())
        row.update(old.prior.old.qualify(c, dest))
        port = a.read(dest/'functional.status.json')
        error = abs(port['current_A_per_um']/port['contact_current_extractor_A_per_um']-1)
        row.update(port_relative=error, qualified=row['qualified'] and port['exit_code']==0 and error<=GATES['port_relative'])
    v.write_same(dest/'result.json', row)
    print(model, c['device'], c['vd'], row['vg'], arm, attempt, row['qualified'],
          'row', row.get('max_row_ratio'), 'Id', row.get('current_A_per_um'), flush=True)
    return row


def native_seed(c, model, index):
    dest = LOCAL/'initial'/model/c['case']/f'vg_{index:03d}'
    if (dest/'freeze.json').exists():
        a.verify(dest/'freeze.json')
        return dest/'state.csv'
    src = LOCAL/'native_exports'/model/c['case']/f'vg_{index:03d}'
    geo, mask = V.m.previous.prior.support(c)
    psi, en, hp, spread = mapping.RAW_READER(src, geo)
    assert spread < 1e-12
    fn = mapping.original.m73.scalar(src/'fields/eQuasiFermiPotential_region0.csv')
    fp = mapping.original.m73.scalar(src/'fields/hQuasiFermiPotential_region0.csv')
    cfg = a.read(Path(c['template']))
    doping = {int(r['node_id']):float(r['donors_cm3'])+float(r['acceptors_cm3']) for r in a.rows(Path(cfg['node_doping_file']))}
    ni = next(m['ni'] for m in a.read(Path(cfg['materials_file']))['materials'] if m['name']=='Si')*1e6
    for i in fn:
        effective = ni*math.exp(delta_eg(doping[i], model)/(2*VT))
        en[i] = effective*math.exp((psi[i]-fn[i])/VT)
        hp[i] = effective*math.exp((fp[i]-psi[i])/VT)
    rows = [dict(node_id=i, psi=psi[i], phin=fn.get(i, 0.), phip=fp.get(i, 0.),
                 electrons_m3=en[i], holes_m3=hp[i]) for i in range(geo.count)]
    a.write_csv(dest/'state.csv', rows)
    d.matrix.freeze(dest/'freeze.json', [dest/'state.csv', Path(cfg['materials_file']), Path(cfg['node_doping_file'])]+list((src/'fields').glob('*.csv')))
    return dest/'state.csv'


def run_group(job):
    c, model, arm = job
    rows = []
    references = a.rows(OUT/'native_points.csv') if arm == 'native' else []
    for index in n.INDICES:
        if arm == 'native':
            ref = next(r for r in references if r['case']==c['case'] and r['model']==model and int(r['index'])==index)
            if ref['native_qualified'] != 'True':
                raise ValueError('Cannot seed from unqualified native target: '+str(ref))
            seed = native_seed(c, model, index)
        else:
            seed = v.LOCAL/'vela'/c['case']/'continuation'/f'vg_{index:03d}'/'attempt_0/state.csv'
        group = [run_attempt(c, model, index, arm, seed, 0)]
        if not group[-1]['qualified'] and (Path(group[-1]['dest'])/'state.csv').exists():
            group.append(run_attempt(c, model, index, arm, Path(group[-1]['dest'])/'state.csv', 1))
        rows += group
    return rows


def run(arm):
    a.verify(OUT/'vela_freeze.json')
    if arm == 'native':
        a.verify(OUT/'export_evidence.json')
    jobs = [(c, m, arm) for m in n.MODELS for c in a.read(OUT/'vela_contract.json')['cases']]
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = [r for group in pool.map(run_group, jobs) for r in group]
    v.csv_union(OUT/(arm+'_attempts.csv'), rows)
    d.matrix.freeze(OUT/('vela_'+arm+'_evidence.json'), [OUT/'vela_freeze.json', OUT/(arm+'_attempts.csv')]+
                    [x for c,m,_ in jobs for x in (LOCAL/'vela'/m/c['case']/arm).rglob('*') if x.is_file()])


def semantics():
    a.verify(OUT/'export_evidence.json')
    rows = []
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo, mask = V.m.previous.prior.support(c)
        cfg = a.read(Path(c['template']))
        base = next(m['ni'] for m in a.read(Path(cfg['materials_file']))['materials'] if m['name']=='Si')
        for model in n.MODELS:
            log = (LOCAL/'native_raw/bundle'/model/c['case']/'console.log').read_text(errors='replace')
            model_lines = [line.strip() for line in log.splitlines() if 'Bandgap narrowing model:' in line]
            for index in n.INDICES:
                src = LOCAL/'native_exports'/model/c['case']/f'vg_{index:03d}'
                fields = lambda name: mapping.original.m73.scalar(src/'fields'/(name+'_region0.csv'))
                psi, en, hp, spread = mapping.RAW_READER(src, geo)
                fn, fp = fields('eQuasiFermiPotential'), fields('hQuasiFermiPotential')
                narrow, ni, nd, na = [fields(name) for name in ('BandgapNarrowing', 'EffectiveIntrinsicDensity', 'DonorConcentration', 'AcceptorConcentration')]
                ids = np.array(sorted(fn))
                x = np.concatenate(([math.log(en[i]/(ni[i]*1e6)) for i in ids], [math.log(hp[i]/(ni[i]*1e6)) for i in ids]))
                y = np.concatenate(([psi[i]-fn[i] for i in ids], [fp[i]-psi[i] for i in ids]))
                vt = float(np.dot(x,y)/np.dot(x,x))
                identity = float(np.max(np.abs(y-vt*x)))
                predicted = {i:delta_eg(nd[i]+na[i], model) for i in ids}
                bgn_error = max(abs(narrow[i]-predicted[i]) for i in ids)
                inferred = np.array([ni[i]*math.exp(-narrow[i]/(2*vt)) for i in ids])
                ni_error = float(np.max(np.abs(inferred/base-1)))
                runtime_ok = len(model_lines)==1 and (('OldSlotboom' in model_lines[0]) if model=='old_slotboom' else ('without bandgap narrowing' in model_lines[0]))
                rows.append(dict(case=c['case'], device=c['device'], vd=c['vd'], model=model,
                    index=index, vg=n.prior.GRID[index], nodes=len(ids), runtime_model=' | '.join(model_lines),
                    native_VT=vt, source_identity_max_V=identity, bgn_max_abs_error_eV=bgn_error,
                    bgn_max_eV=max(narrow.values()), inferred_base_ni_min_cm3=float(np.min(inferred)),
                    inferred_base_ni_max_cm3=float(np.max(inferred)), vela_base_ni_cm3=base,
                    base_ni_max_relative=ni_error, geometry_spread=spread,
                    model_semantics_qualified=runtime_ok and spread<1e-12 and identity<=1e-12 and bgn_error<=1e-12 and ni_error<=1e-5))
    v.csv_union(OUT/'semantics.csv', rows)
    d.matrix.freeze(OUT/'semantics_evidence.json', [OUT/'export_evidence.json', OUT/'semantics.csv', OUT/'vela_freeze.json'])
    print('Model semantics', sum(r['model_semantics_qualified'] for r in rows), '/', len(rows), flush=True)


def jvp():
    a.verify(OUT/'vela_vela_evidence.json')
    rows, files = [], []
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo, mask = V.m.previous.prior.support(c)
        ids = np.where(mask)[0][::3].tolist()
        for model in n.MODELS:
            attempts = [LOCAL/'vela'/model/c['case']/'vela/vg_040'/f'attempt_{i}' for i in (0,1)]
            eligible = [p for p in attempts if (p/'result.json').exists() and a.read(p/'result.json')['qualified']]
            if not eligible:
                rows.append(dict(case=c['case'], model=model, qualification='no_qualified_state'))
                continue
            src = eligible[0]
            dest = LOCAL/'jvp'/model/c['case']
            cfg = a.read(src/'config.json')
            cfg.update(simulation_type='newton_jvp_probe', state_file=str(src/'state.csv'), output_csv=str(dest/'jvp.csv'),
                directions=[dict(name=f'{mode}_{h:.0e}', mode=mode, amplitude_V=h, node_ids=ids, exclude_contacts=True)
                            for mode in ('psi', 'phin', 'phip') for h in (1e-4,3e-5,1e-5)])
            cfg.pop('output_state_file', None)
            v.write_same(dest/'config.json', cfg)
            if not (dest/'freeze.json').exists():
                d.matrix.freeze(dest/'freeze.json', [dest/'config.json', src/'state.csv', OUT/'vela_freeze.json'])
            a.verify(dest/'freeze.json')
            status = V.execute(dest/'config.json', RUNNER, V.environment())
            assert status['exit_code']==0, status
            for row in a.rows(dest/'jvp.csv'):
                for block in ('psi', 'phin', 'phip'):
                    ana, fd = float(row[f'analytic_{block}_norm']), float(row[f'finite_difference_{block}_norm'])
                    absolute = float(row[f'{block}_relative_error'])*max(1., fd)
                    error = absolute/max(ana, fd, 1e-300)
                    weak = (row['mode'], block) in (('phin','phip'), ('phip','phin'))
                    rows.append(dict(case=c['case'], model=model, direction=row['direction'], input_block=row['mode'],
                        output_block=block, amplitude_V=float(row['amplitude_V']), analytic_norm=ana,
                        fd_norm=fd, absolute_error=absolute, true_relative=error,
                        weak_carrier_cross_block=weak, qualified=error<=1e-4,
                        qualification='source_cross_block_needs_independent_precision' if weak else 'strong_block'))
            files += [p for p in dest.rglob('*') if p.is_file()]
    v.csv_union(OUT/'jvp.csv', rows)
    d.matrix.freeze(OUT/'jvp_evidence.json', files+[OUT/'jvp.csv', OUT/'vela_vela_evidence.json'])


def analyze():
    for name in ('vela_vela_evidence', 'vela_native_evidence', 'native_evidence', 'semantics_evidence', 'jvp_evidence'):
        a.verify(OUT/(name+'.json'))
    native = a.rows(OUT/'native_points.csv')
    semantic = a.rows(OUT/'semantics.csv')
    attempts = {arm:a.rows(OUT/(arm+'_attempts.csv')) for arm in ('vela','native')}
    rows = []
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo, mask = V.m.previous.prior.support(c)
        for model in n.MODELS:
            for index in n.INDICES:
                match = lambda r: r['case']==c['case'] and r['model']==model and int(r['index'])==index
                ref, sem = next(r for r in native if match(r)), next(r for r in semantic if match(r))
                selected = {}
                row = dict(case=c['case'], device=c['device'], vd=c['vd'], vg=n.prior.GRID[index], model=model,
                           index=index, sentaurus_Id_A_per_um=float(ref['Id_A_per_um']), sentaurus_qualified=ref['native_qualified']=='True',
                           model_semantics_qualified=sem['model_semantics_qualified']=='True')
                for arm in ('vela','native'):
                    group = [r for r in attempts[arm] if match(r)]
                    selected[arm] = next((r for r in group if r['qualified']=='True'), group[-1])
                    r = selected[arm]
                    row[arm+'_qualified'] = r['qualified']=='True'
                    row[arm+'_Id_A_per_um'] = float(r.get('current_A_per_um') or 'nan')
                    row[arm+'_error_percent'] = 100*(row[arm+'_Id_A_per_um']/row['sentaurus_Id_A_per_um']-1)
                delta = dict(psi_max_V=math.inf, phin_max_V=math.inf, phip_max_V=math.inf, density_max_relative=math.inf)
                if all((Path(r['dest'])/'state.csv').exists() for r in selected.values()):
                    states = [d.ordered(Path(r['dest'])/'state.csv', geo.count) for r in selected.values()]
                    delta = V.m.previous.prior.delta_states(*states, mask)
                error = abs(row['vela_Id_A_per_um']/row['native_Id_A_per_um']-1)
                row.update(**delta, dual_Id_relative=error)
                row['dual_qualified'] = v.dual_qualified(row['vela_qualified'], row['native_qualified'], delta, error)
                row['comparison_qualified'] = row['dual_qualified'] and row['sentaurus_qualified'] and row['model_semantics_qualified']
                rows.append(row)
    v.csv_union(OUT/'comparison.csv', rows)
    increments = []
    for off in [r for r in rows if r['model']=='no_bgn']:
        on = next(r for r in rows if r['model']=='old_slotboom' and r['case']==off['case'] and r['index']==off['index'])
        increments.append(dict(device=off['device'], vd=off['vd'], vg=off['vg'],
            pair_qualified=off['comparison_qualified'] and on['comparison_qualified'],
            native_bgn_change_percent=100*(on['sentaurus_Id_A_per_um']/off['sentaurus_Id_A_per_um']-1),
            vela_bgn_change_percent=100*(on['vela_Id_A_per_um']/off['vela_Id_A_per_um']-1),
            no_bgn_error_percent=off['vela_error_percent'], bgn_error_percent=on['vela_error_percent']))
    v.csv_union(OUT/'increments.csv', increments)
    summary = dict(points=len(rows), native_dc_qualified=sum(r['sentaurus_qualified'] for r in rows),
        semantic_qualified=sum(r['model_semantics_qualified'] for r in rows),
        dual_qualified=sum(r['dual_qualified'] for r in rows), comparison_qualified=sum(r['comparison_qualified'] for r in rows),
        attempts=sum(len(r) for r in attempts.values()), failed_attempts=sum(r['qualified']!='True' for arm in attempts.values() for r in arm),
        qualified_pairs=sum(r['pair_qualified'] for r in increments))
    a.write(OUT/'summary.json', summary)
    d.matrix.freeze(OUT/'comparison_evidence.json', [OUT/(name+'.json') for name in ('vela_vela_evidence','vela_native_evidence','native_evidence','semantics_evidence','jvp_evidence')]+
                    [OUT/'comparison.csv', OUT/'increments.csv', OUT/'summary.json'])
    print(summary, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('prepare','vela','native','semantics','jvp','analyze'))
    action = parser.parse_args().action
    run(action) if action in ('vela','native') else globals()[action]()
