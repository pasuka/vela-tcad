"""Read-only four-point residual, transport-factor and geometry screening ledger."""
import argparse
import copy
from decimal import Decimal
import math
from pathlib import Path
import numpy as np
import validate_simplemos_vector_chain_matrix as matrix
import run_simplemos_fixed_state_formula_audit as fixed
import run_simplemos_m16_transport_factor as factors

a = matrix.a
REPO, ROOT = matrix.REPO, matrix.p.ROOT
LOCAL = REPO / 'build-release/simplemos_calibrated_transport_20260906'
OUT = ROOT / 'calibrated_transport_decomposition'
DOC = REPO / 'docs/validation/simplemos_calibrated_transport_decomposition_2026-09-06.md'
CONTRACT, FREEZE = OUT / 'contract.json', OUT / 'freeze.json'
COLS = fixed.upstream.original.m15.COMPONENT_COLUMNS
SCALE = .001


def ordered(path, count):
    rows = a.rows(path)
    if len(rows) != count or {int(r['node_id']) for r in rows} != set(range(count)):
        raise ValueError(f'Incomplete node support: {path}')
    return sorted(rows, key=lambda r: int(r['node_id']))


def array(rows, cols):
    value = np.array([[float(r[k]) for r in rows] for k in cols])
    if not np.isfinite(value).all():
        raise ValueError('Nonfinite evidence')
    return value


def project(weights, values):
    return -math.fsum((weights * values).ravel())


def allocate(game, names):
    if set(game) != set(range(2**len(names))):
        raise ValueError('Incomplete factorial matrix')
    return factors.shapley_values(game, names)


def partitions(device, geo):
    result = {}
    for depth in (.05, .1, .2):
        masks, _, _ = matrix.spatial.old.m78.supports(device, geo, depth)
        group = {k: masks[k] for k in ('source', 'channel', 'drain', 'substrate')}
        group['nonSi'] = ~masks['all_si']
        result[f'depth_{depth}um'] = group
    triple = masks['interface'] & masks['si_nitride_interface']
    result['interface_tags'] = {'SiSiO2': masks['interface'] & ~triple,
        'SiNitride': masks['si_nitride_interface'] & ~triple, 'triple': triple,
        'otherSi': masks['all_si'] & ~(masks['interface'] | masks['si_nitride_interface']),
        'nonSi': ~masks['all_si']}
    for group in result.values():
        if not np.all(sum(v.astype(int) for v in group.values()) == 1):
            raise ValueError('Overlapping or incomplete spatial partition')
    return result


def prepare():
    a.verify(matrix.FREEZE)
    a.verify(matrix.OUT / 'evidence.json')
    fixed.verify()
    cases, files = [], [Path(__file__).resolve(), matrix.chain.RUNNER,
        matrix.OUT / 'evidence.json', fixed.EVIDENCE,
        REPO / 'src/equation/CoupledDDAssembler.cpp', REPO / 'src/solver/NewtonSolver.cpp',
        REPO / 'tests/regression/test_simplemos_calibrated_transport.py']
    for w in fixed.workflows():
        tag, _, export, mapped = fixed.paths(w)
        c = next(copy.deepcopy(r) for r in a.read(matrix.CONTRACT)['cases'] if r['device'] == w['device'] and float(r['vd']) == float(w['drain_voltage_V']))
        key = c['case']
        root = LOCAL / key
        base = matrix.spatial.precision.LOCAL / key
        adj = (matrix.chain.LOCAL / 'adjoint/adjoint.csv' if key == matrix.REUSED else matrix.LOCAL / key / 'adjoint/adjoint.csv')
        c.update(mapped=str(mapped), export=str(export), adjoint=str(adj))
        c['nodes'] = 1480 if c['device'] == 'n19' else 1482
        c['factor'] = a.read(matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
        deck = a.read(base / 'config.json')
        deck.pop('output_state_file')
        deck['state_file'] = str(base / 'state.csv')
        files += [base / 'config.json', base / 'state.csv', adj, mapped]
        files += [Path(deck[k]) for k in ('mesh_file', 'node_doping_file', 'materials_file')]
        files += list(export.rglob('*.csv'))
        for role, state in (('strict', base / 'state.csv'), ('mapped', mapped)):
            for name, typ, extra in (
                ('functional', 'terminal_current_functional_probe', {'residual_output_csv': str(root / role / 'residual.csv')}),
                ('carrier', 'newton_carrier_term_probe', {'output_csv': str(root / role / 'carrier.csv'), 'carrier_term_probe': {'solved_equation_terms': True}}),
                ('edges', 'sg_edge_flux_probe', {'output_csv': str(root / role / 'edges.csv')})):
                d = dict(deck, state_file=str(state), simulation_type=typ, contact='drain', **extra)
                path = root / role / (name + '.json')
                a.write(path, d)
                files.append(path)
        for name, typ in (('secant', 'electron_transport_secant_factor_probe'), ('four_factor', 'electron_transport_factor_probe')):
            d = dict(deck, simulation_type=typ, replacement_state_file=str(mapped), output_csv=str(root / (name + '.csv')))
            if name == 'secant':
                d['edge_output_csv'] = str(root / 'secant_edges.csv')
            path = root / (name + '.json')
            a.write(path, d)
            files.append(path)
        # Test the corrected Jacobian along a small tangent of the actual state gap.
        old, new = ordered(base / 'state.csv', c['nodes']), ordered(mapped, c['nodes'])
        paths = {}
        for component, prefix in (('psi', None), ('phin', 'electron'), ('phip', 'hole')):
            rows = []
            for u, v in zip(old, new, strict=True):
                uvalue = Decimal(u[component]) if prefix is None else Decimal(u[prefix + '_qf_reference_V']) + Decimal(u[prefix + '_qf_increment_V'])
                value = (Decimal(v[component]) - uvalue) * Decimal(str(SCALE))
                rows.append({'node_id': int(u['node_id']), 'component0': float(value)})
            path = root / 'direction' / (component + '.csv')
            a.write_csv(path, rows)
            paths[component] = str(path)
            files.append(path)
        d = dict(deck, simulation_type='actual_direction_jvp', direction_files=paths,
            plus_state_file=str(mapped), minus_state_file=str(base / 'state.csv'), output_csv=str(root / 'direction/jvp.csv'))
        path = root / 'direction/config.json'
        a.write(path, d)
        files.append(path)
        cases.append(c)
    files += [fixed.upstream.OUT / 'm80_case_ledger.csv', fixed.old.m74.CASES]
    a.write(CONTRACT, {'status': 'frozen_before_execution', 'cases': cases,
        'weights': 'Corrected full DD drain adjoint at strict Vela state; never old mapped-state adjoint.',
        'state_mapping': 'Existing M80b coherent mapping: native potentials, n/p reconstructed with unchanged Vela ni/Vt.',
        'definitions': {'deltaF': 'F(mapped)-F(strict)', 'state_prediction': '-lambda_strict dot deltaF',
            'exact_gap': 'I(strict)-I(native) = state_prediction + nonlinear_remainder + mapped_target_difference',
            'factors': 'Exact Shapley allocation of Vela mixed-state electron transport, not native operator error.',
            'geometry': 'Frozen-state Poisson source screens only: region-local dielectric, electron barycentric-Si/all-cell and signed-Si/all-cell volumes.',
            'proxy': 'Native nodal electron mobility arithmetic edge average at mapped state; include direct current term and residual feedback, not native edge mobility.'},
        'gates': {'reconstruction_relative': 1e-5, 'endpoint_relative': 1e-8,
            'partition_relative': 1e-10, 'factor_closure_relative': 1e-8,
            'jvp_relative': .001, 'duality_relative': 1e-8, 'linear_remainder_fraction': .2},
        'direction_scale': SCALE, 'new_read_only_calls': 36, 'new_nonlinear_solves': 0,
        'new_sentaurus_runs': 0, 'production_changes': False, 'm82_released': False, 'm83_released': False})
    files.append(CONTRACT)
    matrix.freeze(FREEZE, files)
    print('Frozen four strict states, 36 read-only probes and diagnostic gates', flush=True)


def run():
    a.verify(FREEZE)
    for c in a.read(CONTRACT)['cases']:
        root = LOCAL / c['case']
        for role in ('strict', 'mapped'):
            for name in ('functional', 'carrier', 'edges'):
                s = matrix.execute(root / role / (name + '.json'), True)
                assert s['exit_code'] == 0
        for path in (root / 'secant.json', root / 'four_factor.json', root / 'direction/config.json'):
            s = matrix.execute(path, True)
            assert s['exit_code'] == 0
        print(c['case'], 'nine read-only probes complete', flush=True)


def analyze():
    a.verify(FREEZE)
    cfg = a.read(CONTRACT)
    ledgers = {k: [] for k in ('case', 'component', 'spatial', 'factor', 'factor_spatial', 'geometry', 'mobility_proxy', 'jvp')}
    references = {(r['device'], float(r['drain_voltage_V'])): r for r in a.rows(fixed.old.m74.CASES)}
    for c in cfg['cases']:
        key, count = c['case'], c['nodes']
        root = LOCAL / key
        geo = matrix.spatial.m73.Geometry(c['device'])
        groups = partitions(c['device'], geo)
        adjrows = ordered(Path(c['adjoint']), count)
        weights = array(adjrows, ('lambda_poisson', 'lambda_electron', 'lambda_hole'))
        residual = {r: array(ordered(root / r / 'residual.csv', count), ('psi_residual', 'phin_residual', 'phip_residual')) for r in ('strict', 'mapped')}
        terms = {r: ordered(root / r / 'carrier.csv', count) for r in ('strict', 'mapped')}
        status = {r: a.read(root / r / 'functional.status.json') for r in ('strict', 'mapped')}
        states = {r: ordered(Path(a.read(root / r / 'functional.json')['state_file']), count) for r in ('strict', 'mapped')}
        delta = residual['mapped'] - residual['strict']
        bv, sv = (array(states[r], ('psi', 'electrons_m3', 'holes_m3')) for r in ('strict', 'mapped'))
        comp = {}
        pois = {'dielectric': c['factor'] * (geo.matrices['legacy'] @ (sv[0]-bv[0])),
            'electron': c['factor'] * fixed.Q * (sv[1]-bv[1]) * geo.volumes['all_cell'],
            'hole': -c['factor'] * fixed.Q * (sv[2]-bv[2]) * geo.volumes['all_cell'], 'dopant': np.zeros(count)}
        for name, value in pois.items():
            v = np.zeros_like(delta)
            v[0, geo.free] = value[geo.free]
            comp['poisson_' + name] = v
        v = np.zeros_like(delta)
        v[0, geo.contact_nodes] = delta[0, geo.contact_nodes]
        comp['poisson_boundary'] = v
        v = np.zeros_like(delta)
        v[0] = delta[0] - sum(x[0] for x in comp.values())
        comp['poisson_reconstruction_remainder'] = v
        for name, (column, _) in COLS.items():
            v = np.zeros_like(delta)
            v[1 if name.startswith('electron') else 2] = array(terms['mapped'], (column,))[0] - array(terms['strict'], (column,))[0]
            comp[name] = v
        comp['continuity_reconstruction_remainder'] = delta - sum(comp.values())
        denominator = max(float(np.linalg.norm(delta)), 1e-300)
        reconstruction = float(np.linalg.norm(comp['poisson_reconstruction_remainder'] + comp['continuity_reconstruction_remainder'])) / denominator
        ref = references[c['device'], float(c['vd'])]
        iv, im, isn = status['strict']['current_A_per_um'], status['mapped']['current_A_per_um'], float(ref['sentaurus_current_A_per_um'])
        assert abs(iv/float(c['current_A_per_um'])-1)<=1e-8
        mesh = a.read(Path(a.read(root/'strict/edges.json')['mesh_file']))
        edge_identity = []
        for role in ('strict','mapped'):
            erows = a.rows(root/role/'edges.csv')
            currents = fixed.sum_contacts(erows,mesh,
                np.array([float(r['electron_particle_line_flux_per_m_s']) for r in erows],dtype=np.longdouble),
                np.array([float(r['hole_particle_line_flux_per_m_s']) for r in erows],dtype=np.longdouble))
            err=abs(currents['drain']/status[role]['current_A_per_um']-1)
            assert err<=1e-8
            edge_identity.append(err)
        predicted = project(weights, delta)
        actual = iv-isn
        remainder = iv-im-predicted
        checks = []
        for name, value in comp.items():
            influence = -np.sum(weights * value, axis=0)
            contribution = math.fsum(influence)
            ledgers['component'].append({'case': key, 'component': name, 'current_A_per_um': contribution, 'relative_to_native_Id': contribution/isn, 'absolute_node_sum_A_per_um': float(np.sum(abs(influence)))})
            for partition, masks in groups.items():
                region_sum = 0.
                for region, mask in masks.items():
                    val = math.fsum(influence[mask])
                    region_sum += val
                    ledgers['spatial'].append({'case': key, 'component': name, 'partition': partition, 'region': region, 'current_A_per_um': val})
                assert abs(region_sum-contribution) <= 1e-10 * max(float(np.sum(abs(influence))), 1e-300)
        # Exact factor allocation using both stable three-factor and four-factor diagnostics.
        for kind, names in (('secant', ('mobility', 'sg_secant_conductance', 'qf_log_imbalance')), ('four_factor', factors.FACTORS)):
            rows = a.rows(root / (kind + '.csv'))
            assert len(rows) == (2**len(names)) * count
            games = {i: {} for i in range(count)}
            for r in rows:
                node, mask = int(r['node_id']), int(r['mask'])
                assert mask not in games[node]
                games[node][mask] = float(r['electron_flux'])
            influence = {name: np.zeros(count) for name in names}
            endpoint = np.array([games[i][2**len(names)-1]-games[i][0] for i in range(count)])
            target = comp['electron_transport'][1]
            endpoint_err = float(np.linalg.norm(endpoint-target)) / max(float(np.linalg.norm(target)), 1e-300)
            assert endpoint_err <= cfg['gates']['endpoint_relative']
            for node, game in games.items():
                parts = allocate(game, names)
                for name, val in parts.items():
                    influence[name][node] = -weights[1, node] * val
            total = math.fsum(v for values in influence.values() for v in values)
            expected = project(weights, comp['electron_transport'])
            closure = abs(total-expected) / max(abs(expected), 1e-300)
            assert closure <= cfg['gates']['factor_closure_relative']
            checks.append(closure)
            for name, values in influence.items():
                ledgers['factor'].append({'case': key, 'factorization': kind, 'factor': name, 'current_A_per_um': math.fsum(values), 'relative_to_native_Id': math.fsum(values)/isn})
                for partition, masks in groups.items():
                    for region, mask in masks.items():
                        ledgers['factor_spatial'].append({'case': key, 'factorization': kind, 'factor': name, 'partition': partition, 'region': region, 'current_A_per_um': math.fsum(values[mask])})
        # Frozen-state geometric operator differences; no causal release from these screens.
        for name, source in (
            ('region_local_dielectric', c['factor'] * ((geo.matrices['region_local']-geo.matrices['legacy']) @ bv[0])),
            ('electron_barycentric_si', c['factor'] * fixed.Q * bv[1] * (geo.volumes['barycentric_si']-geo.volumes['all_cell'])),
            ('electron_signed_si', c['factor'] * fixed.Q * bv[1] * (geo.volumes['signed_si']-geo.volumes['all_cell']))):
            source[geo.contact_nodes] = 0
            values = -weights[0] * source
            hist = float(ref['candidate_vela_current_A_per_um'])-float(ref['baseline_vela_current_A_per_um']) if name == 'electron_barycentric_si' else None
            ledgers['geometry'].append({'case': key, 'candidate': name, 'predicted_delta_Id_A_per_um': math.fsum(values), 'relative_to_native_Id': math.fsum(values)/isn, 'historical_M74_delta_Id_A_per_um': hist, 'self_consistent_qualified_this_round': False})
            for region, mask in groups['depth_0.05um'].items():
                ledgers['spatial'].append({'case': key, 'component': 'screen_' + name, 'partition': 'depth_0.05um', 'region': region, 'current_A_per_um': math.fsum(values[mask])})
        # Reconstructed native nodal mobility is a declared proxy at the common mapped state.
        edges = a.rows(root / 'mapped/edges.csv')
        mesh = a.read(Path(a.read(root / 'mapped/edges.json')['mesh_file']))
        mu = fixed.scalar(Path(c['export']), 'eMobility')
        dn, physical = [], []
        nodes = np.zeros(count)
        for e in edges:
            i, j = int(e['node0']), int(e['node1'])
            ratio = (.5e-4 * (mu[i]+mu[j]) / float(e['electron_mobility_m2_V_s']) if i in mu and j in mu and float(e['electron_mobility_m2_V_s']) > 0 else 1.)
            change = float(e['electron_flux']) * (ratio-1)
            nodes[i] += change
            nodes[j] -= change
            dn.append(change)
            physical.append(float(e['electron_particle_line_flux_per_m_s']) * (ratio-1))
        constrained = np.array([bool(r['contacts']) for r in adjrows])
        nodes[constrained] = 0
        direct = fixed.sum_contacts(edges, mesh, np.array(physical, dtype=np.longdouble), np.zeros(len(edges)))['drain']
        feedback = -float(weights[1] @ nodes)
        ledgers['mobility_proxy'].append({'case': key, 'direct_A_per_um': direct, 'feedback_A_per_um': feedback, 'total_linear_A_per_um': direct+feedback, 'relative_to_native_Id': (direct+feedback)/isn, 'native_edge_operator_recovered': False})
        # Qualify corrected Jv in a small tangent of the actual mapped-state gap.
        jrows = a.rows(root / 'direction/jvp.csv')
        assert len(jrows) == 36*count
        gradient = array(adjrows, ('dI_dpsi_scaled', 'dI_dphin_scaled', 'dI_dphip_scaled'))
        block = {'psi': 0, 'phin': 1, 'phip': 2}
        for h in (1., .5, .25):
            group = [r for r in jrows if r['mode'] == 'all' and float(r['step']) == h]
            assert len(group) == 3*count
            analytic = math.fsum(weights[block[r['row_block']], int(r['node_id'])] * float(r['analytic']) for r in group)
            fd = math.fsum(weights[block[r['row_block']], int(r['node_id'])] * float(r['fd']) for r in group)
            gdot = math.fsum(gradient[block[r['row_block']], int(r['node_id'])] * float(r['direction_scaled']) for r in group)
            fd_error, duality = abs(analytic-fd)/max(abs(gdot),1e-300), abs(analytic-gdot)/max(abs(gdot),1e-300)
            ledgers['jvp'].append({'case': key, 'step': h, 'gradient_dot_direction_A_per_um': gdot, 'lambda_Jv_A_per_um': analytic, 'lambda_FD_A_per_um': fd, 'relative_error': fd_error, 'duality_relative': duality, 'passed': fd_error<=.001 and duality<=1e-8})
        row = {'case': key, 'device': c['device'], 'vd': float(c['vd']), 'vg': float(c['vg']),
            'strict_Id_A_per_um': iv, 'mapped_Id_A_per_um': im, 'native_Id_A_per_um': isn,
            'actual_gap_A_per_um': actual, 'state_gap_A_per_um': iv-im, 'target_gap_A_per_um': im-isn,
            'predicted_state_gap_A_per_um': predicted, 'nonlinear_remainder_A_per_um': remainder,
            'nonlinear_fraction_of_actual_gap': abs(remainder)/abs(actual), 'reconstruction_relative': reconstruction,
            'max_edge_terminal_relative': max(edge_identity),
            'max_factor_closure_relative': max(checks), 'linear_screen_passed': abs(remainder)<=max(.2*abs(actual),1e-4*abs(isn)*math.log(10))}
        ledgers['case'].append(row)
        print(key, 'remainder fraction', row['nonlinear_fraction_of_actual_gap'], flush=True)
    for name, rows in ledgers.items():
        a.write_csv(OUT / (name + '.csv'), rows)
    # Exact normalized error budget and high-minus-low pairing, including target and remainder.
    pairs = []
    for vd in (.05, 1.):
        lo, hi = [next(r for r in ledgers['case'] if r['vd']==vd and r['device']==d) for d in ('n19','n23')]
        for name in list(COLS) + ['poisson_dielectric','poisson_electron','poisson_hole','poisson_dopant','poisson_boundary','poisson_reconstruction_remainder','continuity_reconstruction_remainder','nonlinear_remainder','target_gap']:
            vals = []
            for c in (lo,hi):
                val = c[name+'_A_per_um'] if name in ('nonlinear_remainder','target_gap') else next(r['current_A_per_um'] for r in ledgers['component'] if r['case']==c['case'] and r['component']==name)
                vals.append(val/c['native_Id_A_per_um'])
            pairs.append({'vd':vd,'component':name,'low_relative_error_contribution':vals[0],'high_relative_error_contribution':vals[1],'high_minus_low_relative_error':vals[1]-vals[0]})
    a.write_csv(OUT / 'pair.csv', pairs)
    a.write(OUT / 'result.json', {'cases':4,'new_read_only_calls':36,'new_nonlinear_solves':0,'new_sentaurus_runs':0,
        'representation_passed':all(r['reconstruction_relative']<=1e-5 for r in ledgers['case']),
        'jvp_passed':all(r['passed'] for r in ledgers['jvp']),
        'linear_screen_passed_count':sum(r['linear_screen_passed'] for r in ledgers['case']),
        'm82_released':False,'m83_released':False})


def seal():
    a.verify(FREEZE)
    matrix.freeze(OUT/'evidence.json',[DOC]+[f for root in (LOCAL,OUT) for f in root.rglob('*') if f.is_file()])


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','run','analyze','seal','verify'))
    action=parser.parse_args().action
    if action=='verify':
        a.verify(FREEZE)
        a.verify(OUT/'evidence.json')
        print('Calibrated transport evidence verified')
    else:
        globals()[action]()
