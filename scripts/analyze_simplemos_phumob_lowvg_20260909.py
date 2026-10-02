"""Unfitted low-Vg mobility comparison and held-out G-floor diagnostics."""
import argparse
import math
from decimal import Decimal as D, localcontext
from pathlib import Path

import numpy as np
import validate_simplemos_phumob_lowvg_20260909 as stage
import analyze_simplemos_phumob_native_cells_20260909 as cells
import audit_simplemos_phumob_cross_derivatives_20260908 as hp
import localize_simplemos_phumob_screening_floor_20260909 as floor

a, d, OUT, LOCAL = stage.a, stage.d, stage.OUT, stage.LOCAL


def calibrate():
    a.verify(stage.NOUT / 'export_evidence.json')
    prior_floor = floor.OUT / 'floor_inference.csv'
    a.verify(floor.OUT / 'evidence.json')
    saved = {r['carrier']: float(r['inferred_effective_G_floor']) for r in a.rows(prior_floor)}
    contract = dict(native_states=16, new_fit=False, production_changed=False,
                    heldout='Use previously frozen high-Vg diagnostic G floors; fit zero parameters to these low-Vg fields or currents.',
                    gates=dict(cell_relative=1e-7, edge_relative=1e-7, native_terminal_relative=1e-6, masetti_reclosure_relative=1e-8),
                    precision=100, q_C=cells.Q, thermal_voltage_V=cells.VT,
                    floor_algorithm_identified=False)
    a.write(OUT / 'calibration_contract.json', contract)
    d.matrix.freeze(OUT / 'calibration_freeze.json', [stage.NOUT / 'export_evidence.json', prior_floor,
        floor.OUT / 'evidence.json', Path(__file__).resolve(), Path(cells.__file__), Path(hp.__file__),
        Path(floor.__file__), OUT / 'calibration_contract.json'])
    cells.s.LOCAL, cells.s.OUT = stage.NLOCAL, stage.NOUT
    records, edges, ports, groups = [], [], [], []
    oldpoints = a.rows(stage.native.original.b.OUT / 'native_points.csv')
    with localcontext() as context:
        context.prec = 100
        for job in a.rows(stage.NOUT / 'native_points.csv'):
            assert job['native_qualified'] == 'True'
            data = cells.load(job)
            f = data['fields']
            nodal, trial, clamped = {}, {}, {}
            for node in f['n']:
                state = {k: f[k][node] for k in ('nd', 'na', 'n', 'p')}
                for carno, car in enumerate(('e', 'h')):
                    if job['model'] == 'masetti_control':
                        value = cells.masetti.formula(state['nd'] + state['na'], cells.masetti.PARAMETERS[car])
                        nodal[node, car] = trial[node, car] = value
                        clamped[node, car] = False
                    else:
                        value, _, active = hp.hp([D.from_float(state[k]) for k in ('nd', 'na', 'n', 'p')], carno)
                        nodal[node, car], clamped[node, car] = float(value), bool(active)
                        trial[node, car] = floor.mobility(state, car, saved[car]) if active else float(value)
            predicted, heldout = {}, {}
            for cid, vs in data['cells'].items():
                nodes = [data['mapping'][int(v['vertex'])] for v in vs]
                volumes = [float(v['measure_um2']) for v in vs]
                volume = math.fsum(volumes)
                assert volume > 0
                weights = [v / volume for v in volumes]
                for car in ('e', 'h'):
                    exact = math.fsum(w * nodal[n, car] for n, w in zip(nodes, weights))
                    held = math.fsum(w * trial[n, car] for n, w in zip(nodes, weights))
                    predicted[cid, car], heldout[cid, car] = exact, held
                    ref = data['native'][car][cid]
                    records.append(dict(key=job['key'], model=job['model'], carrier=car, cell=cid,
                        clamped_vertices=sum(clamped[n, car] for n in nodes), native_mu_cm2_V_s=ref,
                        exact_relative=exact / ref - 1, heldout_diagnostic_relative=held / ref - 1))
            current = {arm: {name: [] for name in data['contacts']} for arm in ('native_mu', 'exact_mu', 'heldout_mu')}
            for (i, j), parts in data['parts'].items():
                coefficients = {
                    'native_mu': {car: math.fsum(g * data['native'][car][cid] for cid, g in parts) for car in ('e', 'h')},
                    'exact_mu': {car: math.fsum(g * predicted[cid, car] for cid, g in parts) for car in ('e', 'h')},
                    'heldout_mu': {car: math.fsum(g * heldout[cid, car] for cid, g in parts) for car in ('e', 'h')}}
                for car in ('e', 'h'):
                    ref = coefficients['native_mu'][car]
                    if ref:
                        edges.append(dict(key=job['key'], model=job['model'], carrier=car, node0=i, node1=j,
                            exact_relative=coefficients['exact_mu'][car] / ref - 1,
                            heldout_diagnostic_relative=coefficients['heldout_mu'][car] / ref - 1))
                eta = (f['psi'][j] - f['psi'][i]) / cells.VT
                en = -cells.VT * f['n'][i] * cells.B(-eta - math.log(f['ni'][j] / f['ni'][i])) * math.expm1((f['fn'][i] - f['fn'][j]) / cells.VT)
                ho = -cells.VT * f['p'][i] * cells.B(eta + math.log(f['ni'][i] / f['ni'][j])) * math.expm1((f['fp'][j] - f['fp'][i]) / cells.VT)
                for arm, coef in coefficients.items():
                    for name, ids in data['contacts'].items():
                        current[arm][name].append(-cells.Q * 1e-4 * (int(i in ids) - int(j in ids)) * (en * coef['e'] - ho * coef['h']))
            ref = float(job['Id_A_per_um'])
            currents = {arm: {name: math.fsum(values) for name, values in cc.items()} for arm, cc in current.items()}
            old = next(r for r in oldpoints if r['model'] == 'old_slotboom' and r['case'] == job['case'] and int(r['index']) == int(job['index']))
            port = dict(key=job['key'], model=job['model'], native_Id_A_per_um=ref, runtime_plot_relative=data['plot_error'],
                        masetti_current_reclosure_relative=ref / float(old['Id_A_per_um']) - 1)
            for arm, values in currents.items():
                port[arm + '_Id_A_per_um'] = values['drain']
                port[arm + '_relative'] = values['drain'] / ref - 1
                port[arm + '_kcl_over_Id'] = abs(math.fsum(values.values())) / abs(values['drain'])
            port['fixed_state_exact_minus_native_mu_A_per_um'] = currents['exact_mu']['drain'] - currents['native_mu']['drain']
            port['fixed_state_heldout_minus_exact_mu_A_per_um'] = currents['heldout_mu']['drain'] - currents['exact_mu']['drain']
            port['native_replay_qualified'] = abs(port['native_mu_relative']) <= 1e-6
            ports.append(port)
            print('Native cell/edge replay', job['model'], job['key'], port['native_mu_relative'], flush=True)
    for model, car, count in sorted({(r['model'], r['carrier'], r['clamped_vertices']) for r in records}):
        rr = [r for r in records if (r['model'], r['carrier'], r['clamped_vertices']) == (model, car, count)]
        groups.append(dict(model=model, carrier=car, clamped_vertices=count, cells=len(rr),
            exact_max_relative=max(abs(r['exact_relative']) for r in rr),
            heldout_diagnostic_max_relative=max(abs(r['heldout_diagnostic_relative']) for r in rr),
            exact_failed_cells=sum(abs(r['exact_relative']) > 1e-7 for r in rr),
            heldout_diagnostic_failed_cells=sum(abs(r['heldout_diagnostic_relative']) > 1e-7 for r in rr)))
    for name, rr in [('native_cells.csv', records), ('native_edges.csv', edges), ('native_ports.csv', ports), ('native_groups.csv', groups)]:
        a.write_csv(OUT / name, rr)
    summary = dict(native_states=len(ports), cell_checks=len(records), edge_checks=len(edges),
        native_replay_failures=sum(not r['native_replay_qualified'] for r in ports),
        native_replay_max_relative=max(abs(r['native_mu_relative']) for r in ports),
        masetti_reclosure_failures=sum(abs(r['masetti_current_reclosure_relative']) > 1e-8 for r in ports if r['model'] == 'masetti_control'),
        exact_phumob_cell_max_relative=max(abs(r['exact_relative']) for r in records if r['model'] == 'phumob'),
        heldout_phumob_cell_max_relative=max(abs(r['heldout_diagnostic_relative']) for r in records if r['model'] == 'phumob'),
        native_algorithm_identified=False, diagnostic_fit_parameters_added=0)
    a.write(OUT / 'calibration_summary.json', summary)
    d.matrix.freeze(OUT / 'calibration_evidence.json', [OUT / 'calibration_freeze.json'] + [OUT / n for n in
        ('native_cells.csv', 'native_edges.csv', 'native_ports.csv', 'native_groups.csv', 'calibration_summary.json')])
    print(summary, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=('calibrate',))
    globals()[parser.parse_args().action]()
