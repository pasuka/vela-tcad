"""Default-off, bounded neighboring-row QF repair; never a device acceptance gate.

Uses existing frozen probes and the unchanged strict Newton validator. No ledger
is written. FD, patch limits and damping are proposal policies, not tolerances.
"""
from dataclasses import asdict, dataclass
from pathlib import Path
import argparse
import json
import math
import time

import numpy as np


@dataclass(frozen=True)
class Policy:
    enabled: bool = False
    max_cycles: int = 3
    max_dofs: int = 8
    neighbor_fraction: float = .5
    fd_step_V: float = 1e-6
    max_qf_step_V: float = 1e-4
    alphas: tuple = (1., .5, .25, .125)
    decrease: float = 1e-4
    predict_boundary: bool = False
    predict_boundary_cycles: int = 1

    def validate(self):
        # Explicit larger diagnostic budgets are allowed; defaults remain 3/8.
        # These bound repair work, never the row or global acceptance criteria.
        if (type(self.enabled) is not bool or type(self.predict_boundary) is not bool
                or type(self.max_cycles) is not int
                or not 1 <= self.max_cycles <= 6 or type(self.max_dofs) is not int
                or not 1 <= self.max_dofs <= 64 or type(self.predict_boundary_cycles) is not int
                or not 1 <= self.predict_boundary_cycles <= 6):
            raise ValueError('invalid bounded repair policy')
        values = (self.neighbor_fraction, self.fd_step_V, self.max_qf_step_V, self.decrease)
        if (not all(math.isfinite(v) and v > 0 for v in values)
                or self.neighbor_fraction > 1 or self.decrease >= 1
                or self.max_qf_step_V > 1e-3 or self.fd_step_V > 1e-4
                or not 1 <= len(self.alphas) <= 4
                or any(not math.isfinite(a) or not 0 < a <= 1 for a in self.alphas)
                or any(a <= b for a, b in zip(self.alphas, self.alphas[1:]))):
            raise ValueError('invalid finite-difference/damping policy')


def row_table(path, floor):
    """All active rows, deliberately at least as strict as qualified-row gating."""
    import csv
    result = {}
    seen = set()
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        for row in csv.DictReader(f):
            node = int(row['node_id'])
            if node in seen:
                raise ValueError('duplicate node in row probe')
            seen.add(node)
            for carrier in ('electron', 'hole'):
                active = row[carrier + '_continuity_active']
                if active not in ('0', '1'):
                    raise ValueError('invalid active flag')
                if active == '0':
                    continue
                r, flux, rec, imp = [float(row[carrier + '_' + k]) for k in
                                     ('residual', 'flux_abs_sum', 'recombination', 'impact')]
                if not all(math.isfinite(v) for v in (r, flux, rec, imp)):
                    raise ValueError('nonfinite row probe')
                scale = max(abs(flux), abs(rec), abs(imp), floor)
                result[(node, carrier)] = dict(residual=r, scale=scale, ratio=abs(r)/scale)
    if not result:
        raise ValueError('empty active row set')
    return result


def select_patch(rows, neighbors, eps, policy):
    bad = {key for key, r in rows.items() if r['ratio'] > eps}
    patch = set(bad)
    frontier = sorted(bad)
    # Include connected near-threshold rows, including neighbors of a selected
    # weak row. A one-ring-only patch can repeatedly export error one edge away.
    while frontier:
        node, carrier = frontier.pop()
        for neighbor in neighbors.get(node, ()):
            key = (neighbor, carrier)
            if (key not in patch and key in rows
                    and rows[key]['ratio'] >= policy.neighbor_fraction * eps):
                patch.add(key)
                frontier.append(key)
    return sorted(patch)


def boundary_candidates(rows, patch, columns, delta, neighbors, eps, policy):
    """Rank same-carrier one-ring rows by predicted full-step residual.

    Uses the current scale for prediction only. Actual trial gates recompute it.
    No node IDs or previous case results enter the decision.
    """
    boundary = {(n, c) for i, c in patch for n in neighbors.get(i, ())}
    ranked = []
    for key in boundary - set(patch):
        if key not in rows:
            continue
        value = rows[key]['residual'] + math.fsum(
            columns[k][key] * float(d) for k, d in zip(patch, delta))
        ratio = abs(value) / rows[key]['scale']
        if not math.isfinite(ratio):
            raise ValueError('nonfinite boundary prediction')
        if ratio >= policy.neighbor_fraction * eps:
            ranked.append(dict(key=key, predicted_ratio=ratio))
    return sorted(ranked, key=lambda r: (-r['predicted_ratio'], r['key']))


def repair(seed, *, rows_at, global_at, trial_at, neighbors, eps, policy):
    """Callbacks are bounded by the caller's deadline; output is a seed only.

    FD probes need not be qualified. Every selected trial must pass the original
    global-condition screen, improve the maximum over *all* active rows, and
    preserve the exact active set. Prefer an all-row-passing damped trial over a
    merely decreasing full step. Final strict validation is the caller's job.
    """
    policy.validate()
    if not math.isfinite(eps) or eps <= 0:
        raise ValueError('invalid row threshold')
    result = dict(role='initial_guess_only', device_qualified=False, policy=asdict(policy),
                  steps=[], state_file='', rows_passed=False)
    if not policy.enabled:
        return dict(result, stop_reason='disabled')
    current = seed
    rows = rows_at(current, 'initial_rows')
    for cycle in range(policy.max_cycles):
        g = global_at(current, f'c{cycle}_current_global')
        if not g['passed']:
            return dict(result, stop_reason='original_global_screen_failed', global_screen=g)
        merit = max(r['ratio'] for r in rows.values())
        if merit <= eps:
            return dict(result, rows_passed=True, state_file=str(current), stop_reason='rows_satisfied')
        patch = select_patch(rows, neighbors, eps, policy)
        if len(patch) > policy.max_dofs:
            return dict(result, stop_reason='patch_budget_exceeded', patch_size=len(patch))
        step = dict(cycle=cycle, patch=patch, before_max_ratio=merit, trials=[])
        result['steps'].append(step)
        columns = {}
        predict = policy.predict_boundary and cycle < policy.predict_boundary_cycles
        # At most one predictive expansion per cycle. Cache columns so expanding
        # the patch never repeats a finite difference at the same state.
        for expansion in range(2 if predict else 1):
            for key in patch:
                if key in columns:
                    continue
                j = len(columns); probes = []
                for sign in (-1, 1):
                    p = trial_at(current, [key], [sign * policy.fd_step_V], f'c{cycle}_fd{j}_{sign}')
                    rr = rows_at(p, f'c{cycle}_fd{j}_{sign}_rows')
                    if rr.keys() != rows.keys():
                        raise ValueError('active rows changed during FD probe')
                    probes.append(rr)
                columns[key] = {k: (probes[1][k]['residual']-probes[0][k]['residual']) /
                                (2*policy.fd_step_V) for k in rows}
            jac = np.array([[columns[j][i] for j in patch] for i in patch])
            scales = np.array([rows[k]['scale'] for k in patch])
            rhs = np.array([rows[k]['residual'] for k in patch])
            try:
                delta = np.linalg.solve(jac / scales[:, None], -rhs / scales)
            except np.linalg.LinAlgError:
                return dict(result, stop_reason='singular_local_jacobian')
            if not np.all(np.isfinite(delta)) or np.max(np.abs(delta)) > policy.max_qf_step_V:
                return dict(result, stop_reason='local_step_guard')
            if not predict or expansion:
                break
            ranked = boundary_candidates(rows, patch, columns, delta, neighbors, eps, policy)
            added = [r['key'] for r in ranked[:policy.max_dofs-len(patch)]]
            step['boundary_prediction'] = dict(initial_patch=patch, ranked=ranked, added=added,
                                               cap_reached=len(ranked)>len(added))
            if not added:
                break
            patch = sorted(set(patch) | set(added))
            step['patch'] = patch
        step['fd_columns'] = list(columns)
        step.update(delta_qf_V=delta.tolist(), jacobian=jac.tolist(),
                    scaled_linear_residual=float(np.max(np.abs((jac @ delta + rhs)/scales))))
        best = None
        candidates = []
        for alpha in policy.alphas:
            p = trial_at(current, patch, alpha*delta, f'c{cycle}_a{alpha}')
            rr = rows_at(p, f'c{cycle}_a{alpha}_rows')
            if rr.keys() != rows.keys():
                raise ValueError('active rows changed in damped trial')
            value = max(r['ratio'] for r in rr.values())
            descent = value < merit and value <= merit*(1-policy.decrease*alpha)
            entry = dict(alpha=alpha, max_ratio=value, global_screen=None,
                         descent=bool(descent), state_file=str(p))
            step['trials'].append(entry)
            if not descent:
                continue
            candidates.append((value, p, rr, alpha, entry))
            if value <= eps:
                entry['global_screen'] = global_at(p, f'c{cycle}_a{alpha}_global')
                if entry['global_screen']['passed']:
                    best = (value, p, rr, alpha)
                    break
        # Cost control: test the global screen only after row-based ranking.
        # Untested trials retain a null screen and can never be selected.
        if best is None:
            for value, p, rr, alpha, entry in sorted(candidates, key=lambda t:t[0]):
                if entry['global_screen'] is None:
                    entry['global_screen'] = global_at(p, f'c{cycle}_a{alpha}_global')
                if entry['global_screen']['passed']:
                    best = (value, p, rr, alpha)
                    break
        if best is None:
            return dict(result, stop_reason='no_globally_feasible_row_descent')
        _, current, rows, alpha = best
        step['selected_alpha'] = alpha
        step['after_max_ratio'] = best[0]
        if best[0] <= eps:
            return dict(result, rows_passed=True, state_file=str(current), stop_reason='rows_satisfied')
    return dict(result, stop_reason='bounded_cycles_exhausted')


def run(base, seed, bias, runner, repo, output, policy=Policy(), max_seconds=600, child_seconds=120):
    """Reusable explicit opt-in adapter. Successful result includes strict proof.

    Supports Boltzmann QF increments and explicit block residual scales only.
    Frozen screens cannot grant convergence. The unchanged full Newton solve,
    all-active row audit and bound file validation remain mandatory.
    """
    import run_slot_ldmos_loadline_continuation as ll
    import run_slot_ldmos_noavalanche_continuation as helper
    import slot_candidate_conservation as conservation
    import slot_state_archive as states
    policy.validate()
    if not policy.enabled:
        return dict(passed=False, stop_reason='disabled', state_file='')
    if not all(math.isfinite(x) and x > 0 for x in (max_seconds, child_seconds)):
        raise ValueError('invalid time budget')
    output = Path(output).resolve(); output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); deadline = started + max_seconds
    base, identity = ll.freeze_inputs(base, repo, runner)
    solver = base['solver']; gate = solver['carrier_row_convergence']
    statistics = solver.get('carrier_statistics', 'boltzmann')
    statistics = statistics.get('model', 'boltzmann') if isinstance(statistics, dict) else statistics
    statistics = solver.get('statistics', statistics)
    if (solver.get('method') != 'newton' or gate.get('mode') != 'enforce'
            or solver['global_continuity_closure'].get('mode') != 'enforce'
            or solver.get('residual_norm', 'block') != 'block'
            or statistics != 'boltzmann'):
        raise ValueError('unsupported strict/Boltzmann/block contract')
    scales = solver.get('residual_scales', {})
    if any(not math.isfinite(scales.get(k, 0)) or scales.get(k, 0) <= 0 for k in ('psi','phin','phip')):
        raise ValueError('explicit positive residual scales required')
    weights = solver.get('residual_weights', dict(psi=1., phin=1., phip=1.))
    if any(weights.get(k, 1.) != 1. for k in ('psi','phin','phip')):
        raise ValueError('non-default residual weights unsupported')
    if (not math.isfinite(bias) or not math.isfinite(gate['scale_floor'])
            or gate['scale_floor'] <= 0 or not math.isfinite(solver['abstol'])
            or solver['abstol'] <= 0):
        raise ValueError('invalid bias/residual scale')
    seed = Path(seed).resolve(); states.read(seed, base, repo)
    protected = {**identity['files'], **ll.hashes([seed, runner, Path(__file__),
        Path(ll.__file__), Path(helper.__file__), Path(states.__file__), Path(conservation.__file__)])}
    ll.atomic_json(output/'identity.json', dict(inputs=protected, policy=asdict(policy), base=base))
    mesh = json.loads(Path(base['mesh_file']).read_text(encoding='utf-8-sig'))
    neighbors = {}
    for tri in mesh['triangles']:
        ids = tri['node_ids'] if isinstance(tri, dict) else tri
        for a in ids:
            neighbors.setdefault(a, set()).update(b for b in ids if b != a)
    runs = []
    def budget():
        left = deadline-time.monotonic()
        if left <= 0:
            raise ll.BudgetExhausted('local row repair wall clock')
        return min(child_seconds, left)
    def execute(doc, folder):
        folder.mkdir(parents=True, exist_ok=False)
        helper.redirect_diagnostics(doc, folder)
        cfg = folder/'config.json'; ll.atomic_json(cfg, doc)
        code, timed, elapsed = helper.execute_config(runner, repo, cfg, folder, budget())
        runs.append(dict(config=str(cfg), exit_code=code, timed_out=timed, elapsed_s=elapsed))
        ll.atomic_json(output/'runs.json', runs)
        if timed:
            raise ll.BudgetExhausted('local row repair child timeout')
        if code != 0:
            raise ll.NumericalFailure('local row repair probe/validation failed')
        return cfg
    def probe_doc(state, kind, target):
        doc = ll.clone(base); driven = doc.pop('sweep')['contact']
        for contact in doc['contacts']:
            if contact['name'] == driven:
                contact['bias'] = bias
        doc.update(simulation_type=kind, state_file=str(state), state_format='hdf5', output_csv=str(target))
        return doc
    def rows_at(state, name):
        folder = output/name
        doc = probe_doc(state, 'newton_carrier_term_probe', folder/'terms.csv')
        doc['carrier_term_probe'] = dict(solved_equation_terms=True)
        execute(doc, folder)
        return row_table(folder/'terms.csv', gate['scale_floor'])
    def global_at(state, name):
        folder = output/name
        doc = probe_doc(state, 'newton_candidate_conservation_probe', folder/'result.json')
        doc['candidate_port_diagnostics'] = True
        execute(doc, folder)
        value = json.loads((folder/'result.json').read_text())
        g = solver['global_continuity_closure']
        closure = conservation.validate_result(value, dict(tolerance=g['tolerance'],
                    source_floor=g['source_floor'], absolute_tolerance=g['absolute_tolerance']))
        blocks = value['frozen_block_residuals']
        norm = math.sqrt(sum((blocks[k]/scales[k])**2 for k in ('psi','phin','phip')))
        absolute = solver.get('block_absolute_convergence', {})
        if absolute.get('mode', 'off') != 'off':
            raise ValueError('authoritative absolute-block contract requires dedicated adapter')
        return dict(passed=bool(closure and math.isfinite(norm) and norm <= solver['abstol']),
                    normalized_residual=norm, abstol=solver['abstol'], closure_satisfied=closure,
                    blocks=blocks)
    def trial_at(state, patch, delta, name):
        fields, meta = states.read(state, base, repo)
        temp = float(solver.get('temperature_K', 300.))
        if not math.isfinite(temp) or temp <= 0:
            raise ValueError('invalid temperature')
        if 'temperature_K' in fields:
            raise ValueError('spatial temperature field unsupported')
        vt = 1.380649e-23*temp/1.602176634e-19
        for (i, carrier), dv in zip(patch, delta):
            prefix = carrier + '_qf_'; physical = 'phin' if carrier == 'electron' else 'phip'
            inc = prefix+'increment_V'; ref = prefix+'reference_V'
            fields[inc][i] = math.fsum([float(fields[inc][i]), float(dv)])
            fields[physical][i] = math.fsum([float(fields[ref][i]), float(fields[inc][i])])
            density = 'electrons_m3' if carrier == 'electron' else 'holes_m3'
            fields[density][i] *= math.exp((-1 if carrier == 'electron' else 1)*float(dv)/vt)
        path = output/(name+'.h5')
        states.write_candidate(path, fields, meta, base, [state], bias, 'local_row_repair_initial_guess')
        return path
    result = dict(passed=False, completed=False, state_file='', runs=runs)
    try:
        proposed = repair(seed, rows_at=rows_at, global_at=global_at, trial_at=trial_at,
                          neighbors=neighbors, eps=gate['eps_row'], policy=policy)
        result['repair'] = proposed
        ll.atomic_json(output/'repair.json', proposed)
        if proposed['rows_passed']:
            folder = output/'strict'; doc = ll.clone(base)
            doc['sweep'].update(initial_state_file=proposed['state_file'], bias_points=[bias],
                write_state_file=str(folder/'final_state.h5'), write_state_every_point_prefix=str(folder/'states/s'))
            doc['output_csv'] = str(folder/'iv.csv')
            cfg = execute(doc, folder)
            backend = ll.SolverBackend(base, repo, runner, timeout_seconds=budget(),
                max_newton_iterations=solver['max_iter'], carrier_corrector_passes=0, restore_poly_equilibrium=False)
            backend.set_deadline(deadline)
            point = backend._qualified(doc, cfg, folder, folder/'final_state.h5', bias, True,
                                        dict(inner_V=bias, state_file=proposed['state_file']))
            backend.validate(point)
            result.update(passed=True, state_file=point['state_file'], device_proof=point)
        result['completed'] = True
    except (ll.NumericalFailure, ll.BudgetExhausted) as exc:
        result.update(completed=True, stop_reason=str(exc))
    except Exception as exc:
        result.update(stop_reason='evidence_or_environment_error', error_type=type(exc).__name__, error=str(exc))
        raise
    finally:
        ll.verify_hashes(protected)
        result.update(elapsed_s=time.monotonic()-started, protected_unchanged=True)
        ll.atomic_json(output/'summary.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', type=Path, required=True); p.add_argument('--state', type=Path, required=True)
    p.add_argument('--bias', type=float, help='Override inner V; default is exact sole base bias point')
    p.add_argument('--runner', type=Path, required=True)
    p.add_argument('--repo', type=Path, default=Path.cwd()); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--enable-neighbor-repair', action='store_true')
    p.add_argument('--predict-boundary', action='store_true', help='Bounded FD-impact expansion in the first repair cycle')
    p.add_argument('--max-seconds', type=float, default=600)
    a = p.parse_args()
    base = json.loads(a.base.read_text(encoding='utf-8-sig'))
    if a.bias is None and len(base['sweep'].get('bias_points', [])) != 1:
        p.error('base must contain exactly one bias point, or specify --bias')
    bias = a.bias if a.bias is not None else base['sweep']['bias_points'][0]
    r = run(base, a.state, bias, a.runner, a.repo,
            a.output, Policy(enabled=a.enable_neighbor_repair, predict_boundary=a.predict_boundary), a.max_seconds)
    print(json.dumps({k:v for k,v in r.items() if k in ('passed','completed','stop_reason','elapsed_s')}))
    raise SystemExit(0 if r['passed'] else 2)
