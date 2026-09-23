"""Fixed-state diagnostic: existing guarded prediction versus its actual parent.

Uses a frozen package in place; writes only to a new output directory. Timings
include fresh processes and are not cross-point-reuse performance measurements.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(a.package/'scripts'))
    import run_templates_ldmos_linked_d5 as linked
    from analyze_templates_ldmos_d5_newton_cost import trace_cost
    runner = a.package/'build-release/vela_example_runner.exe'
    expected = '75d7a34556f6448c4e869707c8661c9717550c3058f9580387ebcc303a7274bc'
    if linked.digest(runner) != expected:
        raise ValueError('Frozen runner mismatch')
    a.output.mkdir(parents=True, exist_ok=False)
    env = os.environ.copy()
    env.update(VELA_LINEAR_SOLVER='umfpack', VELA_LINEAR_THREADS='1',
               VELA_BLAS_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               OMP_DYNAMIC='FALSE', OMP_MAX_ACTIVE_LEVELS='1', VELA_LINEAR_FACTOR_STATS='0')
    result = dict(runner_sha256=expected, mode='fresh_process_fixed_seed_diagnostic', points=[])
    for gate in (4, 8):
        curve = a.package/f'results/t470p/d5_p31_r0_vg{gate}_U1'
        ledger = linked.read(curve/'fixed/ledger.json')
        predicted = [r for r in ledger['runs'] if r.get('outer_predictor')]
        for target in (.8368462059131809, 1.0368462059131809, 20.2, 33.5):
            r = min(predicted, key=lambda r: abs(r['target_V']-target))
            src = curve/r['case']
            original = linked.read(src/'control.json')
            record = dict(gate_V=gate, target_V=r['target_V'], parent_V=r['parent_V'],
                          cap_V=r['cap_V'], source=str(src), variants={})
            result['points'].append(record)
            for variant in ('guarded', 'no_prediction'):
                dest = a.output/f'vg{gate}_{src.name}_{variant}'
                dest.mkdir()
                def relocate(v):
                    if isinstance(v, dict): return {k: relocate(x) for k, x in v.items()}
                    if isinstance(v, list): return [relocate(x) for x in v]
                    if isinstance(v, str):
                        return v.replace(str(src), str(dest)).replace(src.as_posix(), dest.as_posix())
                    return v
                cfg = relocate(original)
                seed = Path(r['parent_state'] if variant == 'guarded' else r['outer_predictor']['current_state'])
                seed_hash = r['parent_sha256'] if variant == 'guarded' else r['outer_predictor']['current_sha256']
                if linked.digest(seed) != seed_hash: raise ValueError('Seed mismatch')
                cfg['sweep']['initial_state_file'] = str(seed)
                linked.write(dest/'control.json', cfg)
                start = time.perf_counter()
                with (dest/'stdout.log').open('w') as log:
                    proc = subprocess.run([str(runner), '--config', str(dest/'control.json')],
                        cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
                wall = time.perf_counter()-start
                status = dict(returncode=proc.returncode, curve=linked.rows(dest/'curve.csv'))
                passed = linked.good(status, dest, r['target_V'])
                trace = linked.rows(dest/'newton_iterations.csv')
                cost, events = trace_cost(trace, cfg['solver']['block_absolute_convergence'],
                                         cfg['solver']['carrier_row_convergence']['eps_row'])
                delta = linked.state_difference(dest/'state.csv', src/'state.csv') if passed else None
                equivalent = passed and max(delta[k]['max_absolute'] for k in ('psi', 'phin', 'phip')) <= 1e-8
                item = dict(directory=str(dest), passed=passed, potential_equivalent_1e8=equivalent,
                    state_difference=delta, seed_sha256=seed_hash, wall_seconds=wall,
                    initial_residual=float(next(x for x in trace if x['event']=='initial')['residual_norm']),
                    config_sha256=linked.digest(dest/'control.json'), events=events, **cost)
                record['variants'][variant] = item
                linked.write(a.output/'summary.json', result)
                print(gate, r['target_V'], variant, cost['updates'], passed, equivalent, flush=True)
                if not passed or not equivalent: raise RuntimeError('Original gate or state comparison failed')
    result['status'] = 'pass'
    linked.write(a.output/'summary.json', result)


if __name__ == '__main__':
    main()
