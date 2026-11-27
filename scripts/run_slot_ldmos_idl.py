"""Bounded cold-start SLOT IDL validation; candidates never count as accepted states.

Input configuration paths are relative to the worktree root. Native exported
mesh/geometry resources remain prerequisites. All output goes to a new directory.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
MATCHED_GATES = [1.98506236911173, 2.58506236911173, 2.59506236911173,
                 2.60506236911173, 2.61506236911173, 2.62506236911173,
                 2.68506236911173, 2.98506236911173, 3.98506236911173, 5.0]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    tmp.replace(path)


def qualify(row, target):
    """Require the complete strict solve's gates, not a successful probe exit."""
    try:
        return (row['converged'] == '1' and row['solver_method'] == 'newton'
                and row['carrier_row_violations'] == '0'
                and row['global_continuity_closure_satisfied'] == '1'
                and row['qf_bounds_violations'] == '0'
                and 0 <= float(row['carrier_row_max_ratio']) <= .001
                and math.isfinite(float(row['current_total_A_per_um']))
                and abs(float(row['bias_V']) - target) < 1e-12)
    except (KeyError, ValueError, TypeError):
        return False


def resource_paths(c):
    return [c['mesh_file'], c['materials_file'], c['node_doping_file'],
            c['mesh_geometry']['external_averagebox_couples_file'],
            c['solver']['mobility']['ialmob']['geometry_file'],
            c['solver']['mobility']['ialmob']['avalanche_geometry_file']]


class Campaign:
    def __init__(self, args):
        self.out = args.output.resolve()
        self.out.mkdir(parents=True, exist_ok=False)
        self.c = json.loads(args.config.read_text())
        def absolute(value):
            if isinstance(value, dict):
                return {k: absolute(v) for k, v in value.items()}
            if isinstance(value, list):
                return [absolute(v) for v in value]
            if isinstance(value, str) and value.startswith(('build/', 'configs/')):
                return str(ROOT / value)
            return value
        self.c = absolute(self.c)
        self.resources = {p: sha(p) for p in resource_paths(self.c)}
        self.resources[str(args.config.resolve())] = sha(args.config)
        for name in ['run_slot_ldmos_idl.py', 'slot_poly_equilibrium.py',
                     'slot_state_archive.py', 'state_archive.py']:
            self.resources[str(ROOT / 'scripts' / name)] = sha(ROOT / 'scripts' / name)
        shutil.copy2(args.runner, self.out / 'runner.exe')
        self.resources[str(self.out / 'runner.exe')] = sha(args.runner)
        self.env = os.environ.copy()
        self.env['PATH'] = 'D:/msys64/ucrt64/bin;D:/msys64/usr/bin;' + self.env['PATH']
        if self.env.get('VELA_LINEAR_SOLVER'):
            raise ValueError('Unexpected VELA_LINEAR_SOLVER override')
        self.t0 = time.monotonic()
        self.budget, self.child_cap = args.budget, args.child_cap
        self.s = dict(completed=False, goal_reached=False, records=[], cases=[],
                      budget_s=args.budget, child_cap_s=args.child_cap,
                      cold_start=True, imported_states=0, resources=self.resources,
                      policy=dict(density_repairs=2, max_halvings=3,
                                  gate_max_step_V=.1, matched_gates=MATCHED_GATES))
        save(self.out / 'base.json', self.c)
        self.persist()

    def persist(self):
        self.s['elapsed_s'] = time.monotonic() - self.t0
        save(self.out / 'summary.json', self.s)

    def execute(self, stage, target, state=None, density=False, poly=False):
        if time.monotonic() - self.t0 > self.budget - self.child_cap - 5:
            raise TimeoutError('insufficient_budget_for_bounded_child')
        from slot_poly_equilibrium import write_poly_equilibrium_restart
        d = self.out / f'{len(self.s["cases"]):03d}_{stage}_{"density" if density else "strict"}'
        d.mkdir()
        c = copy.deepcopy(self.c)
        c['output_csv'] = str(d / 'iv.csv')
        c['solver']['carrier_row_convergence']['diagnostic_csv'] = str(d / 'rows.csv')
        for contact in c['contacts']:
            contact['bias'] = .1 if stage == 'gate' and contact['name'] == 'drain' else 0.
        c['sweep'] = dict(mode='iv', contact='gate' if stage == 'gate' else 'drain',
                          current_contact='drain', bias_points=[target], start=target,
                          stop=target, step=.005, write_vtk=False,
                          write_state_file=str(d / 'state.h5'),
                          diagnostics={'newton_history': dict(enabled=True,
                          attempts_csv_file=str(d / 'attempts.csv'),
                          iterations_csv_file=str(d / 'iterations.csv'),
                          rejected_state_directory=str(d / 'rejected_states'))})
        parent_hash = sha(state) if state else None
        if poly:
            provenance = write_poly_equilibrium_restart(c, state, d / 'candidate.h5',
                                                       ROOT, target)
            save(d / 'poly.json', provenance)
            state = d / 'candidate.h5'
        if state:
            c['sweep']['initial_state_file'] = str(state)
        else:
            c['sweep']['initialization'] = dict(mode='poisson_block',
                diagnostic_csv=str(d / 'initialization.csv'), write_state_file=str(d / 'poisson.h5'))
        if density:
            c.pop('sweep')
            c.update(simulation_type='newton_block_step_probe', block_modes=['carrier_density'],
                     state_file=str(state), output_state_prefix=str(d / 'trial'))
            for contact in c['contacts']:
                if contact['name'] == ('gate' if stage == 'gate' else 'drain'):
                    contact['bias'] = target
        save(d / 'config.json', c)
        t = time.monotonic()
        with (d / 'stdout.log').open('wb') as stdout, (d / 'stderr.log').open('wb') as stderr:
            try:
                rc = subprocess.run([str(self.out / 'runner.exe'), '--config', str(d / 'config.json')],
                    cwd=ROOT, env=self.env, stdout=stdout, stderr=stderr,
                    timeout=self.child_cap).returncode
            except subprocess.TimeoutExpired:
                rc = 124
        result = dict(case=d.name, stage=stage, target_V=target, returncode=rc,
                      elapsed_s=time.monotonic()-t, density_probe=density,
                      parent_sha256=parent_hash, config_sha256=sha(d / 'config.json'), qualified=False)
        if not density and rc == 0 and (d / 'state.h5').exists() and (d / 'iv.csv').exists():
            with (d / 'iv.csv').open() as f:
                rows = list(csv.DictReader(f))
            if rows:
                result['iv'] = rows[-1]
                result['qualified'] = qualify(rows[-1], target)
        save(d / 'exit.json', result)
        self.s['cases'].append({k:v for k,v in result.items() if k != 'iv'})
        self.persist()
        print(json.dumps(self.s['cases'][-1]), flush=True)
        return d, result

    def step(self, stage, target, parent):
        d, result = self.execute(stage, target, parent, poly=stage == 'gate')
        for _ in range(2):
            if result['qualified']:
                break
            rejected = sorted((d / 'rejected_states').glob('*_final.h5'))
            if not rejected:
                break
            b, probe = self.execute(stage, target, rejected[-1], density=True)
            if probe['returncode'] or not (b / 'trial_carrier_density.h5').exists():
                break
            d, result = self.execute(stage, target, b / 'trial_carrier_density.h5', poly=stage == 'gate')
        if result['qualified']:
            record = dict(case=d.name, stage=stage, target_V=target,
                state_sha256=sha(d / 'state.h5'), config_sha256=sha(d / 'config.json'), iv=result['iv'])
            self.s['records'].append(record)
            self.persist()
            return d / 'state.h5'
        return None

    def advance(self, stage, value, target, parent, depth=0):
        state = self.step(stage, target, parent)
        if state:
            return state
        if depth >= 3:
            raise RuntimeError(f'bounded_recovery_exhausted:{stage}:{target}')
        middle = (value + target) / 2
        state = self.advance(stage, value, middle, parent, depth + 1)
        return self.advance(stage, middle, target, state, depth + 1)

    def run(self):
        try:
            state = self.step('equilibrium', 0., None)
            if state is None:
                raise RuntimeError('cold_equilibrium_failed')
            value = 0.
            for target in [.0001, .001, .005, .01] + [i*.005 for i in range(3, 21)]:
                state = self.advance('drain', value, target, state)
                value = target
            value, h = 0., .0005
            while value < 5. - 1e-12:
                target = min(5., value+h, *[x for x in MATCHED_GATES if x > value+1e-12])
                state = self.advance('gate', value, target, state)
                value, h = target, min(.1, h*1.5)
            self.s['goal_reached'] = True
            self.s['stop_reason'] = 'full_idl_path_qualified'
        except Exception as exc:
            self.s['stop_reason'] = str(exc)
            raise
        finally:
            self.s['completed'] = True
            self.s['resources_unchanged'] = all(sha(p) == digest for p,digest in self.resources.items())
            self.persist()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=ROOT/'configs/regression/slot_ldmos_idl_bennett_control.json')
    p.add_argument('--runner', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--budget', type=int, default=7200)
    p.add_argument('--child-cap', type=int, default=180)
    args = p.parse_args()
    if args.child_cap <= 0 or args.budget <= args.child_cap + 5:
        p.error('positive child cap and sufficient bounded total budget required')
    Campaign(args).run()


if __name__ == '__main__':
    main()
