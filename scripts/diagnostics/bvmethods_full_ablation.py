"""Reproducible P0 ablation / qualified P1 supervisor; all host paths are inputs.

Run in the isolated Vela environment with an explicit request JSON. Every solver
is launched in its own process group. P0 negative controls are retained; execution
success never substitutes for numeric qualification. No production defaults change.
"""
from __future__ import annotations
import concurrent.futures
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def put(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False))
    temp.replace(path)


def read(path):
    return json.loads(Path(path).read_text())


def rows(path):
    if not Path(path).exists():
        return []
    with Path(path).open() as stream:
        return list(csv.DictReader(stream))


def accepted(path):
    return [r for r in rows(path) if r.get('converged') == '1' and r.get('outer_voltage_V')]


def gate(row):
    g = json.loads(row['validation_diagnostics'])
    physical = (g['physical_residuals_finite'] and g['terminal_kcl_satisfied']
                and (not g['carrier_row_enforced'] or g['carrier_row_satisfied'])
                and (not g['global_continuity_enforced'] or g['global_continuity_satisfied'])
                and g['actual_load_line_satisfied'])
    if float(row['outer_voltage_V']) == 0:
        return physical
    return (physical and g['density_update_mode'] == 'enforce_joint'
            and g['density_update_valid'] and g['density_update_rms'] < 1
            and g['coupled_device_residual_inf'] <= 1e-8
            and abs(g['actual_load_line_residual_V']) <= 1e-6)


def interpolate(data, current):
    for a, b in zip(data, data[1:]):
        ia, ib = abs(float(a['current_total_A_per_um'])), abs(float(b['current_total_A_per_um']))
        if ia <= current <= ib and ib > ia:
            va, vb = float(a['inner_voltage_V']), float(b['inner_voltage_V'])
            return {'current_A_per_um': current, 'inner_V': va+(vb-va)*(current-ia)/(ib-ia),
                    'bracket_inner_V': [va, vb], 'bracket_current_A_per_um': [ia, ib],
                    'bracket_width_V': abs(vb-va)}
    return None


def hdf_identity(path):
    import h5py
    identity = {}
    with h5py.File(path) as f:
        def visitor(name, item):
            if isinstance(item, h5py.Dataset):
                a = item[()]
                identity[name] = {'dtype': str(a.dtype), 'shape': list(a.shape),
                                  'sha256': hashlib.sha256(a.tobytes() if a.dtype.kind != 'O' else repr(a.tolist()).encode()).hexdigest()}
        f.visititems(visitor)
    return identity


class Supervisor:
    def __init__(self, request):
        self.req = request
        self.out = Path(request['out'])
        self.out.mkdir(parents=True, exist_ok=False)
        self.binary = Path(request['binary'])
        self.standard = Path(request['standard_binary'])
        self.binary_sha = sha(self.binary)
        self.lock = threading.Lock()
        self.state = {'phase': 'P0', 'finished': False, 'groups': {}, 'parity': {},
                      'binary_sha256': self.binary_sha, 'source_commit': request['source_commit'],
                      'full_bv_started': False, 'p0': {}}
        self.env = os.environ.copy()
        self.env.update(OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                        MKL_THREADING_LAYER='SEQUENTIAL', VELA_LINEAR_SOLVER='pardiso_mkl')
        self.env.pop('VELA_CONSTRAINT_ELIMINATION', None)
        self.save()

    def save(self):
        with self.lock:
            put(self.out/'status.json', self.state)

    def process(self, folder, args, timeout, progress=None):
        assert sha(self.binary) == self.binary_sha
        started = time.time()
        with (folder/'solver.log').open('wb') as log:
            p = subprocess.Popen(args, env=self.env, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True, cwd=self.req['source_dir'])
            timed_out = False
            try:
                while p.poll() is None:
                    if progress:
                        progress(started, started+timeout)
                    remaining = started+timeout-time.time()
                    if remaining <= 0:
                        timed_out = True
                        os.killpg(p.pid, signal.SIGTERM)
                        try:
                            p.wait(5)
                        except subprocess.TimeoutExpired:
                            os.killpg(p.pid, signal.SIGKILL)
                            p.wait()
                        break
                    try:
                        p.wait(min(20, remaining))
                    except subprocess.TimeoutExpired:
                        pass
            finally:
                if p.poll() is None:
                    os.killpg(p.pid, signal.SIGTERM)
                    try:
                        p.wait(5)
                    except subprocess.TimeoutExpired:
                        os.killpg(p.pid, signal.SIGKILL)
                        p.wait()
        return {'solver_exit_code': p.returncode, 'timed_out': timed_out,
                'started_unix': started, 'deadline_unix': started+timeout,
                'seconds': time.time()-started}

    def config(self, template, folder, group):
        cfg = copy.deepcopy(template)
        cfg['solver']['population_derivative_recovery'] = 'full' if group == 'G1' else 'automatic'
        cfg['solver']['performance_profiling'] = {'enabled': True, 'json_file': str(folder/'performance.json')}
        cfg['output_csv'] = str(folder/'iv.csv')
        cfg['runtime_log']['file'] = str(folder/'runtime.log')
        sweep = cfg['sweep']
        sweep['vtk_prefix'] = str(folder/'state')
        sweep['write_state_file'] = str(folder/'state.h5')
        sweep['write_state_every_accepted_step_prefix'] = str(folder/'checkpoints/state')
        sweep['external_circuit']['coupled_update_convergence']['diagnostic_csv'] = str(folder/'updates.csv')
        sweep['external_circuit']['coupled_update_convergence']['mode'] = 'enforce_joint'
        sweep['external_circuit']['load_line_convergence_mode'] = 'enforce'
        if 'initialization' in sweep:
            sweep['initialization']['write_state_file'] = str(folder/'poisson.h5')
        return cfg

    def sweep(self, name, template, group, packed=None, audit=True, standard=False, folder=None, timeout=1200, progress=None):
        folder = folder or self.out/name/group
        folder.mkdir(parents=True, exist_ok=True)
        cfg = self.config(template, folder, group)
        put(folder/'config.json', cfg)
        options = {'guard': group == 'G2', 'report_model_quality': audit, 'out': str(folder/'diagnostics')}
        if packed:
            options.update(packed_binary=packed['file'], inner=packed['inner'])
        put(folder/'experiment.json', options)
        args = [str(self.standard if standard else self.binary), '--config', str(folder/'config.json')]
        if not standard:
            args += ['--experiment', str(folder/'experiment.json')]
        result = self.process(folder, args, timeout, progress)
        data = accepted(folder/'iv.csv')
        targets = cfg['sweep']['bias_points']
        result.update(accepted_points=len(data), all_accepted_gates=all(gate(x) for x in data),
                      requested_reached=all(any(abs(float(x['outer_voltage_V'])-t) <= 1e-11 for x in data) for t in targets),
                      attempts=len(rows(folder/'diagnostics/attempts.csv')) if not standard else None,
                      passed=False)
        if (folder/'diagnostics/initial.bin').exists():
            result.update(initial_sha256=sha(folder/'diagnostics/initial.bin'),
                          initial_rhs_sha256=sha(folder/'diagnostics/initial_rhs.bin'))
        result['passed'] = (result['solver_exit_code'] == 0 and not result['timed_out']
                            and result['all_accepted_gates'] and result['requested_reached'])
        if data:
            result['last_accepted'] = {k: data[-1][k] for k in
                                      ['outer_voltage_V', 'inner_voltage_V', 'current_total_A_per_um', 'iterations']}
        if (folder/'diagnostics/partition.json').exists():
            partition = read(folder/'diagnostics/partition.json')
            assert partition['other_ns'] >= 0
            assert sum(partition[k] for k in ['attempt_ns', 'record_point_ns', 'other_ns']) == partition['dc_sweep_total_ns']
            result['partition'] = partition
        put(folder/'qualification.json', result)
        return result

    def checkpoint(self, name, spec, group, audit=True):
        folder = self.out/name/group
        folder.mkdir(parents=True)
        cfg = read(spec['config'])
        cfg['solver']['population_derivative_recovery'] = 'full' if group == 'G1' else 'automatic'
        put(folder/'config.json', cfg)
        req = {k: v for k, v in spec.items() if k in ['packed_binary', 'inner', 'targets', 'budget']}
        req.update(config=str(folder/'config.json'), out=str(folder/'probe'),
                   guard=group == 'G2', report_model_quality=audit)
        put(folder/'request.json', req)
        result = self.process(folder, [str(self.binary), '--checkpoint', str(folder/'request.json')], 1200)
        probe = folder/'probe/report.json'
        if probe.exists():
            report = read(probe)
            result.update(passed=report['passed'], accepted_points=sum(p['strict_gates'] for p in report['points']),
                          attempts=len(report['points']), points=[{**{k:v for k,v in p.items() if k!='trace'},
                              'model_reviews':sum(bool(t['model_review']) for t in p['trace']),
                              'rejections':sum(bool(t['rejected']) for t in p['trace']),
                              'last_model_observation':next((t for t in reversed(p['trace']) if t['model_report_valid']),None)}
                              for p in report['points']],
                          initial_sha256=sha(folder/'probe/initial.bin'),
                          initial_rhs_sha256=sha(folder/'probe/initial_rhs.bin'))
        else:
            result['passed'] = False
        put(folder/'qualification.json', result)
        return result

    def parity(self, window, template, group):
        folder = self.out/'parity_runs'/window/group
        standard = self.sweep(window, template, group, standard=True, audit=False, folder=folder)
        reference = {'iv': rows(folder/'iv.csv'), 'updates': rows(folder/'updates.csv'),
                     'state': hdf_identity(folder/'state.h5')}
        put(folder/'standard_identity.json', reference)
        for audit in [False, True]:
            result = self.sweep(window, template, group, audit=audit, folder=folder)
            current = {'iv': rows(folder/'iv.csv'), 'updates': rows(folder/'updates.csv'),
                       'state': hdf_identity(folder/'state.h5')}
            same = reference == current
            self.state['parity'][window+'_'+group+('_rho_on' if audit else '_rho_off')] = {
                'same_trajectory_and_state_dataset_bits': same, 'standard_passed': standard['passed'],
                'experiment_passed': result['passed']}
            self.save()
            assert same and standard['passed'] and result['passed'], 'entry/report-only parity failed'

    def run_p0(self):
        cases = self.req['cases']
        for name in ['W1', 'W2']:
            self.state['current'] = 'parity_'+name
            self.save()
            for group in ['G3', 'G1']:
                self.state['current'] = 'parity_'+name+'_'+group
                self.save()
                self.parity(name, read(cases[name]['config']), group)
        for index, (name, spec) in enumerate(cases.items()):
            if name == 'B_i1':
                continue
            self.state['p0'][name] = {}
            order = ['G1', 'G2', 'G3']
            order = order[index % 3:]+order[:index % 3]
            for group in order:
                self.state['current'] = name+'_'+group
                self.save()
                if spec['kind'] == 'sweep':
                    result = self.sweep(name, read(spec['config']), group, spec.get('packed'))
                else:
                    result = self.checkpoint(name, spec, group)
                self.state['p0'][name][group] = result
                self.save()
        if 'B_i1' in cases:
            self.state['current'] = 'B_i1_report_only'
            self.save()
            self.state['B_i1_observation'] = self.checkpoint('B_i1', cases['B_i1'], 'G1')
        self.state['initial_identity_agreement'] = {name:len({(c[g].get('initial_sha256'),c[g].get('initial_rhs_sha256')) for g in ['G1','G2','G3']})==1 for name,c in self.state['p0'].items()}
        self.save()
        assert all(self.state['initial_identity_agreement'].values()), 'group initial state/RHS identities differ'
        complete = {g: all(c[g]['passed'] for c in self.state['p0'].values()) for g in ['G1', 'G2', 'G3']}
        costs = {g: {'seconds': sum(c[g]['seconds'] for c in self.state['p0'].values()),
                     'attempts': sum(c[g].get('attempts', 0) for c in self.state['p0'].values())}
                 for g in ['G1', 'G2', 'G3']}
        self.state.update(p0_complete=True, p0_all_windows_passed=complete, p0_costs=costs)
        # Report low-rho failures without asserting a unique cause. Do not silently add candidate selection.
        if complete['G1'] and (not complete['G2'] or
                              (costs['G1']['seconds'] <= costs['G2']['seconds'] and
                               costs['G1']['attempts'] <= costs['G2']['attempts'])):
            candidates = ['G1']
        elif complete['G1'] and complete['G2']:
            candidates = ['G1', 'G2']
        else:
            candidates = []  # Failed-G1 cases require same-state rho/linear-gate review before P1 selection.
        self.state['provisional_p1_candidates'] = candidates
        self.save()
        return candidates

    def p1(self, group):
        folder = self.out/'P1'/group
        folder.mkdir(parents=True)
        def progress(start, deadline):
            snapshot = {'phase': 'full_scan', 'simulation_started_unix': start,
                        'simulation_deadline_unix': deadline, 'elapsed_seconds': time.time()-start}
            try:
                data = accepted(folder/'iv.csv')
                updates = rows(folder/'diagnostics/newton.csv')
                snapshot['accepted_points'] = len(data)
                if data:
                    snapshot['last_accepted'] = {k: data[-1][k] for k in
                                                 ['outer_voltage_V', 'inner_voltage_V', 'current_total_A_per_um', 'iterations']}
                if updates:
                    snapshot['last_completed_attempt'] = updates[-1]
            except (ValueError, KeyError, csv.Error):
                pass  # A row can still be in flight; retain explicit last-complete snapshots next time.
            put(folder/'status.json', snapshot)
        result = self.sweep('P1', read(self.req['full_config']), group, folder=folder,
                            timeout=self.req.get('full_seconds', 28800), progress=progress)
        data = accepted(folder/'iv.csv')
        stop = bool(data and abs(float(data[-1]['current_total_A_per_um'])) >= 1.3e-3)
        bv = interpolate(data, 1e-4)
        relative = abs(bv['inner_V']/6.673833-1) if bv else None
        result.update(stop_current_reached=stop, bv_at_1e_4=bv, bv_reference_V=6.673833,
                      bv_relative_difference=relative,
                      current_table=[interpolate(data, v) for v in [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]])
        result['full_bv_completed'] = bool(result['solver_exit_code'] == 0 and not result['timed_out']
                                           and stop and result['all_accepted_gates'] and relative is not None and relative <= .005)
        result['phase'] = 'finished'
        put(folder/'status.json', result)
        put(folder/'exit.json', {'solver_exit_code': result['solver_exit_code'], 'worker_exit_code': 0 if result['full_bv_completed'] else 4})
        with self.lock:
            self.state['groups'][group] = result
        self.save()
        return result

    def execute(self):
        try:
            candidates = self.run_p0()
            put(self.out/'p0_summary.json', self.state)
            if candidates and self.req.get('allow_p1', False):
                self.state.update(phase='P1', full_bv_started=True)
                self.save()
                with concurrent.futures.ThreadPoolExecutor(max_workers=len(candidates)) as pool:
                    list(pool.map(self.p1, candidates))
            self.state.update(phase='finished', finished=True, worker_exit_code=0,
                              needs_p1_selection_review=not bool(candidates))
            self.save()
            put(self.out/'exit.json', {'worker_exit_code': 0})
        except Exception as error:
            self.state.update(phase='failed', finished=True, worker_exit_code=1, error=str(error))
            self.save()
            put(self.out/'exit.json', {'worker_exit_code': 1})
            raise


if __name__ == '__main__':
    Supervisor(read(sys.argv[1])).execute()
