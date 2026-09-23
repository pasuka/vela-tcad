#!/usr/bin/env python3
"""Independently recount linked D5 work and recheck every accepted state."""
from __future__ import annotations
import argparse
import math
from pathlib import Path

from run_templates_ldmos_linked_d5 import read, rows, digest, good, write
from summarize_templates_ldmos_idvd_ablation import read_curve


def audit(directory: Path) -> dict:
    plan = read(directory/'plan.json')
    ledger = read(directory/'fixed/ledger.json')
    if ledger['status'] != 'completed':
        raise ValueError('A partial or failed run cannot earn a completed audit')
    for name, expected in plan['frozen_files'].items():
        if digest(Path(name)) != expected:
            raise ValueError(f'Changed frozen dependency: {name}')
    updates = 0
    worker_ids = []
    for run in ledger['runs']:
        case = directory/run['case']
        status = read(case/'status.json')
        config = read(case/'control.json')
        trace = rows(case/'newton_iterations.csv')
        profile = read(case/'performance_profile.json')
        count = sum(row['event'] == 'accepted_iteration' for row in trace)
        assert count == run['Newton_updates'] == sum(int(row['newton_iterations']) for row in status['attempts'])
        assert count == profile['counters'].get('newton.updates', 0)
        assert status['config_sha256'] == digest(case/'control.json')
        assert status['seed_sha256'] == run['parent_sha256'] == digest(Path(run['parent_state']))
        assert status['runner_sha256'] == plan['runner_sha256']
        assert config['solver']['block_absolute_convergence'] == plan['original_blocks']
        assert config['solver']['carrier_row_convergence']['eps_row'] == 1e-8
        assert config['solver']['quasi_fermi_update_limit_V'] == run['cap_V']
        if plan.get('execution_mode') == 'dc_worker':
            worker_ids.append((status['pid'],status['worker_request_id']))
        contacts = {c['name']: c['bias'] for c in config['contacts']}
        assert abs(contacts['gate']-contacts['source']-plan['gate_V']) < 1e-13
        assert abs(contacts['drain']-contacts['source']-run['target_V']) < 1e-13
        updates += count
    accepted = ledger['exact_points'] + ledger['transfers']
    accepted += [e['accepted'] for e in ledger['frame_events'] if e['status'] == 'qualified']
    unique = {}
    mesh_bindings = {}
    for point in accepted:
        case = directory/point['case']
        assert good(read(case/'status.json'), case, point['bias_V'])
        assert digest(Path(point['state'])) == point['sha256']
        state_path = Path(point['state'])
        if state_path.suffix == '.h5':
            import state_archive
            cfg = read(case/'control.json')
            mesh_path = Path(cfg['mesh_file'])
            scaling = 1e-6 if cfg.get('scaling',{}).get('mode') == 'unit_scaling' else 1.
            key = (mesh_path, scaling)
            if key not in mesh_bindings:
                mesh = read(mesh_path)
                mesh_bindings[key] = (len(mesh['nodes']), state_archive.mesh_identity(mesh,scaling))
            count,identity = mesh_bindings[key]
            _,metadata = state_archive.read(state_path,count,identity)
            if metadata['mode'] != 'dd' or metadata['potential_origin_V'] != cfg.get('potential_origin_V',0.):
                raise ValueError('Checkpoint mode/frame differs from the accepted control')
        unique[point['state']] = point['sha256']
    assert updates == ledger['total_Newton_updates']
    if worker_ids:
        assert len({pid for pid,_ in worker_ids}) == 1
        assert [request for _,request in worker_ids] == list(range(1,len(ledger['runs'])+1))
    assert [p['bias_V'] for p in ledger['exact_points']] == plan['targets_V']
    reference = dict(read_curve(Path(plan['reference'])))
    errors = []
    for point in ledger['exact_points']:
        row = read(directory/point['case']/'status.json')['curve'][0]
        current = float(row['current_total_A_per_um'])
        target = reference[point['bias_V']]
        errors.append(dict(bias_V=point['bias_V'], current_A_per_um=current,
                           reference_A_per_um=target,
                           relative_error_percent=100*abs(current-target)/abs(target) if point['bias_V'] != 0. and target else None))
    return dict(integrity_pass=True, gate_V=plan['gate_V'], exact_points=len(errors),
                updates=updates, children=len(ledger['runs']), solver_requests=len(ledger['runs']),
                process_count=1 if worker_ids else len(ledger['runs']),transfers=len(ledger['transfers']),
                rollbacks=len(ledger['rollbacks']), accepted_states=len(unique),
                child_wall_seconds=math.fsum(r['wall_seconds'] for r in ledger['runs']),
                child_cpu_seconds=math.fsum(r['cpu_seconds'] for r in ledger['runs']),
                point_comparison=errors)


def main():
    if not __debug__:
        raise RuntimeError('Run without Python -O')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    report = audit(args.directory.resolve())
    write(args.directory/'independent_audit.json', report)
    print({key:value for key,value in report.items() if key != 'point_comparison'})


if __name__ == '__main__':
    main()
