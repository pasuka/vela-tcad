"""Serial Release qualification of the shared Halley default on Linux.

Requires staged, hashed D4/D5 bundles and the prior screening qualification.
The linked sweep's Windows CPU reader is replaced only in this Linux worker
adapter by /proc process CPU counters; physics, targets and gates are unchanged.
No benchmark speedup is inferred from these correctness runs.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

from run_templates_ldmos_screening_candidates import read, save, digest


def linux_worker_cpu(proc):
    # /proc stat fields 14 and 15, after the final ')' enclosing comm.
    fields = Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')', 1)[1].split()
    ticks = os.sysconf('SC_CLK_TCK')
    user, kernel = int(fields[11])/ticks, int(fields[12])/ticks
    return dict(user=user, kernel=kernel, total=user+kernel)


def relocate(value, old, new):
    if isinstance(value, dict):
        return {k: relocate(v, old, new) for k, v in value.items()}
    if isinstance(value, list):
        return [relocate(v, old, new) for v in value]
    if isinstance(value, str):
        converted = value.replace(old.replace('/', '\\'), str(new)).replace(old, str(new))
        return converted.replace('\\', '/') if converted != value else value
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('runner', 'prior', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--d0-evidence', type=Path,
                        help='Reuse completed, native-qualified D0 evidence after a later staging failure')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runner, prior, out = (p.resolve() for p in (args.runner, args.prior, args.output))
    out.mkdir(parents=True, exist_ok=False)
    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', OMP_DYNAMIC='FALSE',
                      VELA_LINEAR_THREADS='1', VELA_LINEAR_SOLVER='umfpack')
    # Standard OpenBLAS startup control works with or without Vela's optional
    # experimental runtime-control API. Do not request an unavailable API.
    os.environ.pop('VELA_BLAS_THREADS', None)
    blas_threads = int(subprocess.check_output([sys.executable, '-c',
        'import ctypes; print(ctypes.CDLL("libopenblas.so.0").openblas_get_num_threads())'], text=True))
    if blas_threads != 1:
        raise ValueError('OpenBLAS startup environment did not select one thread')
    frozen = {str(p): digest(p) for p in [runner, Path(__file__), *root.glob('include/vela/physics/Ial*.h'),
                root/'src/physics/IalMobility.cpp', root/'src/equation/IalTransport.cpp',
                root/'src/simulation/ElectrothermalSimulation.cpp']}
    status = dict(status='running', pid=os.getpid(), started=time.time(), runs=[], joint={}, comparisons=[],
                  frozen_sha256=frozen, expected_curves=7, expected_exact_points=217)
    status['blas_control'] = dict(mode='OPENBLAS_NUM_THREADS=1', fresh_process_threads=blas_threads)
    d0_root = args.d0_evidence.resolve() if args.d0_evidence else out
    if args.d0_evidence:
        previous = read(d0_root/'summary.json')
        if previous.get('joint', {}).get('D0') != 'pass' or len(previous.get('comparisons', [])) != 62:
            raise ValueError('D0 qualification and all state comparisons required')
        if previous['frozen_sha256'].get(str(runner)) != digest(runner):
            raise ValueError('D0 evidence belongs to another executable')
        for p, h in previous['frozen_sha256'].items():
            if digest(Path(p)) != h:
                raise ValueError('Frozen D0 evidence source changed: '+p)
        frozen.update(previous['frozen_sha256'])
        frozen[str(d0_root/'summary.json')] = digest(d0_root/'summary.json')
        status['runs'] = [r for r in previous['runs'] if r['name'] in ('d0_vg4','d0_vg8')]
        status['comparisons'] = previous['comparisons']
        status['reused_d0_evidence'] = str(d0_root)
    publish = lambda: save(out/'summary.json', status)
    def run(name, command):
        status['current'] = name
        start = time.perf_counter()
        with (out/(name+'.log')).open('w') as log:
            child = subprocess.Popen(list(map(str, command)), cwd=root, stdout=log, stderr=subprocess.STDOUT)
            status['child_pid'] = child.pid
            publish()
            rc = child.wait()
        status.pop('child_pid', None)
        row = dict(name=name, returncode=rc, wall_seconds=time.perf_counter()-start)
        status['runs'].append(row)
        publish()
        if rc:
            raise RuntimeError(f'{name}: exit {rc}')
        return row
    publish()
    try:
        from run_templates_ldmos_screening_repeats import inspect_curve, compare_states
        from electrothermal_state import read_bound_record
        from analyze_templates_ldmos_d0 import analyze as analyze_d0
        for gate in (() if args.d0_evidence else (4, 8)):
            name = f'd0_vg{gate}'
            folder = out/name
            folder.mkdir()
            original = prior/f'evidence/full_20260925/vg{gate}_halley'
            cfg = read(original/'input.json')
            cfg.pop('diagnostic_ialmob_screening_method', None)
            cfg['mobility_SI']['ialmob'].pop('screening_method', None)
            save(folder/'input.json', cfg)
            deck = read(original/'deck.json')
            deck.update(input_file=str(folder/'input.json'), output_directory=str(folder/'results'), resume=False)
            save(folder/'deck.json', deck)
            row = run(name, [runner, '--config', folder/'deck.json'])
            ledger, states, totals = inspect_curve(folder, 'halley')
            row['totals'] = totals
            old = read(original/'results/ledger.json')
            for index, (point, state) in enumerate(zip(old['exact_points'], states)):
                if point['bias_V'] != ledger['exact_points'][index]['bias_V']:
                    raise ValueError('Bias changed')
                errors = compare_states(read_bound_record(point['result']), state)
                status['comparisons'].append(dict(gate=gate, point=index, errors=errors))
            publish()
        contract = root/'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json'
        report = analyze_d0(prior/'evidence/native_20260925',
                            {g: d0_root/f'd0_vg{g}/results' for g in (4, 8)}, read(contract), 15)
        save(out/'d0_joint.json', report)
        status['joint']['D0'] = report['status']
        if report['status'] != 'pass':
            raise ValueError('Original D0 gates failed')
        # Recover the frozen G3 control and convert only its old persistence input.
        from migrate_state_seed_to_hdf5 import migrate
        g3 = out/'g3'
        g3.mkdir()
        oldroot = 'D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a'
        oldout = oldroot+'/reference_staging/templates_ldmos_joint_20260914/regression_r6/g3_auger'
        cfg = relocate(relocate(read(root/'g3_original/control.json'), oldout, g3), oldroot, root)
        seed = g3/'seed.h5'
        save(g3/'conversion.json', migrate(Path(cfg['sweep']['initial_state_file']), Path(cfg['mesh_file']), seed, 'unit_scaling', 0.))
        cfg['state_format'] = 'hdf5'
        cfg['sweep'].update(initial_state_file=str(seed), write_state_file=str(g3/'state.h5'))
        save(g3/'control.json', cfg)
        run('g3_full', [runner, '--config', g3/'control.json'])
        reference = relocate(read(root/'g3_original/qualification.json')['reference'], oldroot, root)
        run('g3_acceptance', [sys.executable, root/'scripts/analyze_templates_ldmos_g3_idvg.py',
             '--reference', reference, '--candidate', cfg['output_csv'], '--balance', cfg['sweep']['diagnostics']['srh_balance']['csv_file'],
             '--output-json', g3/'qualification.json', '--output-md', g3/'qualification.md'])
        report = read(g3/'qualification.json')
        if report['status'] != 'pass' or report['metrics']['exact_shared_points'] != 31:
            raise ValueError('Original G3 gates failed')
        status['joint']['G3'] = report['status']
        save(out/'runner.json', dict(runner=str(runner), runner_sha256=digest(runner), linear_solver='umfpack',
              backend='UMFPACK Linux Release 1/1', frozen_sources=frozen))
        from analyze_templates_ldmos_stage4_d5 import analyze as analyze_joint
        for profile in ('D5', 'D4'):
            bundlepath = root/f'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_{profile.lower()}_auger_no_generation_inputs.json'
            bundle = read(bundlepath)
            for gate in (4, 8):
                run(f'{profile.lower()}_vg{gate}', [sys.executable, Path(__file__), '--linked',
                    '--workspace', root, '--physics-profile', profile, '--bundle', bundlepath,
                    '--manifest', out/'runner.json', '--output', out/f'{profile.lower()}_vg{gate}',
                    '--gate', gate, '--points', 31, '--worker', '--reuse-linear-analysis', '--linear-solver', 'umfpack'])
            report = analyze_joint({f'Vg{g}': root/bundle['references'][str(g)] for g in (4, 8)},
                {f'Vg{g}': out/f'{profile.lower()}_vg{g}/score/curve.csv' for g in (4, 8)},
                {f'Vg{g}': out/f'{profile.lower()}_vg{g}/score/terminal_balance.csv' for g in (4, 8)},
                out/f'{profile.lower()}_joint', physics_profile=profile)
            if report['engineering']['status'] != 'pass' or report['final']['status'] != 'pass':
                raise ValueError('Original '+profile+' gates failed')
            status['joint'][profile] = 'pass'
            publish()
        if any(digest(Path(p)) != h for p, h in frozen.items()):
            raise ValueError('Frozen source or executable changed')
        status.update(status='completed', finished=time.time())
        publish()
    except BaseException as error:
        status.update(status='failed', error=repr(error), traceback=traceback.format_exc(), finished=time.time())
        publish()
        raise


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--linked':
        if sys.platform != 'linux' or '--worker' not in sys.argv:
            raise ValueError('This CPU adapter requires a live Linux DC worker')
        import run_templates_ldmos_linked_d5 as linked
        linked.cpu_times = linux_worker_cpu
        del sys.argv[1]
        linked.main()
    else:
        main()
