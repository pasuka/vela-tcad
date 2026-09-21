"""Run an identical relocatable Release package on two Windows hosts."""
import argparse
import ctypes
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

from run_templates_ldmos_linked_d5 import read, write, digest, preflight
from run_templates_ldmos_full_curve_timing import CONFIGS, run_case, check_ready
from analyze_templates_ldmos_full_curve_timing import summarize

ROOT = Path(__file__).resolve().parents[1]
IDS = ('U1', 'T1', 'L1', 'M1', 'S1')


def verify_package(root):
    manifest = read(root/'package_manifest.json')
    for relative, expected in manifest['files'].items():
        path = (root/relative).resolve()
        if not path.is_relative_to(root.resolve()) or digest(path) != expected:
            raise ValueError('Changed or unsafe package file: '+relative)
    return manifest


def schedule():
    return [dict(stage='pilot', profile='D5', points=8, round=0, gate=8,
                 config=config, key=f'd5_p8_r0_vg8_{config}') for config in IDS]


def host_info():
    class PowerStatus(ctypes.Structure):
        _fields_ = [('ac', ctypes.c_byte), ('battery_flag', ctypes.c_byte),
                    ('battery_percent', ctypes.c_byte), ('reserved', ctypes.c_byte),
                    ('life_seconds', ctypes.c_ulong), ('full_life_seconds', ctypes.c_ulong)]
    power = PowerStatus()
    ok = ctypes.windll.kernel32.GetSystemPowerStatus(ctypes.byref(power))
    return dict(hostname=platform.node(), platform=platform.platform(), python=sys.version,
                processor=os.environ.get('PROCESSOR_IDENTIFIER'), logical_processors=os.cpu_count(),
                ac_line_status=power.ac if ok else None,
                power_scheme=subprocess.check_output(['powercfg', '/getactivescheme']).decode(errors='replace'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if os.name != 'nt' or not __debug__:
        raise RuntimeError('Requires Windows and assertions')
    package = verify_package(ROOT)
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/'binary').mkdir()
    shutil.copytree(ROOT/'scripts', out/'source/scripts', ignore=shutil.ignore_patterns('__pycache__'))
    runner = ROOT/'build-release/vela_example_runner.exe'
    dependencies = {str(ROOT/p):sha for p,sha in package['files'].items()}
    for config in IDS:
        write(out/f'binary/{config}.json', dict(runner=str(runner), runner_sha256=digest(runner),
            backend=CONFIGS[config]['backend'], linear_solver=CONFIGS[config]['backend'],
            frozen_sources=dependencies, build_type='Transferred UCRT64 Release -O3 -DNDEBUG'))
    report = dict(status='prepared', schedule=schedule(), configurations={k:CONFIGS[k] for k in IDS},
                  control_config='U1', cases=[], comparisons={}, repeat_comparisons={},
                  joint_qualifications={}, host=host_info(), package_sha256=digest(ROOT/'package_manifest.json'),
                  runner_sha256=digest(runner), scope='Single-round host pilot, Vg8 first eight D5 points, reuse 1+1')
    lock = out/'batch.lock'
    with lock.open('x') as handle:
        handle.write(str(os.getpid()))
    def save():
        write(out/'summary.json', report)
    try:
        bundle = ROOT/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json'
        for config in IDS:
            preflight(bundle, ROOT, out/f'binary/{config}.json', 8)
        report['status'] = 'running'; save()
        for item in report['schedule']:
            if (out/'PAUSE').exists():
                raise KeyboardInterrupt('PAUSE requested')
            run_case(out, report, item, save)
            check_ready(out, report); save()
            write(out/'analysis.json', summarize(out))
        verify_package(ROOT)
        report['status'] = 'pass'
    except KeyboardInterrupt as error:
        report.update(status='paused', error=repr(error))
    except BaseException as error:
        report.update(status='failed', error=repr(error))
        raise
    finally:
        save(); lock.unlink()
        write(out/'analysis.json', summarize(out))
    print('HOST_PILOT_END', report['status'], flush=True)


if __name__ == '__main__':
    main()
