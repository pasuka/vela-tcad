#!/usr/bin/env python3
"""Fresh native-reference acceptance through a verified, unpacked runtime.

The acceptance controller and vela launcher both use the package's Python.
Repository scripts supply reference comparison only; no old candidate is reused.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify_package(root):
    manifest = json.loads((root / 'PACKAGE_MANIFEST.json').read_text())
    for item in manifest['files']:
        path = root / item['path']
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError('package manifest path escapes root')
        if not path.is_file() or path.stat().st_size != item['size'] or sha(path) != item['sha256']:
            raise ValueError('package fingerprint differs: ' + item['path'])
    return manifest


def runtime_environment():
    env = {key: value for key, value in os.environ.items()
           if not key.upper().startswith(('VELA_', 'PYTHON', 'FAKE_RUNNER_'))}
    env['PATH'] = str(Path(os.environ.get('SystemRoot', 'C:/Windows')) / 'System32') if os.name == 'nt' else '/nonexistent'
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    return env


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-dir', type=Path, required=True, help='Unpacked vela directory')
    parser.add_argument('--reference-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--case', nargs='+', choices=('pn2d_0v', 'pn2d_iv', 'pn2d_bv', 'gates'),
                        default=['pn2d_0v', 'pn2d_iv', 'pn2d_bv', 'gates'])
    parser.add_argument('--timeout-seconds', type=float, default=7200)
    args = parser.parse_args(argv)
    root, output = args.package_dir.resolve(), args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = verify_package(root)
    suffix = '.exe' if os.name == 'nt' else ''
    launcher, importer = root / ('bin/vela' + suffix), root / ('bin/sentaurus_import' + suffix)
    python = root / ('bin/python/python.exe' if os.name == 'nt' else 'bin/python/bin/python3')
    if not all(path.is_file() for path in (launcher, importer, python)):
        raise ValueError('package runtime or TDR importer missing')
    script = REPO / 'scripts/sdevice_acceptance.py'
    bootstrap = 'import runpy,sys; sys.path.insert(0,sys.argv.pop(1)); script=sys.argv.pop(1); sys.argv[0]=script; runpy.run_path(script,run_name="__main__")'
    summary = {'package_commit': manifest['commit'], 'platform': manifest['platform'],
               'package_manifest_sha256': sha(root / 'PACKAGE_MANIFEST.json'),
               'controller_sha256': sha(script), 'package_python': str(python),
               'package_launcher': str(launcher), 'host_python_on_path': False, 'fresh_candidates': True, 'cases': []}
    for case in args.case:
        command = [str(python), '-I', '-c', bootstrap, str(REPO / 'scripts'), str(script),
                   '--case', case, '--reference-dir', str(args.reference_dir.resolve()),
                   '--output-dir', str(output / case), '--launcher', str(launcher),
                   '--tdr-importer', str(importer), '--timeout-seconds', str(args.timeout_seconds)]
        with (output / (case + '_controller.log')).open('w', encoding='utf-8') as log:
            done = subprocess.run(command, env=runtime_environment(), stdout=log, stderr=subprocess.STDOUT)
        reports = list((output / case).glob('*/acceptance_report.json'))
        if len(reports) != 1:
            raise ValueError('expected one fresh report for ' + case + '; use a new output directory')
        report = json.loads(reports[0].read_text())
        summary['cases'].append({'case': case, 'exit_code': done.returncode, 'status': report['status'],
                                 'report': str(reports[0]), 'sha256': sha(reports[0])})
        (output / 'package_e2e_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps(summary['cases'][-1]), flush=True)
    verify_package(root)
    summary['package_unchanged'] = True
    summary['status'] = 'pass' if all(item['exit_code'] == 0 and item['status'] == 'pass' for item in summary['cases']) else 'fail'
    (output / 'package_e2e_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return 0 if summary['status'] == 'pass' else 1


if __name__ == '__main__':
    sys.exit(main())
