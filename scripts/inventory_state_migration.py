"""Inventory active source/config state references without scanning historical output.

Administrative inventory only: no solver execution, conversion or file removal.
The input snapshot includes uncommitted files; Git HEAD alone is not its identity.
"""
import argparse
import hashlib
import json
import subprocess
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AREAS = ('src/', 'include/', 'scripts/', 'tests/', 'configs/', 'examples/', 'reference_tcad/',
         '.github/', 'schemas/', 'cmake/')
STATE_KEYS = {'initial_state_file', 'write_state_file', 'output_state_file',
              'previous_state_file', 'state_file', 'seed_file', 'initial_secant_state_file'}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inventory(root):
    names = subprocess.check_output(
        ['git', '-c', 'core.fsmonitor=false', 'ls-files', '-z', '--cached',
         '--others', '--exclude-standard'], cwd=root).decode().split('\0')
    files = sorted({p for p in names if p and (p.startswith(AREAS) or
                   '/' not in p) and (root/p).is_file()})
    references, inline, seeds, dependencies = [], [], [], []
    def walk(value, filename, location=''):
        if isinstance(value, dict):
            for key, item in value.items():
                pointer = location+'/'+key
                if key in ('files', 'files_sha256') and isinstance(item, dict):
                    for name, expected in item.items():
                        if isinstance(expected, str) and re.fullmatch('[0-9a-f]{64}', expected):
                            target = root/name
                            dependencies.append(dict(manifest=filename, path=name, sha256=expected,
                                exists_from_root=target.is_file(),
                                status='requires_manifest_base_resolution' if not target.is_file() else 'present'))
                if key == 'seeds' and isinstance(item, dict):
                    for label, name in item.items():
                        if isinstance(name, str):
                            references.append(dict(config=filename, key=pointer+'/'+label, value=name,
                                resolution='existing_candidate' if (root/name).is_file() else 'runtime_or_missing',
                                existing=[name] if (root/name).is_file() else []))
                if key in STATE_KEYS and isinstance(item, str) and item:
                    candidates = [(root/filename).parent/item, root/item]
                    present = [p for p in candidates if p.is_file()]
                    references.append(dict(config=filename, key=pointer, value=item,
                        resolution='existing_candidate' if present else 'runtime_or_missing',
                        existing=[str(p.resolve().relative_to(root.resolve())).replace('\\', '/')
                                  if p.resolve().is_relative_to(root.resolve()) else str(p.resolve())
                                  for p in present]))
                if key in ('state_interleaved', 'referenced_state_interleaved'):
                    inline.append(dict(config=filename, key=pointer, count=len(item)))
                walk(item, filename, pointer)
        elif isinstance(value, list):
            for i, item in enumerate(value):
                if isinstance(item, (dict, list)):
                    walk(item, filename, location+'/'+str(i))
    for filename in files:
        path = root/filename
        if path.suffix == '.json':
            try:
                walk(json.loads(path.read_text(encoding='utf-8-sig')), filename)
            except (ValueError, UnicodeError):
                pass  # Invalid-input regression fixtures are intentionally not configs.
        if path.suffix == '.csv':
            with path.open(encoding='utf-8-sig', errors='replace') as stream:
                header = stream.readline().strip()
            if header.startswith('node_id,psi,phin,phip,electrons_m3,holes_m3'):
                seeds.append(dict(path=filename, bytes=path.stat().st_size, sha256=digest(path)))
    return dict(schema='vela.state_migration.inventory/1',
                git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root).decode().strip(),
                source_sha256={p:digest(root/p) for p in files},
                state_references=references, inline_electrothermal_states=inline,
                csv_state_files=seeds, manifest_dependencies=dependencies,
                exclusions=['build', 'build-release', 'reference_staging', 'historical docs'],
                note='Unresolved paths may be generated outputs; this is not a missing-input verdict.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = inventory(ROOT)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(data, stream, indent=2)
    print(json.dumps({key:len(data[key]) for key in
                     ('source_sha256', 'state_references', 'inline_electrothermal_states', 'csv_state_files')}))
