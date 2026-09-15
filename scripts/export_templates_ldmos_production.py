"""Export frozen electrothermal inputs into a new local run directory (R7 default).

Only declared file paths change. This prepares decks; it never starts a solver
or transfers the reference qualification to a different binary/environment.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil

ROOT = Path(__file__).resolve().parents[1]
PROFILE = ROOT / 'reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r7.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def child(root, relative):
    path = PurePosixPath(relative)
    if path.is_absolute() or '..' in path.parts or '\\' in relative or ':' in relative:
        raise ValueError('Expected a relative evidence path without traversal')
    result = root.joinpath(*path.parts).resolve()
    if not result.is_relative_to(root.resolve()):
        raise ValueError('Evidence path escapes its root')
    return result


def replace_paths(value, mapping):
    if isinstance(value, dict):
        return {key: replace_paths(item, mapping) for key, item in value.items()}
    if isinstance(value, list):
        return [replace_paths(item, mapping) for item in value]
    return mapping.get(value, value) if isinstance(value, str) else value


def export_bundle(evidence, output, profile_path=PROFILE):
    evidence, output = evidence.resolve(), output.resolve()
    if output.exists():
        raise ValueError('Output already exists; choose a new directory')
    profile = read(profile_path)
    sources = {}
    for relative, digest in profile['files_sha256'].items():
        path = child(evidence, relative)
        if sha(path) != digest:
            raise ValueError('Evidence SHA256 mismatch: ' + relative)
        sources[relative] = path
    # Validate and construct everything before writing an output directory.
    mapping = {}
    copies = {}
    for binding in profile['dependencies']:
        source = sources[binding['evidence_path']]
        target = child(output, binding['bundle_path'])
        mapping[binding['recorded_path']] = str(target)
        copies[target] = source
    generated = {}
    for case in profile['cases']:
        cfg = replace_paths(read(sources[case['input']]), mapping)
        inp = output / ('input_vg%d.json' % case['gate_V'])
        deck = read(sources[case['deck']])
        deck['input_file'] = inp.name
        deck['output_directory'] = 'results_vg%d' % case['gate_V']
        generated[inp] = cfg
        generated[output / ('vg%d.json' % case['gate_V'])] = deck
    output.mkdir(parents=True)
    for target, source in copies.items():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    for path, value in generated.items():
        path.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    manifest = dict(
        schema='vela.ldmos.production_export.v1',
        scope='Frozen input/configuration reproduction only; execution and acceptance are separate',
        evidence_root=str(evidence), output_root=str(output),
        profile_sha256=sha(profile_path),
        qualified_linux_runner_sha256=profile['qualified_linux_runner_sha256'],
        sources_sha256=profile['files_sha256'], path_replacements=mapping,
        files_sha256={str(p.relative_to(output)): sha(p) for p in sorted([*copies, *generated])})
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-root', type=Path, required=True,
                        help='Extracted R7 archive root containing cases_r7 and cases/data')
    parser.add_argument('--output', type=Path, required=True,
                        help='New directory; regenerate after moving it because dependency paths are absolute')
    parser.add_argument('--profile', type=Path, default=PROFILE,
                        help='Explicit frozen profile; defaults to qualified R7')
    args = parser.parse_args()
    export_bundle(args.evidence_root, args.output, args.profile)
    print('Prepared Vg4/Vg8 decks: ' + str(args.output.resolve()))


if __name__ == '__main__':
    main()
