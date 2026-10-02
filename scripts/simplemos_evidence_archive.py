"""Verify or restore the indexed SimpleMOS historical evidence (never solves).

Only manifest-listed paths can be restored. Existing differing files and paths
escaping the destination, including through symlinks, are rejected.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'reference_tcad/simplemos_sentaurus2022/'
MANIFEST = ROOT / PREFIX / 'local_evidence_20260930/manifest.json'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def safe_path(root, name):
    path = PurePosixPath(name)
    if (not name.startswith(PREFIX) or '\\' in name or ':' in name
            or path.is_absolute() or '..' in path.parts
            or path.as_posix() != name):
        raise ValueError('Unsafe evidence path: ' + name)
    root = root.resolve()
    target = (root / name).resolve()
    if not target.is_relative_to(root):
        raise ValueError('Evidence path escapes destination: ' + name)
    return target


def load_manifest(path):
    data = json.loads(path.read_text(encoding='utf-8-sig'))
    if data['schema'] != 'vela.simplemos.local-evidence.v1':
        raise ValueError('Unsupported evidence manifest')
    seen = set()
    for item in data['files']:
        name = item['path']
        safe_path(Path.cwd(), name)
        if name.casefold() in seen:
            raise ValueError('Duplicate evidence path: ' + name)
        seen.add(name.casefold())
        if item['storage'] not in ('git', 'archive') or item['bytes'] < 0:
            raise ValueError('Invalid evidence entry: ' + name)
        digest = item['sha256']
        if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Invalid evidence hash: ' + name)
    if len(seen) != data['file_count']:
        raise ValueError('Manifest count mismatch')
    return data


def verify_files(data, root):
    retained = 0
    cached = 0
    for item in data['files']:
        target = safe_path(root, item['path'])
        if item['storage'] == 'archive' and not target.exists():
            continue
        if (not target.is_file() or target.stat().st_size != item['bytes']
                or sha(target) != item['sha256']):
            raise ValueError('Missing or changed evidence: ' + item['path'])
        if item['storage'] == 'git':
            retained += 1
        else:
            cached += 1
    return {'retained_verified': retained, 'local_archive_copies_verified': cached}


def verify_archive(data, archive):
    if archive.stat().st_size != data['archive']['bytes'] or sha(archive) != data['archive']['sha256']:
        raise ValueError('Archive byte identity mismatch')
    entries = {x['path']: x for x in data['files']}
    with zipfile.ZipFile(archive) as source:
        names = source.namelist()
        if len(names) != len(set(names)) or set(names) != set(entries):
            raise ValueError('Archive member coverage mismatch')
        for info in source.infolist():
            item = entries[info.filename]
            if stat.S_ISLNK(info.external_attr >> 16) or info.file_size != item['bytes']:
                raise ValueError('Invalid archive member: ' + info.filename)
            digest = hashlib.sha256()
            with source.open(info) as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                    digest.update(chunk)
            if digest.hexdigest() != item['sha256']:
                raise ValueError('Archive member hash mismatch: ' + info.filename)
    return len(entries)


def restore(data, archive, destination, selected=None):
    # Validate every member before any filesystem write.
    verified = verify_archive(data, archive)
    entries = {x['path']: x for x in data['files']}
    names = sorted(selected if selected is not None else
                   [x['path'] for x in data['files'] if x['storage'] == 'archive'])
    if len(names) != len(set(names)) or any(name not in entries for name in names):
        raise ValueError('Unknown or duplicate restore selection')
    targets = {}
    for name in names:
        target = safe_path(destination, name)
        if target.exists() and (not target.is_file() or sha(target) != entries[name]['sha256']):
            raise ValueError('Refusing to overwrite differing evidence: ' + name)
        targets[name] = target
    created = 0
    with zipfile.ZipFile(archive) as source:
        for name, target in targets.items():
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            # Recheck containment after creating parents; never overwrite.
            safe_path(destination, name)
            with source.open(name) as inp, target.open('xb') as out:
                for chunk in iter(lambda: inp.read(1024 * 1024), b''):
                    out.write(chunk)
            if sha(target) != entries[name]['sha256']:
                raise ValueError('Restored evidence failed hash check: ' + name)
            created += 1
    return {'archive_members_verified': verified, 'created': created,
            'existing_identical': len(names) - created}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('verify', 'restore'))
    parser.add_argument('--manifest', type=Path, default=MANIFEST)
    parser.add_argument('--root', type=Path, default=ROOT,
                        help='Checkout to verify or destination to restore')
    parser.add_argument('--archive', type=Path,
                        help='Local archive; required for restore, optional for verify')
    parser.add_argument('--path', action='append', dest='paths',
                        help='Exact repo-relative manifest member to restore; repeatable')
    args = parser.parse_args()
    data = load_manifest(args.manifest)
    if args.action == 'restore':
        if args.archive is None:
            parser.error('restore requires --archive; no remote download is implicit')
        result = restore(data, args.archive, args.root, args.paths)
    else:
        if args.paths:
            parser.error('--path is only valid for restore')
        result = verify_files(data, args.root)
        result['archive_verified'] = args.archive is not None
        if args.archive is not None:
            result['archive_members_verified'] = verify_archive(data, args.archive)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
