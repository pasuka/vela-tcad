"""Pack/verify/restore byte-pinned raw assets. Never invokes a solver.

Verification and restoration also support the Sentaurus VM's Python 3.6.
The inventory is an explicit list; packing never follows an implicit file glob.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import zipfile


def digest(stream):
    value = hashlib.sha256()
    for block in iter(lambda: stream.read(1048576), b''):
        value.update(block)
    return value.hexdigest()


def safe(root, name):
    p = PurePosixPath(name)
    if (not name or p.is_absolute() or '..' in p.parts or '.git' in p.parts
            or '\\' in name or ':' in name or p.as_posix() != name):
        raise ValueError('Unsafe archive member: ' + name)
    root = root.resolve()
    if os.name == 'nt' and not str(root).startswith('\\\\?\\'):
        root = Path('\\\\?\\' + str(root))
    base = str(root)
    target = (root / name).resolve()
    if os.path.commonpath([base, str(target)]) != base:
        raise ValueError('Archive member escapes destination: ' + name)
    return target


def pack(root, inventory, archive):
    entries = []
    names = [x['path'] for x in inventory['files']]
    if len(set(n.casefold() for n in names)) != len(names):
        raise ValueError('Duplicate inventory path')
    with zipfile.ZipFile(str(archive), 'x', zipfile.ZIP_DEFLATED, allowZip64=True) as z:
        for number, item in enumerate(inventory['files']):
            source = safe(root, item['path'])
            before = source.stat()
            if before.st_size != item['bytes'] or before.st_mtime_ns != item['mtime_ns']:
                raise ValueError('Asset changed since inventory: ' + item['path'])
            h = hashlib.sha256()
            with source.open('rb') as inp, z.open('files/' + item['path'], 'w', force_zip64=True) as out:
                for block in iter(lambda: inp.read(1048576), b''):
                    h.update(block)
                    out.write(block)
            after = source.stat()
            if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                raise ValueError('Asset changed during packing: ' + item['path'])
            entries.append(dict(item, sha256=h.hexdigest()))
            if number % 1000 == 0:
                print('packed {}/{}'.format(number + 1, len(names)), flush=True)
        manifest = dict(inventory, schema='vela.simplemos.raw-archive.v1', files=entries)
        z.writestr('manifest.json', json.dumps(manifest, ensure_ascii=False, indent=2).encode('utf-8'))
    return manifest


def verify(archive, expected=None, destination=None):
    with archive.open('rb') as stream:
        archive_hash = digest(stream)
    if expected and archive_hash != expected:
        raise ValueError('Archive SHA256 mismatch')
    with zipfile.ZipFile(str(archive)) as z:
        manifest_bytes = z.read('manifest.json')
        data = json.loads(manifest_bytes.decode('utf-8'))
        if data['schema'] != 'vela.simplemos.raw-archive.v1':
            raise ValueError('Unknown schema')
        entries = data['files']
        names = [x['path'] for x in entries]
        if len(set(n.casefold() for n in names)) != len(names):
            raise ValueError('Duplicate manifest path')
        expected_names = {'manifest.json'} | {'files/' + n for n in names}
        if len(z.namelist()) != len(expected_names) or set(z.namelist()) != expected_names:
            raise ValueError('Archive member coverage mismatch')
        for item in entries:
            safe(destination or Path.cwd(), item['path'])
            info = z.getinfo('files/' + item['path'])
            if info.file_size != item['bytes']:
                raise ValueError('Member size mismatch: ' + item['path'])
            with z.open(info) as stream:
                if digest(stream) != item['sha256']:
                    raise ValueError('Member SHA256 mismatch: ' + item['path'])
        created = 0
        if destination is not None:
            # Complete verification and conflict preflight precede all writes.
            for item in entries:
                target = safe(destination, item['path'])
                if target.exists():
                    with target.open('rb') as stream:
                        if digest(stream) != item['sha256']:
                            raise ValueError('Refusing to overwrite: ' + item['path'])
            for item in entries:
                target = safe(destination, item['path'])
                if target.exists():
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                safe(destination, item['path'])
                with z.open('files/' + item['path']) as inp, target.open('xb') as out:
                    shutil.copyfileobj(inp, out, 1048576)
                with target.open('rb') as stream:
                    if digest(stream) != item['sha256']:
                        raise ValueError('Restored file SHA256 mismatch')
                created += 1
    return {'passed': True, 'archive_sha256': archive_hash,
            'archive_bytes': archive.stat().st_size, 'files_verified': len(entries),
            'source_bytes': sum(x['bytes'] for x in entries),
            'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(),
            'restored': created}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['pack', 'verify', 'restore'])
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--inventory', type=Path)
    p.add_argument('--root', type=Path)
    p.add_argument('--sha256')
    p.add_argument('--receipt', type=Path)
    a = p.parse_args()
    if a.action == 'pack':
        if a.inventory is None or a.root is None:
            p.error('pack requires --inventory and --root')
        pack(a.root, json.loads(a.inventory.read_text(encoding='utf-8')), a.archive)
    elif a.action == 'restore' and (a.root is None or not a.sha256):
        p.error('restore requires --root and --sha256')
    result = verify(a.archive, a.sha256, a.root if a.action == 'restore' else None)
    if a.receipt:
        a.receipt.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
