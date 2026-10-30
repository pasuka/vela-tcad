"""Assemble an offline SDE tools package; no Gmsh binary or solver is bundled."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
import zipfile

from package_vela import (distribution_platform, unpack_python, safe_member,
                          collect_linux_runtime_dependencies, normalize_linux_runpaths)

REPO = Path(__file__).resolve().parents[1]
LOCK = REPO/'configs/sde_offline_inputs.json'


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def module_closure(repo=REPO):
    todo,seen = ['sentaurus_gmsh_pipeline'],set()
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        seen.add(name)
        tree = ast.parse((repo/'scripts'/(name+'.py')).read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            imports = ([a.name for a in node.names] if isinstance(node,ast.Import) else
                       [node.module or ''] if isinstance(node,ast.ImportFrom) else [])
            todo.extend(n.split('.')[0] for n in imports if
                        (repo/'scripts'/(n.split('.')[0]+'.py')).is_file())
    return sorted(seen)


def verify_inputs(lock, assets, platform):
    python = lock['python']['python_win' if platform=='win64' else 'python_linux']
    records = [python,*lock['dependencies'][platform]]
    paths = []
    for record in records:
        safe_member(record['filename'])
        path = assets/record['filename'] if record is python else assets/platform/record['filename']
        if not path.is_file() or sha256(path)!=record['sha256']:
            raise ValueError('offline input missing or hash mismatch: '+record['filename'])
        paths.append(path)
    return paths


def install_wheels(paths, vendor):
    vendor.mkdir(parents=True)
    for path in paths:
        with zipfile.ZipFile(path) as archive:
            for item in archive.infolist():
                relative = safe_member(item.filename)
                if (item.external_attr>>16)&0o170000 == 0o120000:
                    raise ValueError('wheel symlink is unsupported')
                if '.data' in relative.parts:
                    raise ValueError('wheel requires unsupported installation scheme')
                target = vendor.joinpath(*relative.parts)
                if item.is_dir():
                    target.mkdir(parents=True,exist_ok=True)
                else:
                    if target.exists():
                        raise ValueError('wheel file collision: '+item.filename)
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with archive.open(item) as source,target.open('wb') as output:
                        shutil.copyfileobj(source,output)


def license_texts(stage, assets, dependencies):
    source = assets/'python-licenses'
    if not (source/'LICENSE.openssl-3.txt').is_file():
        raise ValueError('supply the matching upstream standalone Python license texts')
    shutil.copytree(source,stage/'licenses/python')
    for library in dependencies:
        # Preserve the exact platform copyright records for copied ELF libraries.
        result = subprocess.check_output(['dpkg-query','-S',str(library.resolve())],text=True)
        package = result.split(':',1)[0]
        copyright = Path('/usr/share/doc')/package/'copyright'
        destination = stage/'licenses/native'/package
        destination.mkdir(parents=True,exist_ok=True)
        shutil.copy2(copyright,destination/'copyright')
    if dependencies:
        shutil.copytree('/usr/share/common-licenses',stage/'licenses/common')
    licenses = list((stage/'vendor').glob('*.dist-info/licenses/*'))
    # Some wheels use the pre-PEP639 dist-info/LICENSE layout.
    licenses += [p for d in (stage/'vendor').glob('*.dist-info') for p in d.iterdir()
                 if p.is_file() and ('LICENSE' in p.name.upper() or 'COPYING' in p.name.upper())]
    if len(licenses)<6:
        raise ValueError('dependency wheels lack required upstream license texts')
    (stage/'THIRD_PARTY_NOTICES').write_text(
        'CPython and bundled libraries: licenses/python and bin/python.\n'
        'Python dependencies: vendor/*.dist-info/licenses or LICENSE files.\n'
        'Copied ELF libraries: licenses/native and licenses/common.\n'
        'External Gmsh 4.15.2 is NOT included in this archive.\n'
        'Gmsh license and upstream exception: https://gmsh.info/doc/texinfo/#License\n'
        'Configure and provide that runtime separately; this package does not change its license.\n',encoding='utf-8')


def package(args):
    platform = distribution_platform()
    lock = json.loads(args.lock.read_text(encoding='utf-8'))
    paths = verify_inputs(lock,args.assets,platform)
    for name,digest in lock['python']['python_licenses']['files'].items():
        safe_member(name)
        if sha256(args.assets/'python-licenses'/name)!=digest:
            raise ValueError('upstream Python license hash mismatch: '+name)
    archive = args.output.resolve()
    if archive.exists():
        raise ValueError('refusing to overwrite existing package: '+archive.name)
    archive.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='vela-sde-',dir=archive.parent if os.name=='nt' else None) as temp:
        stage = Path(temp)/'vela-sde';stage.mkdir()
        unpack_python(paths[0],stage/'bin/python',platform)
        install_wheels(paths[1:],stage/'vendor')
        modules = module_closure()
        (stage/'scripts').mkdir()
        for name in modules:
            shutil.copy2(REPO/'scripts'/(name+'.py'),stage/'scripts'/(name+'.py'))
        shutil.copytree(REPO/'schemas',stage/'schemas')
        shutil.copy2(REPO/'scripts/sde_package_entry.py',stage/'run_sde.py')
        shutil.copy2(REPO/'LICENSE',stage/'LICENSE')
        (stage/'lib').mkdir(exist_ok=True)
        dependencies = collect_linux_runtime_dependencies(stage) if os.name!='nt' else []
        runpaths = normalize_linux_runpaths(stage,args.patchelf) if os.name!='nt' else []
        license_texts(stage,args.assets,dependencies)
        interpreter = stage/('bin/python/python.exe' if os.name=='nt' else 'bin/python/bin/python3')
        result = subprocess.run([str(interpreter),'-B','-I','-S','-c',
            'import sys;sys.path.insert(0,sys.argv[1]);from jsonschema import Draft7Validator;'
            'Draft7Validator({"type":"integer"}).validate(1);import ctypes;print(sys.version.split()[0])',str(stage/'vendor')],
            capture_output=True,text=True,check=True)
        revision = subprocess.run(['git','-c','safe.directory='+REPO.as_posix(),'-c','core.fsmonitor=false','rev-parse','HEAD'],cwd=REPO,capture_output=True,text=True)
        dirty=subprocess.run(['git','-c','safe.directory='+REPO.as_posix(),'-c','core.fsmonitor=false','status','--porcelain'],cwd=REPO,capture_output=True,text=True)
        source_commit=revision.stdout.strip() if revision.returncode==0 else args.source_commit
        source_dirty=bool(dirty.stdout.strip()) if dirty.returncode==0 else (
            args.source_tree_state=='dirty' if args.source_tree_state else None)
        import re
        if not source_commit or not re.fullmatch(r'[0-9a-f]{40}',source_commit) or source_dirty is None:
            raise ValueError('source snapshot requires --source-commit and --source-tree-state')
        if args.source_commit and args.source_commit!=source_commit:
            raise ValueError('explicit source commit differs from Git checkout')
        if args.source_tree_state and (args.source_tree_state=='dirty')!=source_dirty:
            raise ValueError('explicit source tree state differs from Git checkout')
        files = []
        for path in sorted(stage.rglob('*')):
            if path.is_file():
                relative = path.relative_to(stage)
                if path.suffix.lower() in {'.tdr','.log','.grd','.dat','.msh','.cmd'} or 'reference_tcad' in relative.parts:
                    raise ValueError('simulation artifact in package: '+str(relative))
                files.append({'path':relative.as_posix(),'sha256':sha256(path)})
        manifest = {'schema':'vela.sde_offline_package.v1','platform':platform,
                    'source_commit':source_commit,
                    'source_tree_dirty':source_dirty,'python_version':result.stdout.strip(),'gmsh_profile':'explicit external 4.15.2 (not bundled)',
                    'source_modules':modules,'inputs':{'python_sha256':lock['python']['python_win' if os.name=='nt' else 'python_linux']['sha256'],
                    'dependency_wheels':lock['dependencies'][platform]},'elf_runpaths':runpaths,'files':files}
        (stage/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        temporary = archive.with_name(archive.name+'.tmp')
        try:
            if os.name=='nt':
                with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED) as output:
                    for path in sorted(stage.rglob('*')):
                        if path.is_file():output.write(path,'vela-sde/'+path.relative_to(stage).as_posix())
            else:
                with tarfile.open(temporary,'w:gz') as output:output.add(stage,arcname='vela-sde')
            temporary.replace(archive)
        finally:
            if temporary.exists():temporary.unlink()
        archive.with_name(archive.name+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
        archive.with_name(archive.name+'.sha256').write_text(sha256(archive)+'  '+archive.name+'\n',encoding='ascii')
    return archive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets',type=Path,required=True,help='Offline wheel/Python/license directory')
    parser.add_argument('--lock',type=Path,default=LOCK)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--patchelf',default='patchelf')
    parser.add_argument('--source-commit',help='Required with tree state for a non-Git source snapshot')
    parser.add_argument('--source-tree-state',choices=['clean','dirty'])
    args = parser.parse_args()
    try:print(package(args))
    except (ValueError,OSError,subprocess.CalledProcessError) as error:parser.exit(2,'package_sde: '+str(error)+'\n')


if __name__=='__main__':main()
