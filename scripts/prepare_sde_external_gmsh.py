"""Prepare a separately supplied, hash-locked Gmsh runtime for qualification."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile
from package_vela import distribution_platform,safe_member
from prepare_sde_offline_assets import download

REPO=Path(__file__).resolve().parents[1]


def prepare(assets,runtime):
    platform=distribution_platform()
    lock=json.loads((REPO/'configs/sde_external_runtimes.json').read_text())
    archive=download(lock[platform],assets)
    if runtime.exists():raise ValueError('use a fresh external runtime directory')
    if os.name=='nt':
        runtime.mkdir(parents=True)
        with zipfile.ZipFile(archive) as source:
            for item in source.infolist():
                safe_member(item.filename)
                if (item.external_attr>>16)&0o170000==0o120000:raise ValueError('SDK symlink unsupported')
            source.extractall(runtime)
        def unique(pattern):
            paths=list(runtime.rglob(pattern))
            if len(paths)!=1:raise ValueError('ambiguous SDK component: '+pattern)
            return paths[0]
        components=[unique('gmsh.py'),unique('gmsh-4.15.dll'),unique('gmsh.exe')]
    else:
        subprocess.run([sys.executable,'-m','venv','--system-site-packages',str(runtime)],check=True)
        subprocess.run([str(runtime/'bin/python'),'-m','pip','install','--no-index','--no-deps',str(archive)],check=True)
        components=[next(runtime.rglob('gmsh.py')),next(runtime.rglob('libgmsh.so.4.15')),runtime/'bin/gmsh']
    result=dict(zip(('VELA_GMSH_PYTHON','VELA_GMSH_LIBRARY','VELA_GMSH_EXECUTABLE'),map(lambda p:str(p.resolve()),components)))
    (runtime/'external-paths.json').write_text(json.dumps(result,indent=2)+'\n')
    if os.environ.get('GITHUB_ENV'):
        with Path(os.environ['GITHUB_ENV']).open('a',encoding='utf-8') as target:
            for key,value in result.items():
                if '\n' in value or '\r' in value:raise ValueError('unsafe environment path')
                target.write(key+'='+value+'\n')
        with Path(os.environ['GITHUB_PATH']).open('a',encoding='utf-8') as target:
            for directory in sorted({components[1].parent.resolve(),components[2].parent.resolve()}):
                target.write(str(directory)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets',type=Path,required=True)
    parser.add_argument('--runtime',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(prepare(args.assets.resolve(),args.runtime.resolve())))


if __name__=='__main__':main()
