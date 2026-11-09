"""Download the locked build inputs; the packager itself is strictly offline."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile
import urllib.request
from package_vela import distribution_platform

REPO=Path(__file__).resolve().parents[1]


def download(record,directory):
    directory.mkdir(parents=True,exist_ok=True)
    name=record['filename']
    if Path(name).name!=name:raise ValueError('invalid input basename')
    target=directory/name
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('existing offline input hash mismatch: '+name)
        return target
    temporary=target.with_name(target.name+'.download')
    if temporary.exists():raise ValueError('preexisting download file: '+temporary.name)
    try:
        with urllib.request.urlopen(record['url'],timeout=120) as source,temporary.open('xb') as out:
            shutil.copyfileobj(source,out)
        if hashlib.sha256(temporary.read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('download hash mismatch: '+name)
        temporary.replace(target)
    finally:
        if temporary.exists():temporary.unlink()
    return target


def prepare(lock,assets,platform):
    key='python_win' if platform=='win64' else 'python_linux'
    download(lock['python'][key],assets)
    for item in lock['dependencies'][platform]:download(item,assets/platform)
    archive=download(lock['python']['python_licenses'],assets)
    licenses=assets/'python-licenses';licenses.mkdir(exist_ok=True)
    with tarfile.open(archive) as source:
        for member in source.getmembers():
            if member.isfile() and len(Path(member.name).parts)==2 and Path(member.name).name.startswith('LICENSE'):
                data=source.extractfile(member).read()
                path=licenses/Path(member.name).name
                if path.exists() and path.read_bytes()!=data:raise ValueError('license text collision')
                path.write_bytes(data)
    if not (licenses/'LICENSE.openssl-3.txt').is_file():raise ValueError('incomplete upstream license source')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assets',type=Path,required=True)
    parser.add_argument('--lock',type=Path,default=REPO/'configs/sde_offline_inputs.json')
    parser.add_argument('--platform',choices=('win64','ubuntu26.04-x86_64'))
    args=parser.parse_args()
    prepare(json.loads(args.lock.read_text()),args.assets,args.platform or distribution_platform())
    print('Locked offline inputs verified')


if __name__=='__main__':main()
