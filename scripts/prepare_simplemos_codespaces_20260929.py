"""Package only the unstarted engineering cases; preserve T470p ownership."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tarfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_engineering_20260929/cloud_submission'
FROZEN=ROOT/'build/simplemos_engineering_20260929/package_v2'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')


def main(out):
    package=out/'package';package.mkdir(parents=True,exist_ok=False)
    shutil.copytree(FROZEN/'data',package/'data')
    for folder,pattern in [('scripts','*.py'),('tests/regression','test_simplemos_engineering_*.py')]:
        target=package/'source'/folder;target.mkdir(parents=True,exist_ok=True)
        for path in (ROOT/folder).glob(pattern):shutil.copy2(path,target/path.name)
    # Full core source comparison, with line-ending normalization reported
    # separately. No inference from a config hash to binary identity.
    core={}
    for folder in ('include','src'):
        for p in sorted((FROZEN/'source'/folder).rglob('*')):
            if p.is_file():
                rel=p.relative_to(FROZEN/'source').as_posix();data=p.read_bytes()
                if (ROOT/rel).read_bytes()!=data:raise ValueError('Local core changed since T470p snapshot: '+rel)
                core[rel]=dict(raw=sha(p),lf=hashlib.sha256(data.replace(b'\r\n',b'\n')).hexdigest())
    write(package/'core_source_identity.json',core)
    cases=json.loads((package/'data/inputs/contract.json').read_text())['cases']
    excluded={'n17_vd_0p05':'completed_on_T470p','n17_vd_1':'failed_bad_alloc_retained_not_retried','n18_vd_0p05':'in_progress_on_T470p'}
    pending=[c['case'] for c in cases if c['case'] not in excluded]
    if len(pending)!=13:raise ValueError('Unexpected unstarted set')
    write(package/'assignment.json',dict(schema='vela.simplemos.assignment.v1',codespace='vela-tcad-compute-69r6rj7pvvjc54g9',
        cloud_base='/workspaces/simplemos-engineering-20260929',cases=pending,points=663,
        retained_T470p=excluded,max_workers=2,reference='M60 default',current_gate_percent=2,
        solver_numerical_changes=False,complete_816_requires_failed_T470p_case=True,
        frozen_source='/workspaces/simplemos-hfs-merged-20260926/source-continuation-v2',
        runner_sha256='87b934c79451110467cee69b62f6ec2e312a8aea51ab85c07bf2c0577ad05ea3',
        importer_sha256='f59df0314ce51f0604e6ece520ab02796b79c996c799b7048b6e56dc4dae7b4f'))
    native=ROOT/'build/simplemos_engineering_20260929/native_fields'
    write(package/'native_archives.json',{f'n{i}.tgz':sha(native/f'n{i}.tgz') for i in range(18,25)})
    manifest={p.relative_to(package).as_posix():sha(p) for p in sorted(package.rglob('*')) if p.is_file()}
    write(package/'package_hashes.json',manifest)
    with tarfile.open(out/'submission.tgz','w:gz') as tar:
        for p in sorted(package.rglob('*')):
            if p.is_file():tar.add(p,arcname=p.relative_to(package).as_posix())
    seal=dict(sha256=sha(out/'submission.tgz'),files=len(manifest),bytes=(out/'submission.tgz').stat().st_size)
    write(out/'transfer_seal.json',seal);print(json.dumps(seal))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUT)
    main(p.parse_args().output)
