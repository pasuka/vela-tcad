"""Package the current checkout and frozen inputs for an isolated T470p build."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_engineering_20260929'


def main():
    package=OUT/'package_v2';package.mkdir(exist_ok=False)
    files=subprocess.check_output(['git','-c','core.fsmonitor=false','ls-files','--cached','--others','--exclude-standard'],cwd=ROOT,text=True).splitlines()
    for rel in sorted(set(files)):
        src=ROOT/rel
        if src.is_file():
            dst=package/'source'/rel
            if os.name=='nt':
                src=Path('\\\\?\\'+str(src.resolve()))
                dst=Path('\\\\?\\'+str(dst.resolve()))
            dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    shutil.copytree(ROOT/'build/hfs_srh_20260926/original_inputs',package/'data/inputs')
    seed_dir=package/'data/seeds';seed_dir.mkdir(parents=True)
    for src in (ROOT/'build/outlier_analysis_20260928/detailed').glob('*.h5'):
        shutil.copy2(src,seed_dir/src.name)
    for src,name in [(ROOT/'build/outlier_analysis_20260928/all_points.csv','frozen_points.csv'),
                     (ROOT/'build/outlier_m60_review_20260928/joined_816.csv','m60_joined.csv'),
                     (OUT/'engineering_contract.json','engineering_contract.json')]:
        shutil.copy2(src,package/'data'/name)
    job=r'''$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;'+$env:Path
$env:VELA_LINEAR_SOLVER='sparselu'
$env:OMP_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:VELA_LINEAR_THREADS='1'
$base=$PSScriptRoot
Set-Location -LiteralPath "$base/source"
try {
  & python -B -X utf8 "$base/verify_package.py"
  if ($LASTEXITCODE -ne 0) { throw 'Package identity check failed' }
  & cmake --preset windows-ucrt64-release *> "$base/configure.log"
  if ($LASTEXITCODE -ne 0) { throw 'Configure failed' }
  & cmake --build --preset windows-ucrt64-release --parallel 2 *> "$base/build.log"
  if ($LASTEXITCODE -ne 0) { throw 'Build failed' }
  $test=Start-Process D:/msys64/ucrt64/bin/python.exe -ArgumentList '-B','-X','utf8','tests/regression/test_simplemos_engineering_audit.py' -WindowStyle Hidden -Wait -PassThru -RedirectStandardOutput "$base/audit_tests.stdout.log" -RedirectStandardError "$base/audit_tests.stderr.log"
  if ($test.ExitCode -ne 0) { throw 'Audit harness tests failed' }
  '0' | Set-Content "$base/build.exit"
} catch {
  $_ | Out-String | Set-Content "$base/build.failure.txt"
  '1' | Set-Content "$base/build.exit"
  exit 1
}
'''
    (package/'build_job.ps1').write_text(job,encoding='utf-8')
    verifier='''import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent
manifest=json.loads((root/'package_hashes.json').read_text())
for rel,expected in manifest.items():
    actual=hashlib.sha256((root/rel).read_bytes()).hexdigest()
    if actual!=expected:raise ValueError(rel)
(root/'package_verified.json').write_text(json.dumps(dict(verified_files=len(manifest))))
print('Package files verified:',len(manifest))
'''
    (package/'verify_package.py').write_text(verifier,encoding='utf-8')
    hashes={p.relative_to(package).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(package.rglob('*')) if p.is_file()}
    (package/'package_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf-8')
    archive=OUT/'t470p_package.tgz'
    with tarfile.open(archive,'w:gz') as tar:
        for p in sorted(package.rglob('*')):
            if p.is_file():tar.add(p,arcname=p.relative_to(package).as_posix())
    seal=dict(sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),bytes=archive.stat().st_size,files=len(hashes))
    (OUT/'transfer_seal.json').write_text(json.dumps(seal,indent=2)+'\n')
    print(json.dumps(seal))


if __name__=='__main__':main()
