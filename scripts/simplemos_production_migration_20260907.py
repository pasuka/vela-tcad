"""Production migration provenance and immutable pre-change evidence access."""
from pathlib import Path
import shutil
import subprocess
import analyze_simplemos_production_consistency_20260907 as previous

a=previous.a; d=previous.d; p=previous.p; b=previous.run.b
LOCAL=p.REPO/'build-release/simplemos_production_migration_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/production_migration_20260907'
RUNNER=p.REPO/'build-release/vela_example_runner.exe'


def snapshot():
    a.verify(previous.OUT/'validation_evidence.json'); a.verify(previous.prior.OUT/'validation_evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=False)
    names=subprocess.check_output(['git','ls-files','src','include','tests','CMakeLists.txt','docs/config_schema.md'],cwd=p.REPO,text=True).splitlines()
    names += ['build-release/libvela_core.a','build-release/vela_example_runner.exe']
    records={}
    for name in names:
        src=p.REPO/name
        if not src.is_file():continue
        dest=LOCAL/'prechange'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(src,dest)
        records[name]=dict(snapshot=a.rel(dest),sha256=a.sha(src))
    a.write(OUT/'prechange.json',dict(status='frozen_before_production_edits',files=records,
        original_tracked_diff='127 insertions / 1 deletion in src/tools/vela_example_runner.cpp',
        previous_evidence=a.rel(previous.OUT/'validation_evidence.json')))
    d.matrix.freeze(OUT/'prechange_evidence.json',[OUT/'prechange.json']+[p.REPO/r['snapshot'] for r in records.values()])
    print('Preserved',len(records),'source/test preimages and original Release library/runner.',flush=True)


def verify_historical(manifest):
    records=a.read(OUT/'prechange.json')['files']
    for name,expected in a.read(manifest)['input_hashes'].items():
        path=p.REPO/name
        if path.is_file() and a.sha(path)==expected:continue
        key=Path(name).as_posix()
        assert key in records and records[key]['sha256']==expected,(name,'unexpected historical drift')
        assert a.sha(p.REPO/records[key]['snapshot'])==expected


if __name__=='__main__':snapshot()
