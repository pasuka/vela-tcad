"""Preserve unsupported reads; repeat the same physics with qualified exports."""
import argparse
import shutil
import tarfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import simplemos_phumob_calibration_native_20260908 as original
import prepare_simplemos_runtime_geometry_only_20260907 as known

a,d=original.a,original.d
LOCAL=original.LOCAL/'supported_export'
OUT=original.OUT/'supported_export'
REMOTE=original.REMOTE+'/supported_export'


def prepare():
    a.verify(original.OUT/'native_freeze.json');a.verify(known.OUT/'freeze.json')
    old_jobs=a.read(known.OUT/'contract.json')['jobs']
    tcl=known.LOCAL/'bundle'/old_jobs[0]['key']/'runtime.tcl'
    script=tcl.read_text()
    assert 'Potential' not in script and 'ReadVector' not in script
    files=[Path(__file__).resolve(),Path(original.__file__),original.OUT/'native_freeze.json',known.OUT/'freeze.json',tcl]
    contract=a.read(original.OUT/'native_contract.json')
    for job in contract['jobs']:
        src=original.LOCAL/'bundle'/job['model']/job['key']
        dest=LOCAL/'bundle'/job['model']/job['key'];dest.mkdir(parents=True,exist_ok=False)
        for path in src.iterdir():
            if path.name=='runtime.tcl':continue
            shutil.copyfile(path,dest/path.name);files += [path,dest/path.name]
        (dest/'runtime.tcl').write_text(script,newline='\n');files.append(dest/'runtime.tcl')
    shutil.copyfile(original.LOCAL/'run.sh',LOCAL/'run.sh');files.append(LOCAL/'run.sh')
    contract.update(remote=REMOTE,prior_failure='Unsupported Edge-RegionWise eMobility request in initial export callback; initial results retained separately.',
        export_amendment='Use previously qualified vertex/element mobility and mesh/box geometry reads only. All physical fields come from supported TDR Plot datasets. Native physics, seeds, solver and gates unchanged.')
    a.write(OUT/'native_contract.json',contract);files.append(OUT/'native_contract.json')
    d.matrix.freeze(OUT/'native_freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print('Prepared identical 16 targets with qualified runtime callback',flush=True)


def collect_failure():
    a.verify(original.OUT/'native_freeze.json')
    dest=original.LOCAL/'native_raw';assert not dest.exists()
    with tarfile.open(original.LOCAL/'results.tgz') as tar:
        for m in tar.getmembers():assert (dest/m.name).resolve().is_relative_to(dest.resolve()) and not m.issym() and not m.islnk()
        tar.extractall(dest,filter='data')
    for path in (original.LOCAL/'bundle').rglob('*'):
        if path.is_file():assert a.sha(path)==a.sha(dest/'bundle'/path.relative_to(original.LOCAL/'bundle'))
    rows=[]
    for job in a.read(original.OUT/'native_contract.json')['jobs']:
        p=dest/'bundle'/job['model']/job['key'];log=(p/'console.log').read_text(errors='replace')
        rows.append(dict(**job,exit_code=int((p/'exit_code.txt').read_text()),qualified=False,
                         unsupported_edge_mobility='undefined Edge-RegionWise eMobility' in log))
    a.write_csv(original.OUT/'failed_attempts.csv',rows)
    d.matrix.freeze(original.OUT/'failed_evidence.json',[original.OUT/'native_freeze.json',original.OUT/'failed_attempts.csv',original.LOCAL/'results.tgz']+[p for p in dest.rglob('*') if p.is_file()])
    print('Retained',len(rows),'failed native export attempts',flush=True)


def unpack():
    original.LOCAL,original.OUT=LOCAL,OUT
    original.unpack()


def export():
    a.verify(OUT/'native_evidence.json')
    jobs=[dict(case=j['key'],index=j['index'],tdr=str(LOCAL/'native_raw/bundle'/j['model']/j['key']/'final_des.tdr'),
               export=str(LOCAL/'exports'/j['model']/j['key'])) for j in a.rows(OUT/'native_points.csv') if j['native_qualified']=='True']
    a.write(OUT/'export_contract.json',dict(jobs=jobs))
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(original.b.exporter.export_one,jobs))
    d.matrix.freeze(OUT/'export_evidence.json',[OUT/'native_evidence.json',OUT/'export_contract.json']+[p for p in (LOCAL/'exports').rglob('*') if p.is_file()])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','collect_failure','unpack','export'))
    globals()[parser.parse_args().action]()
