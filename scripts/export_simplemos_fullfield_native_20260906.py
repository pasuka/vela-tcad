"""Verify downloaded inputs, preserve raw native output, and export target TDRs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import subprocess
import tarfile
import prepare_simplemos_fullfield_native_20260906 as p
from run_simplemos_m81_output_amendment import rows as pltrows

a=p.a;LOCAL=p.LOCAL;OUT=p.OUT


def unpack():
    a.verify(OUT/'native_freeze.json')
    root=(LOCAL/'native_raw').resolve()
    with tarfile.open(LOCAL/'results.tgz','r:gz') as tar:
        for member in tar.getmembers():
            target=(root/member.name).resolve()
            if not target.is_relative_to(root) or member.issym() or member.islnk():raise ValueError('Unsafe archive')
        tar.extractall(root,filter='data')
    for path in (LOCAL/'bundle').rglob('*'):
        if path.is_file():assert a.sha(path)==a.sha(root/'bundle'/path.relative_to(LOCAL/'bundle'))
    jobs=[];points=[]
    for c in a.read(OUT/'native_contract.json')['cases']:
        raw=root/'bundle'/c['case']
        assert int((raw/'exit_code.txt').read_text())==0
        log=(raw/'console.log').read_text(errors='replace')
        assert 'T-2022.03-SP2' in log and 'Good Bye' in log
        curve=pltrows(raw/'IdVg_native_des.plt')
        assert len(curve)==21,(c['case'],len(curve))
        for index,vg in enumerate(c['vg']):
            match=[r for r in curve if abs(r['gate OuterVoltage']-vg)<=1e-10 and abs(r['drain OuterVoltage']-c['vd'])<=1e-10]
            assert len(match)==1,(vg,len(match))
            r=match[0]
            currents={name:r[name+' TotalCurrent'] for name in ('source','drain','gate','substrate')}
            points.append(dict(case=c['case'],device=c['device'],vd=c['vd'],vg=vg,current_A_per_um=currents['drain'],
                kcl_A_per_um=sum(currents.values()),kcl_over_Id=abs(sum(currents.values()))/max(abs(currents['drain']),1e-300),
                native_converged=True))
            source=raw/f'vg_{index:03d}_des.tdr';assert source.exists(),source
            jobs.append(dict(case=c['case'],index=index,tdr=str(source),export=str(LOCAL/'native_exports'/c['case']/f'vg_{index:03d}')))
    a.write_csv(OUT/'native_points.csv',points)
    a.write(OUT/'export_contract.json',dict(jobs=jobs,input_hashes={a.rel(LOCAL/'results.tgz'):a.sha(LOCAL/'results.tgz'),a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())},
        interpretation='Native DC convergence and target voltage verified; terminal KCL reported independently, no silent assertion of Vela acceptance.'))
    print('Downloaded native run verified:',len(jobs),'target states',flush=True)


def export_one(job):
    dest=Path(job['export']);dest.mkdir(parents=True,exist_ok=True)
    if (dest/'field_manifest.json').exists():return
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    done=subprocess.run([str(p.d.REPO/'build-release/sentaurus_import.exe'),'--tdr',job['tdr'],'--export-dir',str(dest),'--inventory-json',str(dest/'inventory.json')],env=env,capture_output=True,text=True)
    (dest/'export.log').write_text(done.stdout+done.stderr)
    assert done.returncode==0,done.stderr[-1000:]
    print(job['case'],job['index'],'exported',flush=True)


def export():
    a.verify(OUT/'export_contract.json')
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(export_one,a.read(OUT/'export_contract.json')['jobs']))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('unpack','export'))
    globals()[parser.parse_args().action]()
