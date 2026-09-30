"""Repeat the temperature controls with a documented post-Load parameter ramp."""
import argparse,math,shutil,tarfile
from pathlib import Path
import simplemos_phumob_floor_temperature_20260909 as old

a,d=old.a,old.d
LOCAL=old.LOCAL/'ramped';OUT=old.OUT/'ramped';REMOTE=old.REMOTE+'/ramped'


def prepare():
    a.verify(old.OUT/'native_evidence.json')
    contract=a.read(old.OUT/'native_contract.json')
    files=[Path(__file__).resolve(),old.OUT/'temperature_preflight_failure.csv',old.OUT/'native_evidence.json']
    for job in contract['jobs']:
        src=old.LOCAL/'bundle/phumob'/job['key'];dest=LOCAL/'bundle/phumob'/job['key']
        dest.mkdir(parents=True,exist_ok=False)
        for path in src.iterdir():
            if path.name!='native_des.cmd':shutil.copy2(path,dest/path.name);files += [path,dest/path.name]
        text=(src/'native_des.cmd').read_text().split('Solve {')[0]
        text=text.replace(f"Temperature={job['temperature_K']}",'Temperature=300')
        text+=f'''Solve {{ Load(FilePrefix="result_050") Coupled {{ Poisson Electron Hole }}
 Quasistationary(InitialStep=0.25 MinStep=1e-6 MaxStep=0.5
 Goal {{ Model=DeviceTemperature Parameter="Temperature" Value={job['temperature_K']} }})
 {{ Coupled {{ Poisson Electron Hole }} }}
 Plot(FilePrefix="final") Save(FilePrefix="final") }}
'''
        (dest/'native_des.cmd').write_text(text,newline='\n');files.append(dest/'native_des.cmd')
    shutil.copy2(old.LOCAL/'run.sh',LOCAL/'run.sh');files.append(LOCAL/'run.sh')
    contract.update(remote=REMOTE,amendment='Start at the saved 300 K state; ramp documented DeviceTemperature/Temperature after Load. Preserve original failed temperature preflight and all ramp points. Final exported temperature remains a mandatory gate.')
    a.write(OUT/'native_contract.json',contract);files.append(OUT/'native_contract.json')
    d.matrix.freeze(OUT/'native_freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')


def unpack():
    a.verify(OUT/'native_freeze.json');dest=LOCAL/'native_raw';assert not dest.exists()
    with tarfile.open(LOCAL/'results.tgz') as tar:
        for m in tar.getmembers():assert (dest/m.name).resolve().is_relative_to(dest.resolve()) and not m.issym() and not m.islnk()
        tar.extractall(dest,filter='data')
    for path in (LOCAL/'bundle').rglob('*'):
        if path.is_file():assert a.sha(path)==a.sha(dest/'bundle'/path.relative_to(LOCAL/'bundle'))
    rows=[]
    for job in a.read(OUT/'native_contract.json')['jobs']:
        path=dest/'bundle/phumob'/job['key'];log=(path/'console.log').read_text(errors='replace')
        code=int((path/'exit_code.txt').read_text());row=dict(job,exit_code=code,native_qualified=False)
        if (path/'native_des.plt').exists():
            points=old.p.b.exporter.pltrows(path/'native_des.plt');point=points[-1]
            currents=[point[n+' TotalCurrent'] for n in ('drain','source','gate','substrate')]
            kcl=abs(math.fsum(currents))/max(abs(currents[0]),1e-300)
            bias=max(abs(point['gate OuterVoltage']-job['vg']),abs(point['drain OuterVoltage']-job['vd']))
            snapshots=sorted(path.glob('runtime_*_vertices.csv'))
            row.update(points=len(points),Id_A_per_um=currents[0],kcl_over_Id=kcl,bias_error_V=bias,
                runtime_prefix=snapshots[-1].name.removesuffix('_vertices.csv') if snapshots else '',
                native_qualified=code==0 and 'T-2022.03-SP2' in log and 'Good Bye' in log and kcl<=1e-8 and bias<=1e-10 and bool(snapshots))
        rows.append(row)
    a.write_csv(OUT/'native_points.csv',rows)
    d.matrix.freeze(OUT/'native_evidence.json',[OUT/'native_freeze.json',OUT/'native_points.csv',LOCAL/'results.tgz']+[p for p in dest.rglob('*') if p.is_file()])
    print(rows,flush=True)


def export():
    old.previous.LOCAL,old.previous.OUT=LOCAL,OUT;old.previous.export()


def analyze():
    import analyze_simplemos_phumob_floor_temperature_20260909 as analysis
    analysis.t.LOCAL,analysis.t.OUT=LOCAL,OUT;analysis.OUT=OUT/'analysis';analysis.main()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','unpack','export','analyze'))
    globals()[parser.parse_args().action]()
