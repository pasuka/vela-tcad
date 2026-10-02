"""Independent native temperature controls for the PhuMob G floor."""
import argparse,shutil,tarfile
from pathlib import Path
import simplemos_phumob_supported_export_20260908 as previous

p=previous.original;a,d=p.a,p.d
LOCAL=p.REPO/'build-release/phumob_floor_20260909'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909'
REMOTE='/tmp/vela_phumob_floor_20260909'


def prepare():
    a.verify(previous.OUT/'native_evidence.json')
    old=next(j for j in a.rows(previous.OUT/'native_points.csv') if j['model']=='phumob' and j['device']=='n19' and int(j['index'])==50 and float(j['vd'])==.05)
    src=previous.LOCAL/'bundle/phumob'/old['key'];files=[Path(__file__).resolve(),previous.OUT/'native_evidence.json'];jobs=[]
    for temperature in (299,300,301,350):
        job=dict(old);job.update(key=f'T{temperature}',temperature_K=temperature,vg=1.,vd=.05,index=50)
        for k in ('native_qualified','runtime_prefix','Id_A_per_um','exit_code','kcl_over_Id','bias_error_V'):job.pop(k,None)
        dest=LOCAL/'bundle/phumob'/job['key'];dest.mkdir(parents=True,exist_ok=False)
        for name in ('input_fps.tdr','result_050_des.sav','result_050_circuit_des.sav','runtime.tcl'):
            shutil.copy2(src/name,dest/name);files += [src/name,dest/name]
        text=(src/'native_des.cmd').read_text()
        assert text.count('Physics { EffectiveIntrinsicDensity')==1
        text=text.replace('Physics { EffectiveIntrinsicDensity',f'Physics {{ Temperature={temperature} EffectiveIntrinsicDensity')
        text=text.replace('Plot { eDensity','Plot { Temperature eDensity')
        (dest/'native_des.cmd').write_text(text,newline='\n');files.append(dest/'native_des.cmd');jobs.append(job)
    shell=(previous.LOCAL/'run.sh').read_text().replace('count % 4','count % 2')
    (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
    a.write(OUT/'native_contract.json',dict(jobs=jobs,remote=REMOTE,
        scope='Only prescribed lattice temperature changes on n19 Vd=.05 Vg=1; same PhuMob/OldSlotboom/SRH, native Math and geometry callback.',
        questions='Check whether inferred G-floor discrepancy persists under independent temperature changes; compare mathematical and open-source root algorithms without fitting production parameters.',
        temperature_export_required=True,gates=dict(temperature_K=1e-8,native_kcl_relative=1e-8,bias_V=1e-10,cell_mobility_relative=1e-7),
        production_parameter_fit_permitted=False))
    d.matrix.freeze(OUT/'native_freeze.json',files+[OUT/'native_contract.json'])
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')


def unpack():
    p.LOCAL,p.OUT=LOCAL,OUT;p.unpack()


def export():
    previous.LOCAL,previous.OUT=LOCAL,OUT;previous.export()


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','unpack','export'))
    globals()[parser.parse_args().action]()
