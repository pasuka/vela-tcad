"""Supplementary zero-drift gate and input guards, frozen before DC results."""
import argparse
import math
from pathlib import Path
import calibrate_simplemos_srh_hotspots_20260908 as s

a=s.a;d=s.d;OUT=s.OUT/'qualification';LOCAL=s.LOCAL/'guards'


def prepare():
    a.verify(s.OUT/'freeze.json');case=a.read(s.OUT/'contract.json')['cases'][0];base=(s.LOCAL/(case['device']+'_sites.txt')).read_text().splitlines();N=int(base[0].split()[0]);files=[Path(__file__).resolve(),s.OUT/'freeze.json'];jobs=[]
    bad={
        'wrong_mesh':[f'{N+1} 4']+base[1:],
        'wrong_count':[f'{N} 3']+base[1:],
        'missing_record':base[:-1],
        'duplicate':base[:2]+[base[1]]+base[3:],
        'zero_ratio':[base[0],base[1].split()[0]+' 0']+base[2:],
        'out_of_range_node':[base[0],f'{N} 1']+base[2:],
        'extra':base+['unexpected'],
        'bad_alpha':base,
        'missing_alpha':base,
    }
    for name,lines in bad.items():
        dest=LOCAL/name;dest.mkdir(parents=True,exist_ok=False);(dest/'sites.txt').write_text('\n'.join(lines)+'\n')
        cfg=a.read(s.LOCAL/case['key']/'fixed/zero/functional.json');cfg['residual_output_csv']=str(dest/'residual.csv');cfg['contact_edge_output_csv']=str(dest/'contact_edges.csv')
        a.write(dest/'config.json',cfg);files.extend([dest/'sites.txt',dest/'config.json']);jobs.append(dict(name=name,dest=str(dest)))
    a.write(OUT/'contract.json',dict(case=case,jobs=jobs,zero_Id_drift_dex=1e-5,combined='Final calibration requires original response/DC gates AND four zero-current-drift gates AND all nine malformed-input guards. No changed numerical threshold.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json'])


def guards():
    a.verify(OUT/'freeze.json');contract=a.read(OUT/'contract.json');case=contract['case'];rows=[]
    for job in contract['jobs']:
        dest=Path(job['dest']);e=s.env(case,0.);e['VELA_NEIGHBORHOOD_SRH_VOLUMES']=str(dest/'sites.txt')
        if job['name']=='bad_alpha':e['VELA_NEIGHBORHOOD_SRH_ALPHA']='2'
        if job['name']=='missing_alpha':e.pop('VELA_NEIGHBORHOOD_SRH_ALPHA')
        st=s.v.execute(dest/'config.json',s.RUNNER,e);r=dict(name=job['name'],exit_code=st['exit_code'],qualified=st['exit_code']!=0);rows.append(r);assert r['qualified'],job
    a.write_csv(OUT/'guards.csv',rows);d.matrix.freeze(OUT/'guards_evidence.json',[OUT/'freeze.json',OUT/'guards.csv']+[f for f in LOCAL.rglob('*') if f.is_file()]);print('Nine malformed-input guards passed',flush=True)


def analyze():
    a.verify(s.OUT/'evidence.json');a.verify(OUT/'guards_evidence.json');rows=[]
    dc=a.rows(s.OUT/'dc.csv')
    for case in a.read(s.OUT/'contract.json')['cases']:
        zero=next(float(r['current_A_per_um']) for r in dc if r['key']==case['key'] and r['label']=='zero');base=a.read(Path(case['baseline'])/'config.status.json')['contact_currents_A_per_um']['drain'];drift=abs(math.log10(zero/base))
        rows.append(dict(key=case['key'],zero_Id_drift_dex=drift,qualified=drift<=1e-5))
    a.write_csv(OUT/'zero_drift.csv',rows);summary=dict(qualified=a.read(s.OUT/'summary.json')['qualified'] and all(x['qualified'] for x in rows),zero_drift_passed=sum(x['qualified'] for x in rows),input_guards_passed=9)
    a.write(OUT/'summary.json',summary);d.matrix.freeze(OUT/'evidence.json',[s.OUT/'evidence.json',OUT/'guards_evidence.json',OUT/'zero_drift.csv',OUT/'summary.json']);print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','guards','analyze'));globals()[parser.parse_args().action]()
