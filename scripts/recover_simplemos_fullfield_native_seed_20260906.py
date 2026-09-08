"""Strict independent reclosure from coherent native states for ascending gaps."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import math
import run_simplemos_fullfield_vela_20260906 as base

a=base.a;d=base.d;LOCAL=base.LOCAL;OUT=base.OUT


def prepare():
    a.verify(OUT/'native_freeze.json');a.verify(OUT/'export_contract.json');a.verify(OUT/'vela_freeze.json')
    previous=a.rows(OUT/'vela_ascending.csv');good={(r['case'],round(float(r['vg']),10)) for r in previous if r['qualified']=='True'}
    jobs=[];files=[Path(__file__).resolve(),OUT/'export_contract.json',OUT/'vela_ascending.csv',base.v.RUNNER]
    for c in a.read(OUT/'vela_contract.json')['cases']:
        geo=d.matrix.spatial.m73.Geometry(c['device'])
        for index in range(21):
            vg=round(.05*index,10)
            if (c['case'],vg) in good:continue
            export=LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
            psi,n,p,spread=d.fixed.upstream.coherent_state(export,geo)
            assert spread<1e-12 and d.matrix.spatial.m73.coordinate_error(export,geo)<=1e-10
            en=d.matrix.spatial.m73.scalar(export/'fields/eQuasiFermiPotential_region0.csv')
            hp=d.matrix.spatial.m73.scalar(export/'fields/hQuasiFermiPotential_region0.csv')
            dest=LOCAL/'native_seed'/c['case']/f'vg_{index:03d}'
            state=[dict(node_id=i,psi=psi[i],phin=en.get(i,0.),phip=hp.get(i,0.),electrons_m3=n[i],holes_m3=p[i]) for i in range(geo.count)]
            a.write_csv(dest/'initial.csv',state)
            cfg=a.read(Path(c['configs'][index]));cfg['state_file']=str(dest/'initial.csv');cfg['output_state_file']=str(dest/'state.csv')
            a.write(dest/'config.json',cfg)
            files += [dest/'initial.csv',dest/'config.json']+[f for f in (export/'fields').glob('*.csv')]
            jobs.append(dict(case=c['case'],vg=vg,config=str(dest/'config.json')))
    a.write(OUT/'recovery_contract.json',dict(status='frozen_before_execution',jobs=jobs,
        comparison_axis='Initial guess only: preserve native psi/phin/phip, reconstruct n/p with unchanged Vela ni/Vt; independently reclose exact Vela equations and original strict acceptance.',
        gates=a.read(OUT/'vela_contract.json')['gates'],ascending_failures_retained=True,
        qualification='Individual DC states; not a claim that the ascending sweep has become continuous or that mapped native input is itself a Vela solution.'))
    files.append(OUT/'recovery_contract.json');d.matrix.freeze(OUT/'recovery_freeze.json',files)
    print('Frozen',len(jobs),'native-seeded recovery states',flush=True)


def one(job):
    row=base.qualify(Path(job['config']).parent);row['initialization']='native_coherent'
    print(row,flush=True);return row


def run():
    a.verify(OUT/'recovery_freeze.json')
    with ThreadPoolExecutor(max_workers=2) as pool:rows=list(pool.map(one,a.read(OUT/'recovery_contract.json')['jobs']))
    a.write_csv(OUT/'vela_recovery.csv',rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'))
    globals()[parser.parse_args().action]()
