"""Export completed EP snapshots; compare same-node fields without interpolation."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_ep_20260928'
OLD=ROOT/'build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst/state_exports'
GATES=(0.,.05,.1,.15,.8)
FIELDS=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','eDensity','hDensity','srhRecombination','eMobility','hMobility')

def csvrows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def scalar(p):return {int(r['node_id']):float(r['component0']) for r in csvrows(p)}
def tag(v):return format(v,'.12g').replace('.','p')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    status=dict(line.split() for line in (OUT/'raw/status.tsv').read_text().splitlines())
    manifest=json.loads((OUT/'manifest.json').read_text())
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin;'+env.get('PATH','')
    importer=ROOT/'build/sentaurus_import.exe'
    exports=[]
    for c in manifest['cases']:
        if status.get(c['name'])!='0':continue
        tdrs=sorted((OUT/'raw/bundle'/c['device']).glob(c['name']+'_state_*_des.tdr'))
        assert len(tdrs)==5,(c,tdrs)
        for vg,tdr in zip(GATES,tdrs):
            dest=OUT/'state_exports'/c['algorithm']/f"{c['device']}_vg_{tag(vg)}"
            if not (dest/'field_manifest.json').exists():
                dest.mkdir(parents=True,exist_ok=True)
                with (dest/'import.stdout.json').open('w') as out,(dest/'import.stderr.txt').open('w') as err:
                    subprocess.run([str(importer),'--tdr',str(tdr),'--export-dir',str(dest)],stdout=out,stderr=err,check=True,env=env)
            exports.append(dict(device=c['device'],algorithm=c['algorithm'],vg=vg,tdr_sha256=sha(tdr),export=dest.relative_to(OUT).as_posix()))
    rows=[]
    for e in exports:
        if e['algorithm']!='default':continue
        a=OUT/e['export']
        candidates=[('EP_direct',OUT/'state_exports/direct'/a.name),
            ('M60_default',OLD/'default'/f"{e['device']}_vd_0p05_vg_{tag(e['vg'])}")]
        for label,b in candidates:
            if not b.exists():continue
            assert csvrows(a/'nodes.csv')==csvrows(b/'nodes.csv'),(a,b)
            for field in FIELDS:
                x,y=scalar(a/'fields'/f'{field}_region0.csv'),scalar(b/'fields'/f'{field}_region0.csv')
                assert x.keys()==y.keys()
                diff=[abs(x[i]-y[i]) for i in x]
                denom=math.fsum(abs(v) for v in y.values())
                rows.append(dict(device=e['device'],vg=e['vg'],comparison='EP_default_minus_'+label,
                    field=field,node_group='all_Si_including_contacts',nodes=len(x),
                    max_abs=max(diff),unweighted_L1_relative=math.fsum(diff)/denom if denom else None,
                    max_relative_nonzero=max((abs(x[i]/y[i]-1) for i in x if y[i]!=0),default=None),
                    zero_or_sign_change_nodes=sum(x[i]*y[i]<=0 and x[i]!=y[i] for i in x)))
    with (OUT/'field_comparison.csv').open('w',newline='') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'field_exports.json').write_text(json.dumps(dict(importer_sha256=sha(importer),exports=exports),indent=2)+'\n')
    print(json.dumps(dict(exported_states=len(exports),comparisons=len(rows))))

if __name__=='__main__':main()
