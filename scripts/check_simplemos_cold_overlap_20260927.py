"""Compare cold-sweep checkpoints with both qualified initialization arms."""
import argparse,json,math
from decimal import Decimal,localcontext
from pathlib import Path
import h5py
import simplemos_hfs_cloud_20260926 as h

def load(path):
    with h5py.File(path) as f:
        fields={k:v[()] for k,v in f['fields'].items()}
        metadata=json.loads(f.attrs['metadata_json'])
    return fields,metadata

def potential(state,field,node):
    fields,meta=state
    prefix={'psi':'packed_psi','phin':'packed_electron_qf_increment','phip':'packed_hole_qf_increment'}[field]
    d=lambda x:Decimal.from_float(float(x))
    value=(d(fields[prefix][node])+d(fields[prefix+'_low'][node]))*d(meta['packed_potential_scale_V'])
    if field!='psi':
        value+=d(fields[('electron' if field=='phin' else 'hole')+'_qf_reference_V'][node])
    return value

def compare(controls,matrix,inputs):
    selected=h.read(controls/'finite/summary/selected.json');out=[]
    for warm in selected:
        case=f"{warm['device']}_vd_{warm['vd']:g}"
        gate=matrix/case/'gate'
        if not gate.exists():continue
        index=round(warm['vg']/.05)
        cold=h.read(gate/f'audit_{index:03d}/result.json')
        path=gate/('state_bias_'+format(warm['vg'],'.6f').replace('.','p')+'.h5')
        prior=Path(warm['dest'])/'state.h5'
        a,b=load(path),load(prior)
        assert a[1]['mesh_sha256']==b[1]['mesh_sha256']
        free=h.read(inputs/warm['device']/'geometry.json')['free_si']
        with localcontext() as ctx:
            ctx.prec=100
            diff={field+'_max_V':max(float(abs(potential(a,field,i)-potential(b,field,i))) for i in free)
                  for field in ('psi','phin','phip')}
        density=max(abs(float(a[0][field][i])/float(b[0][field][i])-1) for i in free for field in ('electrons_m3','holes_m3'))
        current=abs(cold['current_A_per_um']/warm['current_A_per_um']-1)
        ok=(cold['qualified'] and warm['qualified'] and
            all(math.isfinite(v) and v<=1e-6 for v in diff.values()) and
            math.isfinite(density) and density<=1e-4 and math.isfinite(current) and current<=1e-6)
        out.append(dict(case=case,vg=warm['vg'],arm=warm['arm'],**diff,
            density_max_relative=density,Id_relative=current,qualified=ok,
            cold_state_sha256=h.sha(path),warm_state_sha256=h.sha(prior)))
    expected=sum(1 for r in selected if (matrix/f"{r['device']}_vd_{r['vd']:g}"/'gate').exists())
    report=dict(comparisons=len(out),qualified=sum(r['qualified'] for r in out),
        passed=len(out)==expected and len(out)>0 and all(r['qualified'] for r in out))
    h.csvout(matrix/'cold_overlap.csv',out);h.write(matrix/'cold_overlap.json',report)
    assert report['passed'],'Cold/warm overlap failed frozen dual-initialization gates'
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('controls','matrix','inputs'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(compare(a.controls,a.matrix,a.inputs))
