"""Read existing frozen cold/warm states; write only this experiment's outputs."""
import hashlib
import json
import math
from decimal import localcontext
from pathlib import Path
import sys

BASE=Path('/workspaces/simplemos-hfs-merged-20260926')
SOURCE=BASE/'source-continuation-v2'
OUT=BASE/'ep_review_20260928'
sys.path.insert(0,str(SOURCE/'scripts'))
import check_simplemos_cold_overlap_20260927 as check

def main():
    binary=SOURCE/'build/vela_example_runner'
    assert hashlib.sha256(binary.read_bytes()).hexdigest()=='87b934c79451110467cee69b62f6ec2e312a8aea51ab85c07bf2c0577ad05ea3'
    controls=BASE/'continuation_validation_20260927-v2/controls'
    matrix=BASE/'continuation_validation_20260927-v2/original_matrix'
    selected=check.h.read(controls/'finite/summary/selected.json')
    rows=[]
    for warm in selected:
        # Frozen matrix uses 0p05, while the old comparator silently skipped
        # low-Vd directories because it constructed the name with 0.05.
        case=f"{warm['device']}_vd_{warm['vd']:g}".replace('.','p')
        gate=matrix/case/'gate'
        if not gate.is_dir(): raise FileNotFoundError(gate)
        index=round(warm['vg']/.05)
        cold=check.h.read(gate/f'audit_{index:03d}/result.json')
        path=gate/('state_bias_'+format(warm['vg'],'.6f').replace('.','p')+'.h5')
        prior=Path(warm['dest'])/'state.h5'
        a,b=check.load(path),check.load(prior)
        assert a[1]['mesh_sha256']==b[1]['mesh_sha256']
        free=check.h.read(BASE/'original_inputs'/warm['device']/'geometry.json')['free_si']
        with localcontext() as ctx:
            ctx.prec=100
            diff={f+'_max_V':max(float(abs(check.potential(a,f,i)-check.potential(b,f,i))) for i in free) for f in ('psi','phin','phip')}
        density=max(abs(float(a[0][f][i])/float(b[0][f][i])-1) for i in free for f in ('electrons_m3','holes_m3'))
        current=abs(cold['current_A_per_um']/warm['current_A_per_um']-1)
        ok=(cold['qualified'] and warm['qualified'] and all(math.isfinite(v) and v<=1e-6 for v in diff.values()) and
            math.isfinite(density) and density<=1e-4 and math.isfinite(current) and current<=1e-6)
        rows.append(dict(case=case,vg=warm['vg'],arm=warm['arm'],**diff,density_max_relative=density,
            Id_relative=current,qualified=bool(ok),cold_state_sha256=check.h.sha(path),warm_state_sha256=check.h.sha(prior)))
    assert len(rows)==len(selected)==32
    result=dict(comparisons=len(rows),qualified=sum(r['qualified'] for r in rows),
        expected=32,passed=all(r['qualified'] for r in rows),
        coverage='n19/n23 x Vd 0.05/1 x Vg 0/0.2/0.8/1 x two initialization arms')
    check.h.csvout(OUT/'complete_cold_overlap.csv',rows)
    check.h.write(OUT/'complete_cold_overlap.json',result)
    print(json.dumps(result))
    assert result['passed'],'Original cold/warm gates failed; thresholds unchanged'

if __name__=='__main__':main()
