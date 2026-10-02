"""Independent exact-rational check of every C++ projected hi/lo coordinate."""
from fractions import Fraction
import json
from pathlib import Path
import h5py
from prepare_simplemos_three_state_20260928 import OUT, EXPORT, scalar

def read(p):
    with h5py.File(p) as f:
        return {k:v[()] for k,v in f['fields'].items()},json.loads(f.attrs['metadata_json'])
def q(v):return Fraction.from_float(float(v))
def main():
    count=0
    for device in ('n23','n24'):
        original,meta=read(OUT/f'{device}_a.h5')
        scale=q(meta['packed_potential_scale_V'])
        native={k:scalar(EXPORT/f'{device}_vd_0p05_vg_0p05',n) for k,n in
            [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential')]}
        for arm in ('b','c'):
            fields,derived=read(OUT/f'{device}_{arm}.h5')
            for key in ('mesh_sha256','split_mesh_fingerprint','packed_potential_scale_V','potential_origin_V'):
                assert meta[key]==derived[key]
            for block,(physical,packed) in enumerate(zip(('psi','phin','phip'),
                ('packed_psi','packed_electron_qf_increment','packed_hole_qf_increment'))):
                car='electron' if block==1 else 'hole'
                for i in range(len(fields[physical])):
                    value=native[physical].get(i,float(original[physical][i])) if arm=='c' else float(original[physical][i])
                    ref=q(original[car+'_qf_reference_V'][i]) if block else Fraction(0)
                    coordinate=(q(value)-ref)/scale
                    hi=float(coordinate);lo=float(coordinate-q(hi))
                    assert hi==fields[packed][i] and lo==fields[packed+'_low'][i],(device,arm,physical,i)
                    represented=(q(hi)+q(lo))*scale
                    assert float(represented+ref)==fields[physical][i]
                    if block:assert float(represented)==fields[car+'_qf_increment_V'][i]
                    count+=1
    result=dict(coordinates=count,all_exact_rational_rounding_checks_passed=True)
    (OUT.parent/'repack_independent_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
