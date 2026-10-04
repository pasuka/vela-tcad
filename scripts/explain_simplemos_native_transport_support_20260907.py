"""Expose material shares on the largest calibrated-response edges."""
import ast
import math
from pathlib import Path
import validate_simplemos_distributed_transport_20260907 as v
p=v.p;a=p.a


def main():
    a.verify(p.OUT/'independent_check.json')
    contributions=a.rows(p.OUT/'distributed_edge_contributions.csv');rows=[]
    for c in a.read(v.OUT/'contract.json')['cases']:
        geo,xy,cells,vol,K,parts,info,modified=v.mobility.geom.geometry(c['device'])
        edges={int(x['edge_id']):x for x in a.rows(p.prior.LOCAL/'vela'/c['key']/'strict/edges.csv')}
        weights=v.mobility.geometry(c['device'])[1]
        for item in contributions:
            if item['key']!=c['key'] or int(item['rank_absolute_contribution'])>5:continue
            edge=edges[int(item['edge'])];pair=(int(edge['node0']),int(edge['node1']))
            g={mat:math.fsum(x['coefficient'] for x in parts[pair] if x['material']==mat) for mat in ('Si','SiO2','Nitride')}
            g0=float(edge['couple_m'])/float(edge['length_m']);mu0=float(edge['electron_mobility_m2_V_s'])*1e4
            assert g['Si']>0 and mu0>0
            rows.append(dict(key=c['key'],edge=item['edge'],node0=pair[0],node1=pair[1],rank=item['rank_absolute_contribution'],
                x_mid_um=item['x_mid_um'],y_mid_um=item['y_mid_um'],Vela_geometry=g0,native_Si_geometry=g['Si'],native_SiO2_geometry=g['SiO2'],native_Nitride_geometry=g['Nitride'],
                native_all_over_Vela_geometry=math.fsum(g.values())/g0,native_Si_over_Vela_geometry=g['Si']/g0,
                Vela_electron_mobility_cm2_Vs=mu0,native_coefficient_weighted_mobility_cm2_Vs=weights[pair]['e']/g['Si'],
                native_over_Vela_mu_geometry=weights[pair]['e']/(g0*mu0),unit_response_A_per_um=item['unit_total_A_per_um']))
    a.write_csv(p.OUT/'dominant_edge_material_support.csv',rows)
    script=Path(__file__).resolve();ast.parse(script.read_text(encoding='utf-8'))
    a.write(p.OUT/'material_support_evidence.json',dict(status='completed_read_only_support_explanation',
        input_hashes={a.rel(f):a.sha(f) for f in [script,p.OUT/'independent_check.json',p.OUT/'dominant_edge_material_support.csv']},
        interpretation='Per-edge adjoint rankings belong to the calibrated combined direction; individual edge finite differences were not performed.'))
    print('Recorded',len(rows),'dominant edges with explicit material geometry shares',flush=True)


if __name__=='__main__':main()
