"""Locate conservative, adjoint-weighted SG factor contributions on edges."""
import math
from pathlib import Path
import decompose_simplemos_calibrated_transport as d


def main():
    d.a.verify(d.FREEZE)
    names=('mobility','sg_secant_conductance','qf_log_imbalance')
    output,checks=[],[]
    for c in d.a.read(d.CONTRACT)['cases']:
        root=d.LOCAL/c['case']
        adj=d.ordered(Path(c['adjoint']),c['nodes'])
        old={int(r['edge_id']):r for r in d.a.rows(root/'strict/edges.csv')}
        new={int(r['edge_id']):r for r in d.a.rows(root/'mapped/edges.csv')}
        games={}
        for r in d.a.rows(root/'secant_edges.csv'):
            i=int(r['edge_id'])
            games.setdefault(i,{})[int(r['mask'])]=r
        rows=[]
        for edge,variants in games.items():
            assert set(variants)==set(range(8))
            r=variants[0]
            i,j=int(r['node0']),int(r['node1'])
            wi=float(adj[i]['lambda_electron']) if r['node0_constrained']=='0' else 0.
            wj=float(adj[j]['lambda_electron']) if r['node1_constrained']=='0' else 0.
            allocated=d.allocate({mask:float(v['electron_flux']) for mask,v in variants.items()},names)
            for factor,value in allocated.items():
                rows.append({'case':c['case'],'factor':factor,'edge_id':edge,'node0':i,'node1':j,
                    'midpoint_x_um':float(r['midpoint_x']),'midpoint_y_um':float(r['midpoint_y']),
                    'current_A_per_um':-(wi-wj)*value,
                    'strict_mobility_m2_V_s':float(old[edge]['electron_mobility_m2_V_s']),
                    'mapped_mobility_m2_V_s':float(new[edge]['electron_mobility_m2_V_s']),
                    'strict_drive_V_m':float(old[edge]['electron_mobility_drive_V_m']),
                    'mapped_drive_V_m':float(new[edge]['electron_mobility_drive_V_m']),
                    'strict_sg_log_imbalance':float(old[edge]['electron_sg_log_left_over_right']),
                    'mapped_sg_log_imbalance':float(new[edge]['electron_sg_log_left_over_right'])})
        targets=d.a.rows(d.OUT/'factor.csv')
        for factor in names:
            group=[r for r in rows if r['factor']==factor]
            total=math.fsum(r['current_A_per_um'] for r in group)
            target=float(next(r['current_A_per_um'] for r in targets if r['case']==c['case'] and r['factorization']=='secant' and r['factor']==factor))
            err=abs(total-target)/max(abs(target),1e-300)
            assert err<=1e-8
            checks.append({'case':c['case'],'factor':factor,'edge_sum_A_per_um':total,'node_sum_A_per_um':target,'relative_error':err})
            for sign in ('positive','negative'):
                chosen=[r for r in group if (r['current_A_per_um']>0 if sign=='positive' else r['current_A_per_um']<0)]
                for rank,row in enumerate(sorted(chosen,key=lambda r:abs(r['current_A_per_um']),reverse=True)[:10],1):
                    output.append(dict(row,sign=sign,rank=rank))
        d.a.write_csv(root/'edge_factor_localization.csv',rows)
    d.a.write_csv(d.OUT/'top_edges.csv',output)
    d.a.write_csv(d.OUT/'edge_localization_checks.csv',checks)
    d.a.write(d.OUT/'edge_localization_evidence.json',{'input_hashes':{d.a.rel(Path(__file__).resolve()):d.a.sha(Path(__file__).resolve())},'new_runner_calls':0})
    print('Located ten positive and negative edges per factor; all edge/node factor sums close',flush=True)


if __name__=='__main__':main()
