"""Independent decimal block-merit accounting for unchanged rejected trials."""
import decimal
from pathlib import Path
import audit_simplemos_packed_stalls_20260910 as v
a,d=v.a,v.d
D=decimal.Decimal.from_float

def changes(rows,precision):
    with decimal.localcontext() as c:
        c.prec=precision;n=len(rows)//3;block=[]
        for b in range(3):
            group=rows[b*n:(b+1)*n]
            delta=sum(((D(float(r['Rtrial']))-D(float(r['R'])))*(D(float(r['Rtrial']))+D(float(r['R']))) for r in group),decimal.Decimal(0))
            block.append(delta*D(float(group[0]['weight']))/D(float(group[0]['scale']))**2)
        return block,sum(block)

def main():
    a.verify(v.OUT/'replay_evidence.json');records=[];summary=[]
    for control in a.rows(v.OUT/'identity.csv'):
        dest=Path(control['dest']);own=[]
        for file in sorted(dest.glob('trial_*.csv'),key=lambda p:int(p.stem.split('_')[1])):
            rows=a.rows(file);n=len(rows)//3;b,total=changes(rows,90);bb,total2=changes(rows,140)
            assert (total>0)-(total<0)==(total2>0)-(total2<0)
            assert abs(total-total2)<=abs(total2)*decimal.Decimal('1e-70')
            rec={k:control[k] for k in ('device','vd','index')};rec['trial']=file.stem
            rec.update(psi_delta_square=str(bb[0]),electron_delta_square=str(bb[1]),hole_delta_square=str(bb[2]),weighted_delta=str(total2),strict_descent=total2<0)
            for b,name in enumerate(('psi','electron','hole')):
                group=rows[b*n:(b+1)*n]
                rec[name+'_changed']=sum(r['x']!=r['candidate'] for r in group)
                rec[name+'_raw_max']=max(abs(float(r['raw_step'])) for r in group)
                rec[name+'_cap_difference_max']=max(abs(float(r['raw_step'])-float(r['capped_step'])) for r in group)
                rec[name+'_linear_residual_max']=max(abs(float(r['JdxR'])) for r in group)
            records.append(rec);own.append(rec)
        updates=a.rows(dest/'updates.csv');hole=[r for r in updates if r['carrier']=='hole'][-1]
        summary.append(dict(device=control['device'],vd=control['vd'],index=control['index'],trials=len(own),strict_descents=sum(r['strict_descent'] for r in own),psi_increases=sum(decimal.Decimal(r['psi_delta_square'])>0 for r in own),hole_decreases=sum(decimal.Decimal(r['hole_delta_square'])<0 for r in own),node=hole['node_id'],raw_hole_update_V=hole['raw_linear_step_V'],capped_hole_update_V=hole['capped_step_V'],best_trial_hole_update_V=hole['applied_step_V'],current_hole_residual=hole['residual'],best_trial_hole_residual=hole['selected_trial_residual']))
    a.write_csv(v.OUT/'trial_merit.csv',records);a.write_csv(v.OUT/'summary.csv',summary)
    d.matrix.freeze(v.OUT/'analysis_evidence.json',[Path(__file__).resolve(),v.OUT/'replay_evidence.json',v.OUT/'trial_merit.csv',v.OUT/'summary.csv'])
    for r in summary:print(r,flush=True)

if __name__=='__main__':main()
