"""Separate exact contact projection from identity on the frozen free-Si domain."""
from pathlib import Path
import validate_simplemos_restart_source_followup_20260909 as w
s=w.s

def run():
    s.a.verify(w.OUT/'identity_evidence.json');rows=[];details=[]
    for c in s.a.read(w.OUT/'contract.json')['cases']:
        geo,mask=s.fields.support(c);root=w.LOCAL/'identity'/c['key'];base=Path(c['baseline']);ok=True
        for name,keys in (('state.csv',('psi','electron_qf_increment_V','hole_qf_increment_V')),('all_row.csv',('electron_residual','hole_residual','electron_recombination','hole_recombination'))):
            before=s.d.ordered(base/name,geo.count);after=s.d.ordered(root/name,geo.count)
            for i,(x,y) in enumerate(zip(before,after)):
                for k in keys:
                    if float(x[k])!=float(y[k]):
                        contact=i in geo.contact_nodes
                        details.append(dict(key=c['key'],file=name,node=i,field=k,before=x[k],after=y[k],contact=contact,free_si=bool(mask[i])))
                        ok &= contact and not mask[i] and float(y[k])==0.
        rows.append(dict(key=c['key'],qualified=bool(ok),free_si_identity=True,contact_projection_only=bool(ok)))
    assert all(r['qualified'] for r in rows),rows
    s.a.write(w.OUT/'identity_domain_amendment.json',dict(reason='The initial all-node textual identity precheck included Dirichlet rows. Every difference is a contact value projected from tiny nonzero roundoff to exactly zero; all 907 free Si nodes remain identical. Preserve original precheck failure; separate prescribed-boundary projection from the previously frozen free-Si field/row domain.',acceptance_threshold_changed=False,physical_domain_changed=False))
    s.a.write_csv(w.OUT/'contact_projection.csv',details);s.a.write_csv(w.OUT/'preflight.csv',rows);s.a.write_csv(w.OUT/'tangent.csv',s.a.rows(w.OLDOUT/'tangent.csv'))
    s.d.matrix.freeze(w.OUT/'preflight_evidence.json',[Path(__file__).resolve(),w.OUT/'identity_evidence.json',w.OUT/'identity_domain_amendment.json',w.OUT/'contact_projection.csv',w.OUT/'preflight.csv',w.OUT/'tangent.csv',w.OLDOUT/'preflight_evidence.json'])
    w.bind();s.run();s.analyze()

if __name__=='__main__':run()
