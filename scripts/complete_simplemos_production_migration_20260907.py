"""Independent acceptance adapter for the production probe's smaller JSON schema.

The frozen original attempt remains intact. Acceptance is reconstructed from
the exported solved rows and conservative edge fluxes with identical gates.
"""
from validate_simplemos_production_migration_20260907 import *
import validate_simplemos_production_migration_20260907 as original
import json

def edge_config(dest):
    path=dest/'acceptance_edges.json'
    if not path.exists():
        q=a.read(dest/'all_row.json');q.pop('carrier_term_probe');q.update(simulation_type='sg_edge_flux_probe',output_csv=str(dest/'acceptance_edges.csv'))
        a.write(path,q)
    return path

def qualify(c,dest):
    s=a.read(dest/'config.status.json');q=a.read(dest/'all_row.status.json');geo,mask=m.previous.prior.support(c)
    path=edge_config(dest);runner=RUNNER;e=environment()
    if 'simplemos_native_mobility_candidate_20260907' in dest.parts:
        root=p.REPO/'build-release/simplemos_native_mobility_candidate_20260907';runner=root/'runner.exe'
        alpha={'zero':0.,'plus_full':.001,'minus_full':-.001,'plus_half':.0005,'minus_half':-.0005,'replacement':1.,'replacement_native':1.}[dest.name]
        e.update(VELA_CANDIDATE_MOBILITY_ALPHA=format(alpha,'.17g'),VELA_CANDIDATE_MOBILITY_RATIOS=str(root/'ratios'/(c['device']+'.txt')))
    es=execute(path,runner,e);assert es['exit_code']==0
    terms=d.ordered(dest/'all_row.csv',geo.count);edges=a.rows(dest/'acceptance_edges.csv');contacts=set(int(i) for i in geo.contact_nodes)
    rr=[];closure={}
    for r,keep in zip(terms,mask):
        if keep:
            for car in ('electron','hole'):
                scale=max(float(r[car+'_flux_abs_sum']),abs(float(r[car+'_recombination'])),abs(float(r[car+'_impact'])))
                rr.append(abs(float(r[car+'_residual']))/scale if scale else math.inf)
    assert len(rr)==1814
    for car in ('electron','hole'):
        flux=math.fsum((int(int(x['node0']) in contacts)-int(int(x['node1']) in contacts))*float(x[car+'_flux']) for x in edges)
        source=math.fsum(float(r[car+'_recombination'])+float(r[car+'_impact']) for i,r in enumerate(terms) if i not in contacts)
        ratio=abs(flux-source)/max(abs(flux),abs(source),1e-10)
        closure[car]=dict(contact_flux=flux,integrated_source=source,ratio=ratio,qualified=abs(source)>=1e-10,satisfied=abs(source)<1e-10 or ratio<=1e-6)
    bad=sum(x>1e-6 for x in rr);cc=s.get('contact_currents_A_per_um',{});Id=cc.get('drain',math.nan);kcl=abs(math.fsum(cc.values()))/abs(Id)
    result=dict(qualified=s['exit_code']==q['exit_code']==0 and s.get('converged',False) and bad==0 and all(x['satisfied'] for x in closure.values()) and kcl<=1e-8,
        current_A_per_um=Id,max_row_ratio=max(rr),row_violations=bad,kcl_over_Id=kcl,iterations=s.get('iterations'),failure=s.get('failure_reason',''),
        global_electron_ratio=closure['electron']['ratio'],global_hole_ratio=closure['hole']['ratio'])
    record=dict(method='Independent all 1814 rows; unbounded carrier contact edge sum vs original SRH source; unchanged 1e-6 / source_floor 1e-10 qualification.',result=result,closure=closure)
    target=dest/'independent_acceptance.json'
    if target.exists():assert a.read(target)==record
    else:a.write(target,record)
    return result

def post_config(cfg,dest):
    original.post_config(cfg,dest);edge_config(dest)

def addendum():
    a.verify(OUT/'freeze.json');files=[Path(__file__).resolve()]
    for c in a.read(OUT/'contract.json')['cases']:
        for j in c['jobs']:files.append(edge_config(Path(j['dest'])))
    a.write(OUT/'acceptance_addendum.json',dict(status='frozen_before_resume',retained_failure='Original driver expected isolated acceptance JSON fields absent from production carrier-term probe. All three completed DC replays and outputs retained.',
        method='Independent reconstruction of same acceptance from production solved rows plus conservative edge exports. No solver or acceptance threshold changes.',additional_read_only_probes=32))
    d.matrix.freeze(OUT/'acceptance_addendum_freeze.json',files+[OUT/'freeze.json',OUT/'acceptance_addendum.json'])

original.qualify=qualify
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('addendum','run','analyze'));action=parser.parse_args().action
    if action=='addendum':addendum()
    else:
        a.verify(OUT/'acceptance_addendum_freeze.json');getattr(original,action)()
