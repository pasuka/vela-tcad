"""Keep all-row failure and initial-guess dependence separate from qualification."""
from pathlib import Path
from decimal import Decimal
import itertools
import math
import numpy as np
import validate_simplemos_minority_invariance_20260906 as m

a=m.a;d=m.d


def potentials(rows,carrier):
    return np.array([float(Decimal(r[carrier+'_qf_reference_V'])+Decimal(r[carrier+'_qf_increment_V'])) for r in rows])


def complete_qualified_group(group):
    return (len(group)==4
        and {r[0]['initialization'] for r in group}=={'vela','native','hole_plus','hole_minus'}
        and all(r[0]['strict_qualification_with_visibility'] for r in group))


def main():
    a.verify(m.OUT/'freeze.json');a.verify(m.OUT/'screen_freeze.json')
    contract=a.read(m.OUT/'contract.json');runs=[];violations=[];pairs=[];groups={};missing=[]
    files=[Path(__file__).resolve(),m.OUT/'freeze.json',m.OUT/'screen_freeze.json',m.OUT/'screen_runs.csv',m.OUT/'interrupted_processes.json']
    if (m.OUT/'deferred.csv').exists():files.append(m.OUT/'deferred.csv')
    for job in contract['jobs']:
        dest=Path(job['config']).parent
        groups.setdefault((job['case'],job['vg']),[])
        if not (dest/'result.json').exists():
            missing.append(dict(case=job['case'],vg=job['vg'],initialization=job['initialization'],
                status='no_completed_result',interpretation='Deferred by unshifted prerequisite; any interruption is recorded separately.'))
            continue
        result=a.read(dest/'result.json');status=a.read(dest/'config.status.json')
        rows=d.ordered(dest/'state.csv',status['nodes']);terms=d.ordered(dest/'all_row.csv',status['nodes'])
        geo=d.matrix.spatial.m73.Geometry(job['device']);mask=np.zeros(geo.count,dtype=bool)
        masks,_,_=d.matrix.spatial.old.m78.supports(job['device'],geo,.05);mask=masks['all_si'].copy();mask[geo.contact_nodes]=False
        allstatus=a.read(dest/'all_row.status.json');legacy=a.read(dest/'legacy.status.json')
        zero=sum(max(abs(float(r[c+'_flux_abs_sum'])),abs(float(r[c+'_recombination'])),abs(float(r[c+'_impact'])))==0 for i,r in enumerate(terms) if mask[i] for c in ('electron','hole'))
        row={k:result[k] for k in ('case','device','vd','vg','initialization','exit_code','converged','iterations','current_A_per_um','all_row_qualified','legacy_qualified','kcl_over_Id','all_row_violations','all_row_max_ratio','all_row_checked_rows')}
        row.update(failure_reason=status.get('failure_reason',''),zero_scale_active_rows=zero,
            expected_active_rows=2*int(sum(mask)),legacy_post_gate_satisfied=legacy['carrier_row_convergence']['satisfied'] and legacy['global_continuity_closure']['satisfied'] and row['kcl_over_Id']<=1e-8,
            strict_qualification_with_visibility=row['all_row_qualified'] and zero==0)
        assert row['all_row_checked_rows']==row['expected_active_rows'],(job,row)
        row_violations=allstatus['carrier_row_convergence']['violations']
        assert len(row_violations)==row['all_row_violations']
        row.update(electron_violations=sum(v['carrier']=='electron' for v in row_violations),
            hole_violations=sum(v['carrier']=='hole' for v in row_violations),
            channel_hole_violations=sum(v['carrier']=='hole' and bool(masks['channel'][int(v['node_id'])]) for v in row_violations))
        runs.append(row)
        for v in row_violations:
            detail=dict(v);node=int(detail['node_id']);carrier=detail['carrier']
            density=float(rows[node]['electrons_m3' if carrier=='electron' else 'holes_m3'])/1e6
            raw_density=detail.pop('carrier_density_m3')
            assert math.isclose(raw_density,density,rel_tol=1e-11,abs_tol=1e-290),(job,node,raw_density,density)
            violations.append(dict(case=job['case'],vg=job['vg'],initialization=job['initialization'],
                carrier_density_cm3=density,probe_density_raw=raw_density,
                probe_density_note='Original carrier_density_m3 JSON field carries internal cm^-3; checked against exported SI state.',
                channel=bool(masks['channel'][node]),gate_interface=bool(masks['gate_interface'][node]),**detail))
        groups.setdefault((job['case'],job['vg']),[]).append((row,rows,geo,masks))
        files += [dest/name for name in ('config.status.json','state.csv','all_row.csv','all_row.status.json','legacy.status.json','result.json')]
    case_results=[]
    for (key,vg),group in groups.items():
        qualified=complete_qualified_group(group);passes=[]
        for left,right in itertools.combinations(group,2):
            l,sl,geo,masks=left;r,sr,_,_=right
            mask=masks['all_si'].copy();mask[geo.contact_nodes]=False
            pn=float(max(abs(potentials(sl,'electron')[mask]-potentials(sr,'electron')[mask])))
            pp=float(max(abs(potentials(sl,'hole')[mask]-potentials(sr,'hole')[mask])))
            psi=float(max(abs(d.array(sl,('psi',))[0,mask]-d.array(sr,('psi',))[0,mask])))
            nl=d.array(sl,('electrons_m3','holes_m3'))[:,mask];nr=d.array(sr,('electrons_m3','holes_m3'))[:,mask]
            density=float(np.max(abs(nl-nr)/np.maximum(np.maximum(abs(nl),abs(nr)),1e-300)))
            current=abs(l['current_A_per_um']-r['current_A_per_um'])/max(abs(l['current_A_per_um']),abs(r['current_A_per_um']),1e-300)
            okay=max(psi,pn,pp)<=1e-6 and density<=1e-4 and current<=1e-6
            pair=dict(case=key,vg=vg,left=l['initialization'],right=r['initialization'],both_states_qualified=l['strict_qualification_with_visibility'] and r['strict_qualification_with_visibility'],
                psi_max_V=psi,phin_max_V=pn,phip_max_V=pp,density_max_symmetric_relative=density,Id_symmetric_relative=current,
                numerical_invariance=okay,invariance_qualified=okay and qualified)
            pairs.append(pair);passes.append(okay)
        case_results.append(dict(case=key,vg=vg,initializations=len(group),qualified_states=sum(r[0]['strict_qualification_with_visibility'] for r in group),
            all_states_qualified=qualified,numerical_invariance=bool(passes) and all(passes),invariance_qualified=qualified and bool(passes) and all(passes)))
    for name,data in (('qualification',runs),('violations',violations),('pairs',pairs),('cases',case_results),('not_completed',missing)):
        if data:a.write_csv(m.OUT/(name+'.csv'),data)
    a.write(m.OUT/'result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},
        planned_states=len(contract['jobs']),completed_states=len(runs),not_completed_states=len(missing),
        all_row_qualified=sum(r['strict_qualification_with_visibility'] for r in runs),
        legacy_post_gate_satisfied=sum(r['legacy_post_gate_satisfied'] for r in runs),invariance_qualified_cases=sum(r['invariance_qualified'] for r in case_results),
        boundary='All completed original runs are retained, including shifted runs completed before the prerequisite screen took effect. Differences between failed states do not establish nonuniqueness of converged physical solutions. No acceptance criterion or previous failure record is replaced.',production_changes=False))
    print(case_results,flush=True)


if __name__=='__main__':main()
