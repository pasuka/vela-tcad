"""Read-only remaining-failure and zero-RHS backward-error interpretation."""
from pathlib import Path
import math
import validate_simplemos_linear_refinement_20260906 as m

a=m.a;d=m.d


def main():
    a.verify(m.OUT/'result.json');failures=[];boundaries=[]
    files=[Path(__file__).resolve(),m.OUT/'result.json',m.OUT/'comparison.csv',m.OUT/'last_trials.csv']
    comparisons={r['tag']:r for r in a.rows(m.OUT/'comparison.csv')}
    for job in a.read(m.OUT/'contract.json')['jobs']:
        dest=Path(job['config']).parent;data=a.rows(dest/'updates.csv');iteration=max(int(r['iteration']) for r in data)
        rows=[r for r in data if int(r['iteration'])==iteration];geo=d.matrix.spatial.m73.Geometry(job['device'])
        contacts=set(geo.contact_nodes);identity=[]
        for r in rows:
            if int(r['node_id']) in contacts and float(r['jacobian_diagonal'])==1 and float(r['residual'])==0 and float(r['raw_linear_residual'])!=0:
                # Contact rows are identity rows by assembler constrained-row construction.
                # Unit scaling uses thermal voltage for the potential scale; stored updates are physical volts.
                assert math.isclose(float(r['raw_linear_step_V']),float(r['raw_linear_residual'])*d.fixed.upstream.VT,rel_tol=1e-12,abs_tol=1e-290),r
                identity.append(r)
        if identity:
            r=max(identity,key=lambda r:abs(float(r['raw_linear_residual'])))
            boundaries.append(dict(tag=job['tag'],iteration=iteration,identity_zero_rhs_nonzero_defect_rows=len(identity),carrier=r['carrier'],node=r['node_id'],
                max_absolute_identity_defect=float(r['raw_linear_residual']),step_V=float(r['raw_linear_step_V']),
                backward_ratio=1.,interpretation='For F=0 and J(row,:)=unit row, abs(Jdx+F)/(abs(F)+abs(J)abs(dx)) is 1 for any nonzero roundoff step, independent of its absolute size.'))
        c=comparisons[job['tag']]
        if c['refined_qualified']=='False':
            gate=a.read(dest/'all_row.status.json')['carrier_row_convergence'];worst=max(gate['violations'],key=lambda r:float(r['ratio'])) if gate['violations'] and 'ratio' in gate['violations'][0] else None
            if worst is None:
                terms=a.rows(dest/'all_row.csv');candidates=[]
                mask=d.matrix.spatial.old.m78.supports(job['device'],geo,.05)[0]['all_si'].copy();mask[geo.contact_nodes]=False
                for t in terms:
                    node=int(t['node_id'])
                    if not mask[node]:continue
                    for carrier in ('electron','hole'):
                        scale=max(float(t[carrier+'_flux_abs_sum']),abs(float(t[carrier+'_recombination'])),abs(float(t[carrier+'_impact'])))
                        candidates.append(dict(node_id=node,carrier=carrier,ratio=abs(float(t[carrier+'_residual']))/scale))
                worst=max(candidates,key=lambda r:r['ratio'])
            target=next(r for r in rows if r['carrier']==worst['carrier'] and int(r['node_id'])==int(worst['node_id']))
            trials=[r for r in a.rows(m.OUT/'last_trials.csv') if r['tag']==job['tag']]
            last=a.rows(dest/'linear_summary.csv')[-1]
            failures.append(dict(tag=job['tag'],iteration=iteration,worst_carrier=worst['carrier'],worst_node=worst['node_id'],
                worst_raw_step_V=float(target['raw_linear_step_V']),worst_capped_step_V=float(target['capped_step_V']),
                capped_carrier_rows=sum(float(r['raw_linear_step_V'])!=float(r['capped_step_V']) for r in rows),
                worst_nonlinear_ratio=float(c['refined_max_ratio']),last_linear_max=float(last['max_carrier_linear_ratio']),
                last_linear_rows_above_1e_minus_8=int(last['rows_above_1e_minus_8']),trials=len(trials),
                all_finite_valid=all(r['finite']=='True' and r['caller_valid']=='True' for r in trials),
                all_rejected_non_decrease=all(r['accepted']=='False' and r['reason']=='line_search_non_decrease' for r in trials),
                poisson_merit_fraction=float(c['last_poisson_merit_fraction'])))
        files += [dest/'updates.csv',dest/'all_row.csv',dest/'all_row.status.json',dest/'linear_summary.csv']
    a.write_csv(m.OUT/'remaining_failures.csv',failures);a.write_csv(m.OUT/'zero_rhs_identity_rows.csv',boundaries)
    a.write(m.OUT/'detail_result.json',dict(input_hashes={a.rel(p):a.sha(p) for p in sorted(set(files))},failed_states=len(failures),
        identity_boundary_examples=len(boundaries),scope='Read-only interpretation; no new solve or change of monitor/acceptance.',
        boundary='Zero-RHS identity rows explain why the ALL-MATRIX relative backward maximum can remain 1 while absolute constraint defects are tiny; no rows are removed from that reported maximum.'))
    print(failures,flush=True);print(boundaries,flush=True)


if __name__=='__main__':main()
