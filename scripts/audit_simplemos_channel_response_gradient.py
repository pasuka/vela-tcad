"""Contract the frozen current gradient with the actual signed DC state response."""
from decimal import Decimal
import math
import validate_simplemos_vela_channel_charge as run


def audit():
    run.audit.verify(run.FREEZE)
    path = run.native.upstream.LOCAL / 'm65_n23_vd_0p050000_endpoint/adjoint.csv'
    gradients = {int(r['node_id']): r for r in run.audit.rows(path)}
    ledger = run.audit.rows(run.OUT / 'calibration.csv')
    results = []
    nodes = []
    for suffix, amp in (('1e13',1e13),('5e12',5e12)):
        states = [{int(r['node_id']):r for r in run.audit.rows(run.LOCAL/name/'state.csv')}
                  for name in ('plus_'+suffix,'minus_'+suffix)]
        terms = {'psi':[], 'phin':[], 'phip':[]}
        for i,g in gradients.items():
            p,m = (s[i] for s in states)
            for k, prefix in (('psi',None),('phin','electron'),('phip','hole')):
                if prefix:
                    # Same convention as packReferencedSolution; preserve small increments.
                    ref = prefix+'_qf_reference_V'; inc = prefix+'_qf_increment_V'
                    delta = float(Decimal(p[ref])-Decimal(m[ref])+Decimal(p[inc])-Decimal(m[inc]))/2
                else:
                    delta = float(Decimal(p[k])-Decimal(m[k]))/2
                term = float(g[f'dI_d{k}_A_per_um_per_V'])*delta
                terms[k].append(term)
                nodes.append({'amplitude_cm_3':amp,'node_id':i,'block':k,'state_central_response_V':delta,'gradient_contribution_A_per_um':term})
        current = next(r for r in ledger if float(r['amplitude_cm_3'])==amp)
        actual = float(current['central_delta_A_per_um'])
        adjoint = float(current['ifm_delta_A_per_um'])
        predicted = math.fsum(v for part in terms.values() for v in part)
        results.append({'amplitude_cm_3':amp,'actual_fd_A_per_um':actual,'gradient_dot_actual_state_A_per_um':predicted,
            'gradient_fd_relative_error':abs(predicted/actual-1),
            'adjoint_source_prediction_A_per_um':adjoint,
            'gradient_minus_adjoint_A_per_um':predicted-adjoint,
            **{k+'_contribution_A_per_um':math.fsum(v) for k,v in terms.items()}})
    run.audit.write_csv(run.OUT/'actual_direction_gradient.csv',results)
    run.audit.write_csv(run.OUT/'actual_direction_gradient_nodes.csv',nodes)
    for r in results: print(r)


if __name__=='__main__':
    audit()
