"""Exact Vela Poisson bookkeeping, conservative edge identities and local SG proxies."""
from pathlib import Path
from decimal import Decimal
import math
import numpy as np
import prepare_simplemos_masetti_local_20260907 as p
import audit_simplemos_local_electrochemical_fields as f

a=p.a;d=p.d;OUT=p.OUT;LOCAL=p.LOCAL


def qf(row,c='electron'):
    key=c+'_qf_reference_V'
    return Decimal(row[key])+Decimal(row[c+'_qf_increment_V']) if key in row else Decimal(row['phin' if c=='electron' else 'phip'])


def main():
    a.verify(OUT/'vela_freeze.json');cases=[];components=[];nodes=[];edges=[];checks=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge'];Q=d.fixed.Q
    anchor=d.matrix.spatial.m73.Geometry('n23')
    for c in a.read(OUT/'vela_contract.json')['cases']:
        root=LOCAL/'vela'/c['key'];geo=d.matrix.spatial.m73.Geometry(c['device']);N=geo.count
        masks,_,xy=d.matrix.spatial.old.m78.supports(c['device'],geo,.05)
        def match(k):
            distances={i:math.dist(v,anchor.coords[k]) for i,v in geo.coords.items()};node=min(distances,key=distances.get)
            # Cross-device correspondence only: interface coordinates differ by
            # ~9e-10 um. This is not the same-mesh export identity gate.
            assert distances[node]<1e-8;return node
        focus={name:tuple(match(k) for k in ids) for name,ids in [('channel',(320,324)),('drain_contact',(1091,1092))]}
        interface=match(338)
        state={role:d.ordered(Path(c[role]),N) for role in ('strict','mapped')}
        ss={role:a.read(root/role/'functional.status.json') for role in state}
        rr={role:d.array(d.ordered(root/role/'residual.csv',N),('psi_residual','phin_residual','phip_residual')) for role in state}
        ee={role:a.rows(root/role/'edges.csv') for role in state}
        mesh=a.read(Path(a.read(Path(c['config']))['mesh_file']))
        for role in state:
            cc=d.fixed.sum_contacts(ee[role],mesh,np.array([float(e['electron_particle_line_flux_per_m_s']) for e in ee[role]],dtype=np.longdouble),np.array([float(e['hole_particle_line_flux_per_m_s']) for e in ee[role]],dtype=np.longdouble))
            error=abs(cc['drain']/ss[role]['current_A_per_um']-1);assert error<=1e-8
            checks.append(dict(key=c['key'],check='edge_terminal_'+role,error=error,passed=True))
        strictI=ss['strict']['current_A_per_um'];mappedI=ss['mapped']['current_A_per_um'];nativeI=c['native_Id_A_per_um']
        ref=a.read(Path(c['base'])/'result.json')['Id_A_per_um'];assert abs(strictI/ref-1)<=1e-8
        adj=a.read(root/'adjoint.status.json');assert adj['adjoint_relative_residual']<=1e-10
        weights=d.array(d.ordered(root/'adjoint.csv',N),('lambda_poisson','lambda_electron','lambda_hole'))
        deltaF=rr['mapped']-rr['strict'];prediction=d.project(weights,deltaF)
        cases.append(dict(key=c['key'],device=c['device'],vd=c['vd'],vg=c['vg'],strict_Id_A_per_um=strictI,native_Id_A_per_um=nativeI,
            mapped_functional_A_per_um=mappedI,mapped_target_relative=mappedI/nativeI-1,
            screening_state_prediction_A_per_um=prediction,finite_gap_remainder_A_per_um=strictI-mappedI-prediction,
            relative_error=strictI/nativeI-1,adjoint_relative_residual=adj['adjoint_relative_residual'],
            screening_only=True,channel_node0=focus['channel'][0],channel_node1=focus['channel'][1],interface_node=interface))
        delta=d.array(state['strict'],('psi','electrons_m3','holes_m3'))-d.array(state['mapped'],('psi','electrons_m3','holes_m3'))
        dr=(rr['strict'][0]-rr['mapped'][0])/factor
        parts={}
        for carrier,source in [('electron',Q*delta[1]*geo.volumes['all_cell']),('hole',-Q*delta[2]*geo.volumes['all_cell'])]:
            for name in ('source','channel','drain','substrate'):
                parts[carrier+'_'+name]=geo.solve('legacy',source*masks[name])
        boundary=np.zeros(N);boundary[geo.contact_nodes]=delta[0,geo.contact_nodes]
        boundary+=geo.solve('legacy',geo.matrices['legacy']@boundary);parts['boundary']=boundary
        parts['operator_residual']=geo.solve('legacy',-dr)
        predicted=sum(parts.values());err=np.linalg.norm((predicted-delta[0])[geo.free])/np.linalg.norm(delta[0,geo.free]);assert err<=1e-5
        checks.append(dict(key=c['key'],check='poisson_conditional_reconstruction',error=float(err),passed=True))
        selection=(*focus['channel'],interface)
        for name,value in parts.items():
            for node in selection:components.append(dict(key=c['key'],component=name,node=node,potential_V=float(value[node])))
        for node in sorted(set(selection+focus['drain_contact'])):
            u,v=state['strict'][node],state['mapped'][node]
            nodes.append(dict(key=c['key'],node=node,x_um=geo.coords[node][0],y_um=geo.coords[node][1],is_interface=bool(masks['interface'][node]),
                strict_minus_native_psi_V=float(u['psi'])-float(v['psi']),strict_minus_native_phi_n_V=float(qf(u)-qf(v)),
                strict_minus_native_phi_p_V=float(qf(u,'hole')-qf(v,'hole')),electron_density_dex=math.log10(float(u['electrons_m3'])/float(v['electrons_m3'])),
                reconstructed_delta_psi_V=float(predicted[node]),operator_residual_potential_V=float(parts['operator_residual'][node])))
        for name,(i,j) in focus.items():
            es=next(e for e in ee['strict'] if (int(e['node0']),int(e['node1']))==(i,j));em=next(e for e in ee['mapped'] if int(e['edge_id'])==int(es['edge_id']))
            assert float(es['electron_mobility_m2_V_s'])==float(em['electron_mobility_m2_V_s'])
            vals={}
            for role,e in [('strict',es),('mapped',em)]:
                u,v=state[role][i],state[role][j];kap,z=f.secant(u['psi'],v['psi'],qf(u),qf(v),e['ni0_m3'],e['ni1_m3'],d.fixed.upstream.VT,float(e['couple_m'])/float(e['length_m']))
                vals[role]=(kap,z)
            ks,zs=vals['strict'];km,zm=vals['mapped']
            edges.append(dict(key=c['key'],focus=name,edge=int(es['edge_id']),node0=i,node1=j,mobility_cm2_Vs=1e4*float(es['electron_mobility_m2_V_s']),
                strict_qf_drop_V=float(qf(state['strict'][j])-qf(state['strict'][i])),native_qf_drop_V=float(qf(state['mapped'][j])-qf(state['mapped'][i])),
                strict_over_mapped_SG_kappa=ks/km,strict_over_mapped_log_imbalance=zs/zm,
                strict_over_mapped_flux=float(es['electron_flux'])/float(em['electron_flux']),
                native_edge_operator_exported=False))
    for name,values in [('case_ledger',cases),('poisson_components',components),('local_nodes',nodes),('local_sg',edges),('checks',checks)]:a.write_csv(OUT/(name+'.csv'),values)
    a.write(OUT/'analysis_scope.json',dict(status='completed_read_only_ledger',input_hashes={a.rel(Path(__file__).resolve()):a.sha(Path(__file__).resolve())},
        finite_gap_weights='Fresh but only screened here; source calibration is recorded separately. No arbitrary finite-state-gap linearization claim.',
        poisson='The residual term is retained. Exact reconstruction is bookkeeping, not proof of native Poisson equivalence or of charge causing the discrepancy.',
        sg='Mapped-native potentials evaluated with common Vela ni/Vt/geometry; actual native edge SG operator has not been exported.',new_nonlinear_solves=0))
    print('Eight local/Poisson ledgers, 16 conservative port identities passed',flush=True)


if __name__=='__main__':main()
