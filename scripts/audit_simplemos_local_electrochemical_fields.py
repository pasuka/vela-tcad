"""Native field identities and common-geometry SG secant proxies near two edges."""
from decimal import Decimal, localcontext
from pathlib import Path
import math
import numpy as np
import decompose_simplemos_calibrated_transport as d
import validate_simplemos_local_conservative_flux as out


def scalar(path):
    return {int(r['node_id']):r['component0'] for r in d.a.rows(path)}


def secant(psi0,psi1,phi0,phi1,ni0,ni1,vt,geometry):
    with localcontext() as context:
        context.prec=60
        p0,p1,f0,f1,n0,n1,t,g=map(Decimal,map(str,(psi0,psi1,phi0,phi1,ni0,ni1,vt,geometry)))
        eta=(p1-p0)/t
        def b(x):return Decimal(1) if x==0 else x/(x.exp()-1)
        left=n0*((p0-f0)/t).exp()*b(-eta)
        right=n1*((p1-f1)/t).exp()*b(eta)
        z=(left/right).ln()
        mean=left if z==0 else (left-right)/z
        return float(t*g*mean),float(z)


def main():
    d.a.verify(d.FREEZE);d.a.verify(d.OUT/'evidence.json')
    source=d.a.rows(d.fixed.upstream.OUT/'m80b_source_identity.csv')
    edges_out,nodes_out=[],[]
    for c in d.a.read(d.CONTRACT)['cases']:
        if c['device']!='n23':continue
        key=c['case'];folder=Path(c['export'])/'fields'
        fields={name:scalar(folder/(name+'_region0.csv')) for name in ('ElectrostaticPotential','eQuasiFermiPotential','eDensity','EffectiveIntrinsicDensity')}
        state={int(r['node_id']):r for r in d.a.rows(d.m.spatial.precision.LOCAL/key/'state.csv')} if hasattr(d,'m') else {int(r['node_id']):r for r in d.a.rows(d.matrix.spatial.precision.LOCAL/key/'state.csv')}
        geo=d.matrix.spatial.m73.Geometry('n23')
        masks,_,xy=d.matrix.spatial.old.m78.supports('n23',geo,.05)
        interface=np.flatnonzero(masks['interface'])
        edge_rows=d.a.rows(d.LOCAL/key/'strict/edges.csv')
        si_vt=next(r['source_inferred_thermal_voltage_V'] for r in source if r['device']=='n23' and float(r['drain_voltage_V'])==float(c['vd']))
        vela_vt=d.fixed.upstream.VT
        for eid in (2392,965):
            edge=next(r for r in edge_rows if int(r['edge_id'])==eid)
            i,j=int(edge['node0']),int(edge['node1'])
            selected={i,j}
            for r in edge_rows:
                if int(r['node0']) in (i,j) or int(r['node1']) in (i,j):selected.update((int(r['node0']),int(r['node1'])))
            centre=(xy[i]+xy[j])/2
            nearest=int(interface[np.argmin(np.linalg.norm(xy[interface]-centre,axis=1))])
            selected.add(nearest)
            native_phi=[fields['eQuasiFermiPotential'][k] for k in (i,j)]
            native_psi=[fields['ElectrostaticPotential'][k] for k in (i,j)]
            vela_phi=[str(Decimal(state[k]['electron_qf_reference_V'])+Decimal(state[k]['electron_qf_increment_V'])) for k in (i,j)]
            vela_psi=[state[k]['psi'] for k in (i,j)]
            ni=[str(Decimal(fields['EffectiveIntrinsicDensity'][k])*Decimal('1e6')) for k in (i,j)]
            geometry=float(edge['couple_m'])/float(edge['length_m'])
            ks,zs=secant(*native_psi,*native_phi,*ni,si_vt,geometry)
            kv,zv=secant(*vela_psi,*vela_phi,*ni,vela_vt,geometry)
            length=float(edge['length_m'])
            ns=float(Decimal(native_phi[1])-Decimal(native_phi[0]));nv=float(Decimal(vela_phi[1])-Decimal(vela_phi[0]))
            sensitivity=(math.ulp(float(native_phi[0]))+math.ulp(float(native_phi[1])))/max(abs(ns),1e-300)
            edges_out.append({'case':key,'edge_id':eid,'node0':i,'node1':j,
                'native_qf_drop_V':ns,'vela_qf_drop_V':nv,'native_qf_gradient_V_m':ns/length,'vela_qf_gradient_V_m':nv/length,
                'native_double_ulp_fraction_of_drop':sensitivity,'native_SG_state_conductance_proxy':ks,'vela_SG_state_conductance_proxy':kv,
                'native_over_vela_conductance':ks/kv,'native_log_imbalance':zs,'vela_log_imbalance':zv,
                'native_psi_drop_V':float(Decimal(native_psi[1])-Decimal(native_psi[0])),
                'vela_psi_drop_V':float(Decimal(vela_psi[1])-Decimal(vela_psi[0])),
                'nearest_SiSiO2_node':nearest,'distance_to_interface_um':float(np.linalg.norm(xy[nearest]-centre)),
                'native_edge_conductance_exported':False})
            for node in sorted(selected):
                if node not in fields['eDensity']:continue
                npotential=Decimal(fields['ElectrostaticPotential'][node]);nphi=Decimal(fields['eQuasiFermiPotential'][node])
                vphi=Decimal(state[node]['electron_qf_reference_V'])+Decimal(state[node]['electron_qf_increment_V'])
                nodes_out.append({'case':key,'focus_edge':eid,'node_id':node,'x_um':xy[node,0],'y_um':xy[node,1],
                    'is_SiSiO2_interface':bool(masks['interface'][node]),'native_psi_V':float(npotential),'vela_psi_V':float(state[node]['psi']),
                    'native_minus_vela_psi_V':float(npotential-Decimal(state[node]['psi'])),
                    'native_minus_vela_eQF_V':float(nphi-vphi),
                    'native_density_cm3':float(fields['eDensity'][node]),'vela_density_cm3':float(state[node]['electrons_m3'])*1e-6})
    d.a.write_csv(out.OUT/'local_fields.csv',edges_out)
    d.a.write_csv(out.OUT/'local_nodes.csv',nodes_out)
    d.a.write(out.OUT/'field_audit_evidence.json',{'input_hashes':{d.a.rel(Path(__file__).resolve()):d.a.sha(Path(__file__).resolve())},
        'definition':'Common Vela edge geometry and Boltzmann SG state secant using each solver thermal convention, excluding mobility. Native proxy is not an exported native edge operator.',
        'conductance_units':'particle line flux per mobility (m2/V/s) per dimensionless log imbalance; geometry couple/length dimensionless, density m^-3, thermal voltage V',
        'new_solver_calls':0})
    print('Native/local field audit complete',flush=True)


if __name__=='__main__':main()
