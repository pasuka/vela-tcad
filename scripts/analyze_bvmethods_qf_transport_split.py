#!/usr/bin/env python3
"""Separate inverse-statistics/contact changes from frozen-population transport.

The counterfactuals retain populations, mobility, degeneracy secants and Vela
sources, and shift effective band drift with QF. They are diagnostic operators,
not admissible device states or a certificate of native SG support. Every full
RHS uses the baseline Vela chained diagonal and the native original diagonal;
no new-state Jacobian is normalized using an old population transform.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

import numpy as np
from compare_sentaurus_debug_matrix import scalar_fields, iterations

VT = .025851999786435535
NREF = 14824351000.


def divergence(nodes, left, right, flux):
    """Oriented flux contributes + at left and - at right, including contacts."""
    left, right = np.asarray(left, dtype=int), np.asarray(right, dtype=int)
    flux = np.asarray(flux, dtype=float)
    if (left.shape != right.shape or left.shape != flux.shape or
            np.any(left < 0) or np.any(right < 0) or np.any(left >= nodes) or
            np.any(right >= nodes) or np.any(~np.isfinite(flux))):
        raise ValueError("invalid edge incidence or flux")
    out = np.zeros(nodes)
    np.add.at(out, left, flux)
    np.add.at(out, right, -flux)
    return out


def qf_components(psi, ni, dos, native_band, native_qf, eta, native_eta, native_vt, sign):
    """Exact additive split; DOS, voltage convention, anchor and native approximation."""
    if sign not in (-1, 1) or np.any(np.asarray(ni) <= 0) or np.any(np.asarray(dos) <= 0):
        raise ValueError("invalid carrier sign, intrinsic density or DOS")
    anchor = psi - sign * VT * np.log(ni / dos) + native_band
    dos_shift = sign * VT * (eta - native_eta)
    thermal = sign * (VT - native_vt) * native_eta
    native_inverse = -native_band + sign * native_vt * native_eta - native_qf
    return np.array([anchor, dos_shift, thermal, native_inverse])


def residual_partition(baseline, contact, interior, all_qf, native):
    """A telescoping vector identity; L2 norms of terms do not add."""
    arrays = [np.asarray(x, dtype=float) for x in (baseline, contact, interior, all_qf, native)]
    if any(x.shape != arrays[0].shape or np.any(~np.isfinite(x)) for x in arrays):
        raise ValueError("incompatible/nonfinite residual vectors")
    a,b,c,d,n = arrays
    return np.array([a-b, b-c, c-d, d-n])


def stats(x):
    x = np.asarray(x)
    return dict(zip(('median','p95','max'), map(float,np.quantile(np.abs(x),[.5,.95,1]))))


def relative(x, ref):
    den = float(np.linalg.norm(ref))
    return float(np.linalg.norm(x)/den) if den else None


def read(path):
    return np.genfromtxt(path, delimiter=',', names=True, encoding='utf-8')


def helper(root, lines):
    env=os.environ.copy()
    if os.name=='nt': env['PATH']='D:/msys64/ucrt64/bin;'+env['PATH']
    exe=root/('fermi_sg.exe' if os.name=='nt' else 'fermi_sg')
    p=subprocess.run([str(exe)], input='\n'.join(lines)+'\n', text=True,
                     capture_output=True, check=True, env=env, timeout=60)
    return np.loadtxt(p.stdout.splitlines())


def inverse(root, ratios):
    return np.atleast_2d(helper(root, ['r '+format(x,'.17g') for x in ratios]))[:,0]


def analyze(root):
    mesh=json.loads((root/'inputs/mesh.json').read_text()); nodes=len(mesh['nodes'])
    contacts=np.array(sorted({i for c in mesh['contacts'] if c['name']!='gate' for i in c['node_ids']}))
    contact_mask=np.zeros(nodes,dtype=bool); contact_mask[contacts]=True
    doping=read(root/'inputs/doping.csv')
    c0=max(abs(doping['donors_cm3']-doping['acceptors_cm3']))
    materials=json.loads((root/'inputs/materials.json').read_text())['materials']
    mu0=max(m.get(k,0.) for m in materials for k in ['electron_mobility_cm2_per_V_s','hole_mobility_cm2_per_V_s'])
    # SI n/mu exports -> internal cm^-3/cm^2; mesh is um, field factor is 1e4.
    scale=100./(c0*mu0*VT)
    native_sources=None
    source_root=root/'inputs/native_source_controls'
    if source_root.exists():
        native_sources={a:next(iterations((source_root/a/'debug').open())) for a in
                        ['default','no_srh','no_btbt','no_avalanche']}
        primary=native_sources['default']['variables']
        for arm,item in native_sources.items():
            if item['variables'] != primary:
                raise ValueError('native source controls have different primary states: '+arm)
    report={'method':__doc__, 'continuity_SI_kernel_to_internal_factor':scale,
            'qualified_native_full_residual':False, 'snapshots':{}}
    for label in ['zero','first','second']:
        inp=root/'inputs'/label; fields=inp/'native_fields'; mapping=np.load(inp/'mapping.npz')
        nv=float(mapping['native_potential_scale']); common=read(inp/'rows.csv')
        populations={a:read(root/a/label/'populations/rows.csv') for a in ['baseline','all_population','native_qf']}
        terms={a:read(root/a/label/'terms/rows.csv') for a in populations}
        rhs={a:read(root/a/label/'matrix/rhs.csv') for a in populations}
        edges={a:read(root/a/label/'sg/rows.csv') for a in populations}
        rec={'qf_mapping':{}, 'carriers':{}, 'full_rhs':{}, 'internal_closure':{}}
        # All common rows share fixed left scales across every diagnostic arm.
        rb=common['block'].astype(int); rn=common['node'].astype(int)
        diag=common['vela_chained_diagonal']; native_rhs=common['native_rhs_over_diagonal']
        full={a:rhs[a]['rhs'][rb*nodes+rn]/diag for a in rhs}
        for b,name in enumerate(['psi','electron','hole']):
            mask=rb==b
            rec['full_rhs'][name]={a+'_difference_L2':float(np.linalg.norm(v[mask]-native_rhs[mask])) for a,v in full.items()}
        for bi,(carrier,qf_name,pop_name,dos_name,band_name,sign) in enumerate([
            ('electron','eQuasiFermiPotential','eDensity','eEffectiveStateDensity','ConductionBandEnergy',-1),
            ('hole','hQuasiFermiPotential','hDensity','hEffectiveStateDensity','ValenceBandEnergy',1)]):
            b=bi+1; ids=rn[rb==b]; native=scalar_fields(fields,qf_name)
            native_pop=scalar_fields(fields,pop_name); native_dos=scalar_fields(fields,dos_name)
            native_band=scalar_fields(fields,band_name)
            p=populations['all_population']; target=p['input_'+carrier+'_density_m3'][ids]
            # Use debug-matched target population in both inverse maps, never a post-update frame.
            eta=inverse(root,target/p[('Nc_m3' if bi==0 else 'Nv_m3')][ids])
            native_eta=inverse(root,target/1e6/np.array([native_dos[i] for i in ids]))
            native_qf=np.array([native[i] for i in ids]); band=np.array([native_band[i] for i in ids])
            comp=qf_components(p['psi'][ids],p['ni_eff_m3'][ids],p[('Nc_m3' if bi==0 else 'Nv_m3')][ids],
                               band,native_qf,eta,native_eta,nv,sign)
            diff=p[carrier+'_required_qf_V'][ids]-native_qf
            rec['qf_mapping'][carrier]={'difference_V':stats(diff),
                'components_V':{k:stats(x) for k,x in zip(['band_anchor','DOS','thermal_voltage','native_inverse_approximation'],comp)},
                'component_closure_max_V':float(max(abs(comp.sum(axis=0)-diff))),
                'population_errors':{a:{'free_max_relative':float(max(abs(v['reconstructed_'+carrier+'_density_m3'][ids]/target-1))),
                    'contact_max_relative':float(max(abs(v['reconstructed_'+carrier+'_density_m3'][contacts]/v['input_'+carrier+'_density_m3'][contacts]-1)))} for a,v in populations.items()},
                'all_population_contact_QF_violation_max_V':float(max(abs(p[('phin' if bi==0 else 'phip')][contacts]-populations['baseline'][('phin' if bi==0 else 'phip')][contacts])))}
            s=edges['all_population']; left=s['node0'].astype(int); right=s['node1'].astype(int)
            active=(s[carrier+'_mobility_m2_V_s']>0)&(s[carrier+'_density0_m3']>0)&(s[carrier+'_density1_m3']>0)
            ii=left[active]; jj=right[active]; e=s[active]
            qcol='phin' if bi==0 else 'phip'
            old0=e[qcol+'0_V']; old1=e[qcol+'1_V']
            if bi==0:
                eta0=e['electron_eta0']; eta1=e['electron_eta1']; drift=e['electron_drift_potential_V']
            else:
                eta0=(old0-e['psi0_V'])/VT+np.log(e['ni0_m3']/p['Nv_m3'][ii])
                eta1=(old1-e['psi1_V'])/VT+np.log(e['ni1_m3']/p['Nv_m3'][jj])
                drift=e['psi1_V']-e['psi0_V']+VT*np.log((e['ni0_m3']/p['Nv_m3'][ii])/(e['ni1_m3']/p['Nv_m3'][jj]))
            coef=e[carrier+'_mobility_m2_V_s']*VT*e['couple_m']/e['length_m']*scale
            native0=np.array([native[i] for i in ii]); native1=np.array([native[i] for i in jj])
            counter={}; edge_table=[]
            for arm in ['replay','native_interior','native_all']:
                if arm=='replay': q0,q1=old0,old1
                elif arm=='native_interior':
                    q0=np.where(contact_mask[ii],old0,native0); q1=np.where(contact_mask[jj],old1,native1)
                else: q0,q1=native0,native1
                d=drift+(q1-q0)-(old1-old0)
                args=zip(e[carrier+'_density0_m3'],e[carrier+'_density1_m3'],eta0,eta1,d,q0,q1,np.full(len(e),VT),coef)
                flux=helper(root,[('e' if bi==0 else 'h')+' '+' '.join(format(float(v),'.17g') for v in row) for row in args])
                allflux=np.zeros(len(s)); allflux[active]=flux
                counter[arm]=divergence(nodes,left,right,allflux)
                edge_table.append(allflux)
            out=root/'counterfactuals'; out.mkdir(exist_ok=True)
            with (out/(label+'_'+carrier+'_edges.csv')).open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['edge_id','node0','node1','contact_adjacent','production_flux','replay_flux','native_interior_flux','native_all_flux'])
                w.writerows(zip(s['edge_id'].astype(int),left,right,contact_mask[left]|contact_mask[right],s[carrier+'_flux'],*edge_table))
            t=terms['all_population']; af=t[carrier+'_flux']; ar=t[carrier+'_residual']
            # Never recover a weak source by subtracting two large transport sums.
            sources=sum(t[carrier+'_'+k] for k in ['recombination','impact','gauge','boundary'])
            cf={a:sources+flux for a,flux in counter.items()}
            mask=rb==b; d=diag[mask]; nr=native_rhs[mask]
            bb=full['baseline'][mask]; pp=full['all_population'][mask]
            ci=-cf['native_interior'][ids]/d; ca=-cf['native_all'][ids]/d
            parts=residual_partition(bb,pp,ci,ca,nr)
            replay=counter['replay'][ids]; flux_diff=replay-af[ids]
            rec['carriers'][carrier]={'positive_mobility_edges':int(active.sum()),
                'kernel_replay_relative_L2':relative(edge_table[0]-s[carrier+'_flux'],s[carrier+'_flux']),
                'divergence_replay_relative_L2':relative(flux_diff,af[ids]),
                'partition_closure_relative_L2':relative(parts.sum(axis=0)-(bb-nr),bb-nr),
                'RHS_difference_L2':{a:float(np.linalg.norm(v-nr)) for a,v in [('baseline',bb),('all_population',pp),('frozen_native_interior',ci),('frozen_native_all',ca)]},
                'partition_L2':{k:float(np.linalg.norm(v)) for k,v in zip(['contact_population','interior_QF_transport','contact_QF_transport','unresolved_transport_sources'],parts)},
                'frozen_sources_L2_over_baseline_diag':{k:float(np.linalg.norm(t[carrier+'_'+k][ids]/d)) for k in ['recombination','impact','gauge','boundary']},
                'population_relative_RHS_error':{a:stats((v-nr)/(np.array([native_pop[i] for i in ids])/NREF)) for a,v in [('baseline',bb),('all_population',pp),('frozen_native_all',ca),('native_qf_state',full['native_qf'][mask])]}}
            if rec['carriers'][carrier]['kernel_replay_relative_L2'] > 1e-10:
                raise ValueError('frozen SG replay failed: '+label+' '+carrier)
            if label=='first' and native_sources is not None:
                base=native_sources['default']; keys=[(bi+2,int(i)) for i in ids]
                nd=common['native_diagonal'][mask]
                native_fixed=np.array([base['rhs1'][k] for k in keys])/nd
                native_parts={a:np.array([base['rhs1'][k]-native_sources[a]['rhs1'][k] for k in keys])/nd
                              for a in ['no_srh','no_btbt','no_avalanche']}
                native_transport=native_fixed-sum(native_parts.values())
                vela_transport=-counter['native_all'][ids]/d
                vela_recomb=-t[carrier+'_recombination'][ids]/d
                vela_impact=-t[carrier+'_impact'][ids]/d
                native_recomb=native_parts['no_srh']+native_parts['no_btbt']
                native_impact=native_parts['no_avalanche']
                original=np.load(inp/'native_vectors.npz')
                ov={(int(bb),int(nn)):xx for bb,nn,xx in zip(original['block'],original['node'],original['variables'])}
                primary_relative=max(abs(base['variables'][k]/ov[k]-1) for k in keys)
                rec['carriers'][carrier]['native_source_control']={
                    'scope':'first pre-update Save round trip; identical primary variables across four native arms; finite subtraction resolution retained; frozen Vela QF transport and sources do not form a re-solved state',
                    'save_population_roundtrip_max_relative_error':primary_relative,
                    'save_rhs_vs_original_normalized_difference_L2':float(np.linalg.norm(native_fixed-nr)),
                    'native_transport_L2':float(np.linalg.norm(native_transport)),
                    'frozen_transport_difference_L2':float(np.linalg.norm(vela_transport-native_transport)),
                    'recombination_difference_L2':float(np.linalg.norm(vela_recomb-native_recomb)),
                    'avalanche_difference_L2':float(np.linalg.norm(vela_impact-native_impact)),
                    'native_avalanche_L2':float(np.linalg.norm(native_impact)),
                    'frozen_vela_avalanche_L2':float(np.linalg.norm(vela_impact)),
                    'full_difference_L2_against_save':float(np.linalg.norm(ca-native_fixed)),
                    'full_partition_closure_relative_L2':relative((vela_transport-native_transport)+
                        (vela_recomb-native_recomb)+(vela_impact-native_impact)-(ca-native_fixed),ca-native_fixed),
                    'source_subtraction_unresolved_rows':{a:int(sum(abs(base['rhs1'][k]-native_sources[a]['rhs1'][k]) <=
                        64*max(abs(np.spacing(base['rhs1'][k])),abs(np.spacing(native_sources[a]['rhs1'][k]))) for k in keys)
                        ) for a in native_parts}}
                with (out/(label+'_'+carrier+'_native_source_rows.csv')).open('w',newline='') as f:
                    w=csv.writer(f);w.writerow(['node','native_fixed_rhs','native_transport_rhs','native_SRH_rhs','native_BTBT_rhs','native_avalanche_rhs',
                        'frozen_vela_transport_rhs','frozen_vela_recombination_rhs','frozen_vela_avalanche_rhs'])
                    w.writerows(zip(ids,native_fixed,native_transport,native_parts['no_srh'],native_parts['no_btbt'],native_impact,
                                    vela_transport,vela_recomb,vela_impact))
            # Every common carrier row is retained, including the minority tail.
            with (out/(label+'_'+carrier+'_rows.csv')).open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['node','contact_adjacent','native_rhs_over_diagonal','baseline_rhs_over_fixed_diagonal','all_population_rhs_over_fixed_diagonal','frozen_native_interior_rhs','frozen_native_all_rhs','native_qf_state_rhs',
                    'contact_population','interior_QF_transport','contact_QF_transport','unresolved_transport_sources'])
                adj=np.zeros(nodes,dtype=bool)
                hit=contact_mask[left]|contact_mask[right];adj[left[hit]]=True;adj[right[hit]]=True
                w.writerows(zip(ids,adj[ids],nr,bb,pp,ci,ca,full['native_qf'][mask],*parts))
            rec['internal_closure'][carrier]={}
            for a in populations:
                tt=terms[a]; es=edges[a]
                div=divergence(nodes,es['node0'],es['node1'],es[carrier+'_flux'])
                term_sum=sum(tt[carrier+'_'+k] for k in ['flux','recombination','impact','gauge','boundary'])
                rec['internal_closure'][carrier][a]={'edge_to_flux_relative_L2':relative(div[ids]-tt[carrier+'_flux'][ids],tt[carrier+'_flux'][ids]),
                    'term_to_residual_relative_L2':relative(term_sum[ids]-tt[carrier+'_residual'][ids],tt[carrier+'_residual'][ids]),
                    'residual_to_full_RHS_relative_L2':relative(tt[carrier+'_residual'][ids]+rhs[a]['rhs'][b*nodes+ids],tt[carrier+'_residual'][ids])}
        geometry_path=root/'inputs/native_charge_geometry_rows.csv'
        if geometry_path.exists():
            geometry=read(geometry_path)
            ratio=np.ones(nodes)
            gi=geometry['node_id'].astype(int)
            ratio[gi]=geometry['native_semiconductor_volume_um2']/geometry['barycentric_total_um2']
            mask=rb==0;ids=rn[mask];d=diag[mask];p=populations['baseline']
            charge=p['electron_charge']+p['hole_charge']+p['doping_charge']
            raw=p['dielectric_flux']+charge+p['fixed_interface_charge']
            corrected=p['dielectric_flux']+charge*ratio+p['fixed_interface_charge']
            native=-native_rhs[mask]
            rec['poisson_box_counterfactual']={
                'scope':'residual only, supplied populations and native Si box replacement; fixed baseline diagonal; no Jacobian or production geometry change',
                'probe_vs_full_residual_relative_L2':relative(raw[ids]-rhs['baseline']['residual'][ids],rhs['baseline']['residual'][ids]),
                'baseline_difference_L2':float(np.linalg.norm(raw[ids]/d-native)),
                'native_Si_box_difference_L2':float(np.linalg.norm(corrected[ids]/d-native))}
            cell_path=root/'inputs/native_cell_geometry_rows.csv'
            if cell_path.exists():
                cells=np.genfromtxt(cell_path,delimiter=',',names=True,dtype=None,encoding='utf-8')
                left=cells['node_id'].astype(int)
                right=np.array([mesh['triangles'][int(r['cell_id'])]['node_ids'][(int(r['local'])+1)%3] for r in cells])
                eps={m['name']:m['eps_r'] for m in materials}; epsref=max(eps.values())
                conductance=np.array([eps[m] for m in cells['material']])*cells['native_coefficient']/(epsref*VT)
                native_flux=divergence(nodes,left,right,conductance*(p['psi'][left]-p['psi'][right]))
                # Calibrate one Poisson charge MATRIX coefficient, never fit a residual.
                # It determines q/epsilon in native conventions, not q and epsilon separately.
                fit=geometry[geometry['node_id']==243]
                if len(fit)!=1: raise ValueError('missing independent bulk charge coefficient')
                coeff=(geometry['native_normalized_charge_coefficient']*
                       geometry['native_laplace_diag_relative_epsilon']/geometry['native_semiconductor_volume_um2'])
                k=float(fit['native_normalized_charge_coefficient'][0]*fit['native_laplace_diag_relative_epsilon'][0]/fit['native_semiconductor_volume_um2'][0])
                nominal=1.602176634e-19*NREF*1e-6/(8.8541878128e-12*nv)
                native_scale=k/nominal
                native_geometry=native_flux+charge*ratio+p['fixed_interface_charge']
                native_operator=native_flux+charge*ratio*native_scale+p['fixed_interface_charge']
                rec['poisson_box_counterfactual'].update(
                    native_cell_coefficients_and_box_difference_L2=float(np.linalg.norm(native_geometry[ids]/d-native)),
                    native_matrix_charge_scale_over_vela=native_scale,
                    matrix_charge_scale_fit_node=243,
                    heldout_charge_matrix_coefficient_max_relative_error=float(max(abs(coeff/k-1))),
                    native_cell_box_and_matrix_scale_difference_L2=float(np.linalg.norm(native_operator[ids]/d-native)),
                    native_cell_flux_change_L2=float(np.linalg.norm((native_flux-p['dielectric_flux'])[ids]/d)))
                with (out/(label+'_poisson_native_operator_rows.csv')).open('w',newline='') as f:
                    w=csv.writer(f);w.writerow(['node','native_residual_over_diagonal','native_cell_box_residual_over_fixed_diagonal','native_cell_box_matrix_scale_residual_over_fixed_diagonal'])
                    w.writerows(zip(ids,native,native_geometry[ids]/d,native_operator[ids]/d))
            out=root/'counterfactuals'
            with (out/(label+'_poisson_rows.csv')).open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['node','native_residual_over_diagonal','baseline_residual_over_fixed_diagonal','native_box_residual_over_fixed_diagonal','native_box_over_global_bary_volume'])
                w.writerows(zip(ids,native,raw[ids]/d,corrected[ids]/d,ratio[ids]))
        report['snapshots'][label]=rec
    (root/'analysis.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('root',type=Path)
    args=parser.parse_args();result=analyze(args.root.resolve())
    print(json.dumps({s:{c:r['RHS_difference_L2'] for c,r in a['carriers'].items()} for s,a in result['snapshots'].items()},indent=2))
