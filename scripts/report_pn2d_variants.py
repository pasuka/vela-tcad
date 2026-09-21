#!/usr/bin/env python3
"""Render actual PN2D campaign evidence; absent experiments remain visibly absent."""
from __future__ import annotations
import argparse
import csv
import math
import shutil
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from pn2d_variants import DEFAULT_ROOT, SPEC, FIXTURE, read_json, write_json, tag
from compare_pn2d_variants import node_field, rows
from analyze_pn2d_sections import analyze as analyze_nodal_sections
from compare_genius_bjt_transport_fields import read_vtk_point_data


def report(root):
    result=read_json(root/'results.json')
    out=root/'report'; out.mkdir(exist_ok=True)
    coverage=[]
    for entry in read_json(root/'campaign.json')['entries']:
        path=root/entry['directory']/'inputs/mesh.json'
        if not path.exists(): continue
        mesh=read_json(path)
        for contact in mesh['contacts']:
            name=contact['name'];boundary_x=0 if name=='Anode' else 2
            fraction=entry['parameters']['anode_fraction'] if name=='Anode' else 1
            low,high=.25*(1-fraction),.25*(1+fraction)
            expected={n['id'] for n in mesh['nodes'] if abs(n['x']-boundary_x)<1e-10 and low-1e-10<=n['y']<=high+1e-10}
            actual=set(contact['node_ids'])
            if expected!=actual: raise ValueError(f'Incomplete or excessive electrode node coverage: {entry["directory"]}/{name}')
            coverage.append({'id':entry['id'],'mesh':entry['mesh'],'contact':name,'expected_nodes':len(expected),'actual_nodes':len(actual),'status':'pass'})
    write_json(root/'boundary_coverage.json',{'description':'Exact electrode node sets match the requested intervals; remaining left-boundary nodes have no electrode','results':coverage})
    nodal_sections=analyze_nodal_sections(root)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    available=[r for r in result['results'] if 'metrics' in r['terminal']]
    grids=list(read_json(SPEC)['meshes'])
    rank={m:i for i,m in enumerate(grids)}
    refinement_pair=result['refinement']['additional_pair']
    mesh_sets=[('M0/M1',result['mesh']),*result.get('historical_refinement_mesh',{}).items(),
               ('/'.join(refinement_pair),result.get('refined_mesh',{})),
               ('/'.join(result['refinement']['contact_edge_control']['pair']),result.get('contact_edge_mesh',{}))]
    preferred={}
    contact_states=[]
    for r in available:
        if r['id'] not in preferred or rank[r['mesh']]>rank[preferred[r['id']]['mesh']]: preferred[r['id']]=r
    if available:
        fig,axes=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
        for r in available:
            if r is not preferred[r['id']]: continue
            for j,branch in enumerate(['forward','reverse']):
                data=[p for p in r['terminal']['metrics'] if p['contact']=='Anode' and p['component']=='total' and p['branch']==branch and p['bias_V']!=0]
                x=[p['bias_V'] for p in data]
                line,=axes[0,j].semilogy(x,[abs(p['reference_A_per_um']) for p in data],label=f'{r["id"]}/{r["mesh"]} Sentaurus')
                axes[0,j].semilogy(x,[abs(p['candidate_A_per_um']) for p in data],'--',color=line.get_color())
                axes[1,j].plot(x,[100*p['relative_error'] for p in data],label=f'{r["id"]}/{r["mesh"]}')
                axes[0,j].set(title=branch,ylabel='|Anode current| [A/um]')
                axes[1,j].set(xlabel='Anode voltage [V]',ylabel='Absolute relative error [%]')
                axes[0,j].legend(fontsize=7); axes[1,j].grid(alpha=.2)
        fig.suptitle('Solid: Sentaurus; dashed: independently solved Vela')
        fig.savefig(out/'terminal_curves.png',dpi=160); plt.close(fig)
    plot_terminal_components(preferred,out)
    for r in available:
        if r is not preferred[r['id']] or 'points' not in r['spatial']: continue
        d=root/r['id']/r['mesh']; mesh=read_json(d/'inputs/mesh.json')
        nodes=mesh['nodes']; ids=[n['id'] for n in nodes]
        index={v:i for i,v in enumerate(ids)}
        doping={int(p['node_id']):(float(p['donors_cm3']),float(p['acceptors_cm3'])) for p in rows(d/'inputs/doping.csv')}
        x=np.array([n['x'] for n in nodes]); y=np.array([n['y'] for n in nodes])
        triangles=[[index[i] for i in t['node_ids']] for t in mesh['triangles']]
        triangulation=mtri.Triangulation(x,y,triangles)
        fig,axes=plt.subplots(5,3,figsize=(12,10),constrained_layout=True)
        for row,p in enumerate(r['spatial']['points']):
            bias=p['bias_V']; branch='reverse' if bias<0 else 'forward'
            fields=d/'fields'/f'{branch}_{tag(bias)}'/'fields'
            _,scalars,vectors=read_vtk_point_data(d/p['vtk'])
            ref=np.array(node_field(fields/'ElectrostaticPotential_region0.csv',ids))
            actual=np.array(scalars['Potential'])
            density=np.log10(np.maximum(np.array(scalars['Electrons']),1.0))
            native={name:node_field(fields/f'{field}_region0.csv',ids) for name,field in
                    [('n','eDensity'),('p','hDensity'),('ni','EffectiveIntrinsicDensity'),
                     ('qn','eQuasiFermiPotential'),('qp','hQuasiFermiPotential')]}
            candidate={name:scalars[field] for name,field in [('n','Electrons'),('p','Holes'),('ni','EffectiveIntrinsicDensity'),
                       ('qn','ElectronQuasiFermi'),('qp','HoleQuasiFermi')]}
            for contact in mesh['contacts']:
                selected=[index[i] for i in contact['node_ids']]
                item={'id':r['id'],'mesh':r['mesh'],'bias_V':bias,'contact':contact['name'],'nodes':len(selected),
                      'potential_difference_max_V':max(abs(actual[i]-ref[i]) for i in selected)}
                for tool,values in [('reference',native),('candidate',candidate)]:
                    item[tool+'_mass_action_relative_max']=max(abs(values['n'][i]*values['p'][i]/values['ni'][i]**2-1) for i in selected)
                    item[tool+'_neutrality_relative_max']=max(abs(values['n'][i]-values['p'][i]-(doping[ids[i]][0]-doping[ids[i]][1]))/max(sum(doping[ids[i]]),1) for i in selected)
                    voltage=bias if contact['name']=='Anode' else 0
                    item[tool+'_qf_bias_error_max_V']=max(abs(values[q][i]-voltage) for q in ['qn','qp'] for i in selected)
                    for key in ['n','p','ni']:
                        item[tool+'_'+key+'_min_cm3']=min(values[key][i] for i in selected)
                        item[tool+'_'+key+'_max_cm3']=max(values[key][i] for i in selected)
                for key in ['n','p','ni']:
                    item[key+'_agreement_relative_max']=max(abs(candidate[key][i]-native[key][i])/max(abs(native[key][i]),1e-300) for i in selected)
                contact_states.append(item)
            for col,(values,title) in enumerate([(actual,f'{bias:g} V: Vela potential [V]'),((actual-ref)*1e3,'Vela - Sentaurus [mV]'),(density,'Vela log10(n / cm^-3)')]):
                artist=axes[row,col].tripcolor(triangulation,values,shading='gouraud',cmap='coolwarm' if col==1 else 'viridis')
                axes[row,col].set(title=title,xlabel='x [um]',ylabel='y [um]',aspect='equal')
                fig.colorbar(artist,ax=axes[row,col],shrink=.65)
        fig.savefig(out/f'{r["id"]}_{r["mesh"]}_states.png',dpi=150); plt.close(fig)
        if r['id'].startswith('G'):
            plot_legacy_contact_current(r,d,mesh,ids,index,x,y,triangulation,out)
    write_json(root/'contact_state_audit.json',{'scope':'All electrode nodes, including minority populations below the global carrier mask; formula diagnostics without new acceptance thresholds','records':contact_states})
    fig,axes=plt.subplots(1,3,figsize=(12,3.6),constrained_layout=True)
    for ax,identifier in zip(axes,['P1','P2','P3']):
        candidates=result['model_effects'][identifier]
        key=next((m for m in reversed(grids) if 'metrics' in candidates.get(m,{})),None)
        if key:
            points=[p for p in candidates[key]['metrics'] if p['contact']=='Anode' and p['component']=='total']
            points=sorted(points,key=lambda p:p['bias_V'])
            for tool,style in [('reference','-'),('candidate','--')]:
                ax.plot([p['bias_V'] for p in points],[p[tool+'_delta_dex'] for p in points],style,label=tool)
            ax.set(title=f'{identifier}/{key} relative to P0',xlabel='Anode voltage [V]',ylabel='Delta log10 |I| [dex]')
            ax.legend();ax.grid(alpha=.2)
        else: ax.text(.5,.5,'Not run',ha='center',transform=ax.transAxes)
    fig.savefig(out/'model_effects.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,6),constrained_layout=True)
    for ax,identifier in zip(axes.flat,['P0','P1','P2','P3']):
        points=[p for p in nodal_sections if p['id']==identifier and p['bias_V']==.8 and p['component']=='total']
        if not points: continue
        for tool,label in [('reference','Sentaurus nodal'),('candidate','Vela nodal')]:
            ax.plot([p['x_um'] for p in points],[100*p[tool+'_relative_to_port'] for p in points],'.-',label=label)
        ax.set(title=f'{identifier}/{points[0]["mesh"]}, 0.8 V',xlabel='x cut [um]',ylabel='Nodal integral vs port error [%]')
        ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Reconstructed current sections: diagnostic, not discrete conservative flux')
    fig.savefig(out/'nodal_section_diagnostics.png',dpi=160);plt.close(fig)
    plot_spreading_metrics(preferred,out,root)
    fig,axes=plt.subplots(1,len(mesh_sets),figsize=(5*len(mesh_sets),3.5),constrained_layout=True)
    for identifier in [e['id'] for e in read_json(SPEC)['experiments']]:
        for ax,(title,gates) in zip(axes,mesh_sets,strict=True):
            points=[p for p in gates.get(identifier,{}).get('terminal_metrics',[]) if p['tool']=='reference' and p['contact']=='Anode' and p['component']=='hole' and p['branch']=='forward' and p['bias_V']!=0]
            if points: ax.plot([p['bias_V'] for p in points],[100*p['relative_error'] for p in points],label=identifier)
            ax.set(title=title,xlabel='Forward bias [V]',ylabel='Anode hole-current mesh sensitivity [%]')
            ax.axhline(2,color='black',ls=':',lw=1);ax.grid(alpha=.2)
    for ax in axes:
        if ax.get_legend_handles_labels()[0]: ax.legend()
    fig.savefig(out/'mesh_sensitivity.png',dpi=160);plt.close(fig)
    resolution=[]
    fig,axes=plt.subplots(1,2,figsize=(10,3.6),constrained_layout=True)
    for grid in grids:
        if grid=='original': continue
        d=root/'P0'/grid; fields=d/'fields/forward_p0p20/fields'
        if not (fields/'srhRecombination_region0.csv').exists(): continue
        mesh=read_json(d/'inputs/mesh.json');nodes=mesh['nodes'];ids=[n['id'] for n in nodes]
        srh=node_field(fields/'srhRecombination_region0.csv',ids)
        ni=node_field(fields/'EffectiveIntrinsicDensity_region0.csv',ids)
        doping={int(n['node_id']):n for n in rows(d/'inputs/doping.csv')}
        profile=sorted([{'x_um':n['x'],'SRH_cm3_s':srh[i],'ni_eff_cm3':ni[i],
                         'donors_cm3':float(doping[n['id']]['donors_cm3']),'acceptors_cm3':float(doping[n['id']]['acceptors_cm3'])}
                        for i,n in enumerate(nodes) if abs(n['y']-.25)<1e-10 and abs(n['x']-1)<.02],key=lambda p:p['x_um'])
        if not profile: continue
        resolution.append({'mesh':grid,'bias_V':.2,'y_um':.25,'centerline':profile})
        axes[0].plot([p['x_um'] for p in profile],[p['SRH_cm3_s'] for p in profile],'.-',label=grid)
        axes[1].plot([p['x_um'] for p in profile],[p['ni_eff_cm3'] for p in profile],'.-',label=grid)
    for ax,title,ylabel in [(axes[0],'Resolved samples of the narrow SRH peak','SRH [cm^-3 s^-1]'),(axes[1],'Reported compensated junction plane','Effective ni [cm^-3]')]:
        ax.set(title=title,xlabel='x [um], y=0.25 um, 0.2 V',ylabel=ylabel);ax.legend();ax.grid(alpha=.2)
    write_json(root/'junction_resolution.json',resolution)
    fig.savefig(out/'junction_resolution.png',dpi=160);plt.close(fig)
    lines=['# PN2D incremental validation results','',f'Overall status: **{result["overall_status"]}**. Stage A accepted: **{result["stage_A_accepted"]}**.','',
           'Sentaurus T-2022.03-SP2; Vela MSYS2 UCRT64 GNU 16.2.0 Release; HDF5 enabled; ordinary DD uses Eigen SparseLU/COLAMD (SPQR/UMFPACK also compiled).',
           '','Silicon, 300 K, 2 x 0.5 um, straight junction x=1 um. Boltzmann DD, no avalanche. AreaFactor=1 and 1 um depth; currents in A/um. Anode is swept; Cathode is grounded.',
           '',f'Frozen contract SHA256: `{result["contract_sha256"]}`. See `../contract.json` and `../results.json` for every gate and metric.',
           '',f'Completed paired runs: {sum(r["status"]=="pass" for r in result["results"])}. Required final pairs: {sum(r.get("required_for_final_qualification",False) for r in result["results"])}. The remaining generated entries are optional controls and retain their explicit not-run status.',
           '','| Case | Mesh | Required | Input | Ports | Spatial | Conservative transport | Result |','|---|---|---|---|---|---|---|---|']
    for r in result['results']:
        required='yes' if r.get('required_for_final_qualification',False) else 'no'
        lines.append(f'| {r["id"]} | {r["mesh"]} | {required} | {r["input"]["status"]} | {r["terminal"]["status"]} | {r["spatial"]["status"]} | {r["spatial"].get("transport_source_and_section_status","not_run")} | {r["status"]} |')
    lines+=['','## Mesh and model effects','','| Case | '+' | '.join(label for label,_ in mesh_sets)+' |',
            '|---|'+'---|'*len(mesh_sets)]
    lines += ['| '+i+' | '+' | '.join(gates.get(i,{}).get('status','not_run') for _,gates in mesh_sets)+' |' for i in result['mesh']]
    lines += ['', 'Final qualification requires P0/original and all eight configurations on '+ '/'.join(refinement_pair)+
              ', plus P0/G1/G2 on the endpoint control pair '+ '/'.join(result['refinement']['contact_edge_control']['pair'])+
              ', with the unchanged current and source-integral mesh thresholds. Earlier coarse controls and their failures remain visible. Unrun coarse B/C controls are not required to qualify the final refined pair.']
    lines+=['','## Numerical correction','',
            'The original PN2D IV template uses source-only continuity row scaling. Independent P0 scans stalled near +1.12 mV and -1.08 mV with rejected continuity closure. Tightening tolerances and jumping directly to 0.02 V did not resolve the issue. The campaign uses flux_fraction=1 in left row scaling. It changes conditioning, retaining equations, local eps_row=0.001 and global closure tolerance=0.01. The original template and global defaults are unchanged. Failed probes remain under P0/original/initial_source_scaling.',
            '', 'P3 exposed loss of sub-femtovolt conductive increments in the homogeneous-ni SG path: separate quasi-Fermi exponentials rounded to equal values. The homogeneous functions now use the existing expm1-factorized law at equal ni. A physical small-signal conductance test failed before the fix and passes afterward. The 26 SG tests (221 assertions), 13 Gummel cases, and the six selected PN2D/import/reference regressions passed. Failed P3 scans remain under P3/M0/initial_unstable_constant_ni. Repaired independent solves are recorded separately by executable fingerprint.',
            '', 'D1/D2 initially failed low-current port gates with the native double-precision reference. Their input, spatial and conservative checks passed. The D1/M2 EP128 diagnostic reduces native equilibrium current from 4.3e-20 to about 1e-34 A/um; all port components at -0.1 V agree within 0.057%. Separate EP80 controls pass D1/M3 zero/-0.1 V and both D1/M3 and D2/M3 zero/0.02 V with strict RHS. B/C therefore use EP80 (long double), chosen before G results; EP128 cost controls remain separate. Digits=8 and every frozen gate remain unchanged. Each exact Goal starts with a full segment and may shrink adaptively; required spatial states remain saved. All A inputs, all Vela inputs and all mesh/material inputs remain byte-identical. Original double-precision failures and the frozen exported-state roundoff diagnostic are retained separately. Only complete repeated sweeps qualify the final results.',
            '', 'The first EP128 full sweep additionally exposed premature native stopping at 0.02 V: update error was still 811 when RHS=1.03e-6 met default RHSMin=1e-5. Both EP128/M2 and EP80/M3 strict-RHS probes reduce all 0.02 V component errors below 0.057%. B/C therefore explicitly use RHSMin=1e-20, retaining Digits=8 and all physical parameters. Complete native sweeps are rerun; no reference points are patched from these probes. The prior complete and interrupted attempts are archived under initial_EP128_loose_rhs.',
            '', '## Scope and conventions','',
            'The historical source sweeps to 10 V; the IV template defaults to 20 V. The checked-in legacy fixture IV contract covers 0.2--0.3 V; the separate high-voltage forward guard has anchors at 1, 2, 5, 10, 15, 20 V. This campaign solves 0--0.8 V and 0---1 V and does not requalify those high-voltage anchors for Sentaurus 2022.',
            '', 'The shipped material laws are evaluated at 300 K. SRH uses fixed 1e-5/3e-6 s lifetimes with doping dependence off. Constant mobility is 1417/470.5 cm2/V/s. P3 removes OldSlotboom and its associated -0.01595 eV band-gap reference correction. No parameters are fit to current.',
            '', 'Vela terminal electron/hole columns are q times particle inflow. Conventional current is electron minus hole; the hole column is negated for component comparison. A +x section equals minus the Cathode terminal current. VTK densities use cm^-3; restart CSV density columns serialize SI m^-3 and are read through the matching restart API.',
            '', 'M0/M1 add four left-boundary vertices and local endpoint refinement; M1 halves mesh length limits. The original unsegmented P0 mesh is retained as a separate control. P1/P2/P3 copy the corresponding P0 binary mesh. D1/D2 export fresh node doping. Uncontacted left-boundary segments have no electrode and use insulating natural conditions.',
            '', 'Nodal current vectors are reconstruction diagnostics, not conservative fluxes. The numerical acceptance uses production SG edge cuts at x=0.25, 0.75, 1, 1.25, 1.75 um, alongside terminal KCL and global continuity closure. SRH integrals use common triangular area lumping; endpoint and junction region statistics accompany nodal maxima.',
            '', 'nodal_section_diagnostics.json additionally integrates both tools\' piecewise-linear reconstructed electron, hole and total current fields on the same five sections. These are not native discrete face fluxes and are not declared conserved: at 0.8 V the stage-A junction-line deviations reach 2.67% in the native vector and 9.08% in the Vela vector, while the same Vela states\' exact SG fluxes pass. At P2 reverse leakage, the native nodal integral can differ by 75% (about 5e-19 A/um), and the Vela nodal integral by 13.5%, despite agreeing terminal currents. The raw absolute residuals remain visible; reconstruction errors are not counted as port or SG conservation passes.',
            '', 'The raw solver closure ratio is qualified only when its integrated source exceeds the configured source floor. For P2, zero source can produce a roundoff/roundoff ratio of one. Unqualified points instead receive an explicit electron/hole net-port-current check using the same frozen absolute floor and 1% through-current tolerance. Raw ratios, qualification flags and physical port-net currents are retained in conservation.csv; no threshold is relaxed.',
            '', 'Missing data are not passes. B/C simulations are gated on stage A. Binary TDR/PLT, VTK, states and logs are kept in the ignored build directory.',
            '', 'Legacy current_spreading maps and the left panel of current_spreading_metrics show the original nodal quasi-Fermi-gradient reconstruction, including its large G1/G2 discrepancy. They are retained as diagnostics; use the separately exported SG recovery maps and table for the transport-consistent spreading comparison.',
            '', 'The SG y partitions sum fluxes between node sets on opposite sides of a coordinate threshold. Their staircase dual faces can include longitudinal-current contributions, including in P0; they are not flat-plane integrals of Jy. The recovered SG/native field comparison supports the physical spreading conclusion.',
            '','## Figures','']
    lines += [f'![{p.stem}]({p.name})' for p in sorted(out.glob('*.png'))]
    (out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return out


def plot_legacy_contact_current(r,d,mesh,ids,index,x,y,triangulation,out):
    p=next(p for p in r['spatial']['points'] if p['bias_V']==0.8)
    _,scalars,vectors=read_vtk_point_data(d/p['vtk'])
    ref=node_field(d/'fields/forward_p0p80/fields/TotalCurrentDensity_region0.csv',ids,2)
    paired=[('Sentaurus nodal',ref),('Vela legacy qF-gradient',vectors['SentaurusTotalCurrentDensityVector'])]
    bounds=[(min(0.0,min(v[col] for _,field in paired for v in field)),max(v[col] for _,field in paired for v in field)) for col in [0,1]]
    jy_max=max(abs(v) for v in bounds[1]);bounds[1]=(-jy_max,jy_max)
    anode=next(c for c in mesh['contacts'] if c['name']=='Anode')
    contact_y=[y[index[i]] for i in anode['node_ids']]
    fig,axes=plt.subplots(2,2,figsize=(11,5),constrained_layout=True)
    fig.suptitle(f'{r["id"]}/{r["mesh"]}: legacy reconstruction diagnostic; Anode length {max(contact_y)-min(contact_y):g} um\nCommon column scales; see SG recovery maps for the transport-consistent comparison')
    for row,(tool,field) in enumerate(paired):
        for col in [0,1]:
            artist=axes[row,col].tripcolor(triangulation,[v[col] for v in field],shading='gouraud',cmap='coolwarm',vmin=bounds[col][0],vmax=bounds[col][1])
            axes[row,col].set(title=f'{tool} J{"xy"[col]} [A/cm2] at 0.8 V',xlabel='x [um]',ylabel='y [um]',aspect='equal')
            axes[row,col].plot([min(x),min(x)],[min(contact_y),max(contact_y)],color='black',lw=3)
            axes[row,col].plot([max(x),max(x)],[min(y),max(y)],color='black',lw=3)
            fig.colorbar(artist,ax=axes[row,col],shrink=.7)
    fig.savefig(out/f'{r["id"]}_current_spreading.png',dpi=170); plt.close(fig)


def plot_spreading_metrics(preferred,out,root):
    spreading=[]
    fig,axes=plt.subplots(1,2,figsize=(10,3.5),constrained_layout=True)
    for identifier in ['P0','G1','G2']:
        r=preferred.get(identifier)
        if not r or 'points' not in r['spatial']: continue
        point=next(p for p in r['spatial']['points'] if p['bias_V']==.8)
        vector=next(v for v in point['current_vectors'] if v['field']=='TotalCurrentDensity')
        for tool in ['reference','candidate']:
            value=vector.get(tool+'_Jy_over_Jx_l2')
            if value is None: continue
            spreading.append({'id':identifier,'mesh':r['mesh'],'bias_V':.8,'tool':tool,'Jy_over_Jx_l2':value})
            axes[0].plot(identifier,value,'o' if tool=='reference' else 'x',color='C0' if tool=='reference' else 'C1')
        cuts=point.get('transverse_sg_sections',[])
        axes[1].plot([p['y_um'] for p in cuts],[p['total_A_per_um'] for p in cuts],'.-',label=identifier)
    axes[0].set(title='Legacy reconstruction diagnostic at 0.8 V\nCircle: Sentaurus; cross: Vela qF-gradient',ylabel='Area-weighted ||Jy|| / ||Jx||')
    axes[1].set(title='Vela SG dual-face y partitions at 0.8 V',xlabel='Node partition threshold y [um]',ylabel='Current toward upper partition [A/um]')
    if axes[1].get_legend_handles_labels()[0]: axes[1].legend()
    for ax in axes: ax.grid(alpha=.2)
    write_json(root/'current_spreading.json',spreading)
    fig.savefig(out/'current_spreading_metrics.png',dpi=160);plt.close(fig)


def plot_terminal_components(preferred,out):
    # Expose every measured port component, with the same tool convention as
    # the exact-point tables. Zero-bias absolute gates remain in the tables.
    for branch in ['forward','reverse']:
        fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
        for color,(identifier,r) in enumerate(preferred.items()):
            for row,contact in enumerate(['Anode','Cathode']):
                for col,component in enumerate(['total','electron','hole']):
                    points=sorted([v for v in r['terminal']['metrics'] if v['branch']==branch and v['contact']==contact and v['component']==component and v['bias_V']!=0],key=lambda v:v['bias_V'])
                    ax=axes[row,col]
                    for tool,style in [('reference','-'),('candidate','--')]:
                        ax.semilogy([v['bias_V'] for v in points],[max(abs(v[tool+'_A_per_um']),1e-300) for v in points],style,color=f'C{color}',label=identifier if tool=='reference' else None)
                    ax.set(title=f'{contact}: {component}',xlabel='Anode voltage [V]',ylabel='Absolute current [A/um]')
                    ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=7,ncol=2)
        fig.suptitle(f'{branch}: solid Sentaurus, dashed independent Vela; signed values and zero-bias gates in terminal_points.csv')
        fig.savefig(out/f'terminal_components_{branch}.png',dpi=150);plt.close(fig)


def export_evidence(root, out):
    """Publish compact numerical evidence locally, never binary solver states."""
    data=read_json(root/'results.json')
    dest=FIXTURE/'variants/results';dest.mkdir(exist_ok=True)
    def table(name, records):
        if not records: return
        fields=list(dict.fromkeys(k for r in records for k in r))
        with (dest/name).open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(records)
    summary={k:v for k,v in data.items() if k not in ['results','mesh','refined_mesh','historical_refinement_mesh','contact_edge_mesh','model_effects']}
    summary['cases']=[]; ports=[]; numerical=[]; fields=[]; mesh=[]; effects=[]
    for r in data['results']:
        d=root/r['id']/r['mesh']; key={'id':r['id'],'mesh':r['mesh']}
        item={**key,'status':r['status'],'required_for_final_qualification':r.get('required_for_final_qualification',False),
              'input':r['input'],'evidence_sha256':r.get('evidence_sha256',{}),
              'terminal_status':r['terminal']['status'],'spatial_status':r['spatial']['status'],
              'conservation_status':r['spatial'].get('transport_source_and_section_status','not_run')}
        for name in ['sentaurus_run.json','mesh_run.json','mesh_input_comparison.json','vela_forward_run.json','vela_reverse_run.json']:
            if (d/name).exists(): item[name]=read_json(d/name)
        summary['cases'].append(item)
        ports.extend({**key,**p} for p in r['terminal'].get('metrics',[]))
        numerical.extend({**key,**p} for p in r['terminal'].get('numerical',[]))
        fields.append({**key,**r['spatial']})
    def mesh_evidence(gates,pair):
        compact={}
        for identifier,g in gates.items():
            compact[identifier]={k:v for k,v in g.items() if k!='terminal_metrics'}
            mesh.extend({'pair':pair,'id':identifier,**p} for p in g.get('terminal_metrics',[]))
        return compact
    summary['mesh']=mesh_evidence(data['mesh'],'M0/M1')
    summary['refined_mesh']=mesh_evidence(data.get('refined_mesh',{}),'/'.join(data['refinement']['additional_pair']))
    summary['contact_edge_mesh']=mesh_evidence(data.get('contact_edge_mesh',{}),'/'.join(data['refinement']['contact_edge_control']['pair']))
    summary['historical_refinement_mesh']={pair:mesh_evidence(gates,pair) for pair,gates in data.get('historical_refinement_mesh',{}).items()}
    summary['model_effects']={}
    for identifier,grids in data['model_effects'].items():
        summary['model_effects'][identifier]={m:g['status'] for m,g in grids.items()}
        for m,g in grids.items(): effects.extend({'id':identifier,'mesh':m,**p} for p in g.get('metrics',[]))
    write_json(dest/'summary.json',summary);write_json(dest/'spatial_metrics.json',fields)
    table('terminal_points.csv',ports);table('conservation.csv',numerical);table('mesh_metrics.csv',mesh);table('model_effects.csv',effects)
    for name in ['parameter_audit.json','environment.json','verification.json','checkpoint.json','junction_resolution.json','current_spreading.json','continuity_semantics.json','boundary_coverage.json','nodal_section_diagnostics.json','contact_state_audit.json']:
        if (root/name).exists(): shutil.copyfile(root/name,dest/name)
    for name in ['native_precision_repair.json','native_precision_input_preservation.json','initial_native_double_failures.json','native_stopping_repair.json','EP128_cost_controls.json']:
        if (root/name).exists(): shutil.copyfile(root/name,dest/name)
    for source,name in [('diagnostics/native_precision/D1_M2_EP128/comparison.json','native_precision_probe.json'),
                        ('diagnostics/native_precision/D1_M2_EP128/run.json','native_precision_probe_run.json'),
                        ('diagnostics/native_precision/D1_M3_EP80/comparison.json','native_precision_80_probe.json'),
                        ('diagnostics/native_precision/D1_M3_EP80/run.json','native_precision_80_probe_run.json'),
                        ('diagnostics/native_precision/D1_M2_EP128_forward_rhs1e-20/comparison.json','native_strict_rhs_128_probe.json'),
                        ('diagnostics/native_precision/D1_M3_EP80_forward_rhs1e-20/comparison.json','native_strict_rhs_80_probe.json'),
                        ('diagnostics/native_precision/D2_M3_EP80_forward_rhs1e-20/comparison.json','native_strict_rhs_80_D2_probe.json'),
                        ('diagnostics/native_precision/EP128_EP80_M2_overlap.json','native_precision_overlap.json'),
                        ('diagnostics/equilibrium_roundoff/D1_M2/result.json','frozen_equilibrium_roundoff.json')]:
        if (root/source).exists(): shutil.copyfile(root/source,dest/name)
    figures=dest/'figures';figures.mkdir(exist_ok=True)
    selected=['terminal_curves.png','terminal_components_forward.png','terminal_components_reverse.png','model_effects.png','mesh_sensitivity.png','junction_resolution.png','current_spreading_metrics.png','nodal_section_diagnostics.png']
    for identifier in [e['id'] for e in read_json(SPEC)['experiments']]:
        candidate=next((out/f'{identifier}_{m}_states.png' for m in reversed(list(read_json(SPEC)['meshes'])) if (out/f'{identifier}_{m}_states.png').exists()),None)
        if candidate: selected.append(candidate.name)
    selected.extend(p.name for p in out.glob('G*_current_spreading.png'))
    for name in selected:
        if (out/name).exists(): shutil.copyfile(out/name,figures/name)
    recovery_rows=[]
    for path in (root/'diagnostics/current_reconstruction').glob('*/result.json'):
        recovery=read_json(path)
        for state in recovery['states']:
            for item in state['recoveries']:
                if 'Total' in item['field']:
                    recovery_rows.append(f"| {recovery['id']}/{recovery['mesh']} | {state['state']} | {item['field']} | {100*item['regions']['all']['relative_l2']:.4f}% | {item['regions']['all']['Jy_over_Jx_l2']:.6f} |")
        name=f'current_reconstruction_{recovery["id"]}_{recovery["mesh"]}.json'
        shutil.copyfile(path,dest/name)
        if recovery.get('figure') and (path.parent/recovery['figure']).exists():
            name=recovery['figure'];shutil.copyfile(path.parent/name,figures/name);selected.append(name)
    text=(out/'report.md').read_text(encoding='utf-8')
    text=text.split('## Figures')[0]
    text=text.replace('`../contract.json` and `../results.json`','`../acceptance.json`, `summary.json`, `spatial_metrics.json` and the adjacent exact-point CSV tables')
    if recovery_rows:
        text+='## Current-vector recovery diagnosis\n\nThe legacy nodal quasi-Fermi-gradient current is a reconstruction diagnostic. Its large local-contact discrepancy is reproduced on the unchanged exported native state. Existing SG recoveries agree much more closely without solving again, fitting currents, changing defaults, or changing gates. The original fields remain recorded; the additional SG maps use shared symlog scales. These recoveries are not native internal face fluxes.\n\n| Case | Frozen state | Recovery | Weighted L2 error | Jy/Jx norm |\n|---|---|---|---|---|\n'+'\n'.join(recovery_rows)+'\n\n'
    text+='## Figures\n\n'+'\n\n'.join(f'![{Path(name).stem}](figures/{name})' for name in selected if (figures/name).exists())+'\n'
    (dest/'technical_report.md').write_text(text,encoding='utf-8')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--export-evidence',action='store_true')
    args=parser.parse_args();root=args.root.resolve();out=report(root)
    if args.export_evidence: export_evidence(root,out)
