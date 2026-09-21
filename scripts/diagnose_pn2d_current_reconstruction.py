#!/usr/bin/env python3
"""Compare existing current recoveries on frozen Vela and exported native states."""
import argparse
import csv
import math
from pathlib import Path
from pn2d_variants import DEFAULT_ROOT, REPO, read_json, write_json, run_cmd, sha, tag
from compare_pn2d_variants import node_field, read_vtk_point_data, spatial_regions


def metrics(reference, candidate, weights, mask):
    selected=[i for i,keep in enumerate(mask) if keep]
    area=sum(weights[i] for i in selected)
    norm=sum(weights[i]*sum(x*x for x in reference[i]) for i in selected)
    error=sum(weights[i]*sum((candidate[i][j]-reference[i][j])**2 for j in range(2)) for i in selected)
    transverse=sum(weights[i]*candidate[i][1]**2 for i in selected)
    longitudinal=sum(weights[i]*candidate[i][0]**2 for i in selected)
    return {'nodes':len(selected),'area_um2':area,'relative_l2':math.sqrt(error/max(norm,1e-300)),
            'rms_error_A_per_cm2':math.sqrt(error/max(area,1e-300)),
            'squared_error_integral':error,'Jy_over_Jx_l2':math.sqrt(transverse/max(longitudinal,1e-300))}


def plot_recovery(mesh, fields, out, identifier, grid, bias):
    ids=[n['id'] for n in mesh['nodes']]
    positions={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    index={node:i for i,node in enumerate(ids)}
    triangulation=mtri.Triangulation([n['x'] for n in mesh['nodes']],[n['y'] for n in mesh['nodes']],[[index[i] for i in t['node_ids']] for t in mesh['triangles']])
    reference=node_field(fields/'TotalCurrentDensity_region0.csv',ids,2)
    _,_,recovered=read_vtk_point_data(out/'vela.vtk');candidate=recovered['CellFirstSgTotalCurrentDensityVector']
    fig,axes=plt.subplots(2,2,figsize=(11,5.5),constrained_layout=True)
    for col in range(2):
        values=[v[col] for field in [reference,candidate] for v in field];low,high=min(values),max(values)
        norm=matplotlib.colors.SymLogNorm(linthresh=max(abs(low),abs(high))*.01 or 1.,vmin=low,vmax=high)
        for row,(label,field) in enumerate([('Sentaurus nodal output',reference),('Vela cell-first SG recovery',candidate)]):
            ax=axes[row,col];art=ax.tripcolor(triangulation,[v[col] for v in field],shading='gouraud',cmap='coolwarm',norm=norm)
            ax.set(title=f'{label}: J{"xy"[col]} [A/cm2]',xlabel='x [um]',ylabel='y [um]',aspect='equal')
            for contact in mesh['contacts']:
                points=[positions[i] for i in contact['node_ids']]
                ax.plot([points[0][0]]*2,[min(v[1] for v in points),max(v[1] for v in points)],color='black',lw=3)
            fig.colorbar(art,ax=ax,shrink=.8)
    fig.suptitle(f'{identifier}/{grid} at {bias:g} V: existing cell-first SG recovery; shared symlog scales')
    figure=out/f'{identifier}_SG_current_spreading.png';fig.savefig(figure,dpi=160);plt.close(fig)
    return figure


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--id',choices=['P0','G1','G2'],required=True)
    p.add_argument('--mesh',default='E1')
    p.add_argument('--bias',type=float,default=.8)
    a=p.parse_args();root=a.root.resolve();d=root/a.id/a.mesh
    out=root/'diagnostics/current_reconstruction'/f'{a.id}_{a.mesh}_{tag(a.bias)}';out.mkdir(parents=True,exist_ok=True)
    mesh=read_json(d/'inputs/mesh.json');ids=[n['id'] for n in mesh['nodes']];positions={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    area={i:0. for i in ids}
    for t in mesh['triangles']:
        x,y,z=[positions[i] for i in t['node_ids']]
        value=abs((y[0]-x[0])*(z[1]-x[1])-(z[0]-x[0])*(y[1]-x[1]))/6
        for i in t['node_ids']:area[i]+=value
    weights=[area[i] for i in ids]
    regions=spatial_regions(mesh,read_json(root/'contract.json'))
    regions.update(all=[True]*len(ids),outside_contact_edges=[not x for x in regions['contact_edges']],
                   x_ge_0p1_um=[n['x']>=.1 for n in mesh['nodes']])
    branch='reverse' if a.bias<0 else 'forward';fields=d/'fields'/f'{branch}_{tag(a.bias)}'/'fields'
    native={name:node_field(fields/f'{field}_region0.csv',ids) for name,field in
            [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('electrons_m3','eDensity'),('holes_m3','hDensity')]}
    frozen=out/'native_exported_state.csv'
    with frozen.open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['node_id',*native])
        for j,i in enumerate(ids):w.writerow([i,*[v[j]*(1e6 if k.endswith('_m3') else 1) for k,v in native.items()]])
    artifacts=read_json(d/'spatial_artifacts.json');own=d/artifacts[f'{a.bias:g}']['state_file']
    result={'scope':'Diagnostic recovery only: unchanged frozen states; native replay means exported fields, not native internal fluxes. No transport solve, current fitting, default or threshold changes.',
            'id':a.id,'mesh':a.mesh,'bias_V':a.bias,'runner_sha256':sha(REPO/'build-release/vela_example_runner.exe'),
            'states':[],'reference_field_sha256':{q.name:sha(q) for q in fields.glob('*.csv')}}
    for label,state in [('vela',own),('native_exported',frozen)]:
        cfg=read_json(d/f'vela_{branch}.json');cfg.pop('sweep');cfg.pop('output_csv',None)
        for key in ['mesh_file','node_doping_file','materials_file']:cfg[key]=str(d/cfg[key])
        cfg.update(simulation_type='write_dd_state_vtk',state_file=str(state),output_vtk=str(out/f'{label}.vtk'),
                   output_diagnostics={'cell_first_sg_current_recovery':True})
        for c in cfg['contacts']:c['bias']=a.bias if c['name']=='Anode' else 0.
        config=out/f'{label}.json';write_json(config,cfg)
        rc=run_cmd([REPO/'build-release/vela_example_runner.exe','--config',config],REPO,out/f'{label}.log')
        if rc:raise RuntimeError(f'recovery failed: {label}')
        count,_,vectors=read_vtk_point_data(Path(cfg['output_vtk']));assert count==len(ids)
        record={'state':label,'state_sha256':sha(state),'config_sha256':sha(config),'vtk_sha256':sha(Path(cfg['output_vtk'])),'recoveries':[]}
        for carrier,native_name in [('Electron','eCurrentDensity'),('Hole','hCurrentDensity'),('Total','TotalCurrentDensity')]:
            ref=node_field(fields/f'{native_name}_region0.csv',ids,2)
            for prefix in ['Sentaurus','DualFaceSg','CellFirstSg']:
                key=prefix+carrier+'CurrentDensityVector';v=vectors[key]
                record['recoveries'].append({'field':key,'regions':{name:metrics(ref,v,weights,mask) for name,mask in regions.items()}})
        result['states'].append(record)
        print(a.id,label,[(x['field'],x['regions']['all']['relative_l2']) for x in record['recoveries'] if 'Total' in x['field']],flush=True)
    result['figure']=plot_recovery(mesh,fields,out,a.id,a.mesh,a.bias).name
    write_json(out/'result.json',result)


if __name__=='__main__':main()
