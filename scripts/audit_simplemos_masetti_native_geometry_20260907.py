"""Use native box primitives; keep unresolved mobility averaging tests diagnostic."""
from pathlib import Path
import math
import numpy as np
from scipy.sparse import coo_matrix
import prepare_simplemos_masetti_local_20260907 as p
import run_simplemos_m34_sentaurus_interface_box_probe as box
import audit_simplemos_subset_mobility_export_20260907 as mobility

a=p.a;d=p.d;OUT=p.OUT;LOCAL=p.LOCAL;sc=d.matrix.spatial.m73.scalar


def geometry(device):
    src=LOCAL/'native_exports'/device
    xy={int(x['id']):(float(x['x_um']),float(x['y_um'])) for x in a.rows(src/'nodes.csv')}
    elements={int(x['id']):dict(nodes=[int(x['node'+str(k)]) for k in range(3)],material=x['material']) for x in a.rows(src/'elements.csv')}
    debug=(LOCAL/'native_raw/bundle'/device/'MeasureCoefficients.debug').read_text()
    coeff=box.parse_debug_block(debug,'Coefficients');measure=box.parse_debug_block(debug,'Measure')
    assert set(elements)=={k for k,v in coeff.items() if v['type']==2 and v['des']>=0}
    # Use the same 40 Si/SiO2 interface cells as the independently checked M34
    # numbering witness. Other regions contain actual native box modifications:
    # requiring every acute cell to equal an unmodified cotangent was false.
    adjacency={}
    for i,e in elements.items():
        for k in range(3):adjacency.setdefault(tuple(sorted((e['nodes'][(k+1)%3],e['nodes'][(k+2)%3]))),[]).append(i)
    anchors={i for adjacent in adjacency.values() if {elements[j]['material'] for j in adjacent}=={'Si','SiO2'} for i in adjacent}
    assert len(anchors)==40
    cp,ce=box.infer_input_to_debug_local_permutation(anchors,elements,xy,coeff)
    vp,ve=box.infer_input_to_debug_measure_permutation(anchors,elements,xy,measure)
    assert ce<1e-8 and ve<1e-12
    modified=[]
    for i,e in elements.items():
        pts=[xy[n] for n in e['nodes']];expected=box.signed_average_box_measures(pts)
        cerror=max(abs(coeff[i]['values'][cp[k]]-box.raw_coefficient(pts,k)) for k in range(3))
        verror=max(abs(measure[i]['values'][vp[k]]-expected[k]) for k in range(3))
        modified.append(dict(device=device,cell=i,material=e['material'],coefficient_max_absolute_difference=cerror,
            measure_max_absolute_difference_um2=verror,measure_sum_minus_triangle_area_um2=sum(measure[i]['values'])-sum(expected)))
    geo=d.matrix.spatial.m73.Geometry(device);coordinate=max(math.dist(xy[i],geo.coords[i]) for i in xy);assert coordinate<=1e-12
    mat=a.read(Path(a.read(Path(a.read(p.prev.OUT/'vela_contract.json')['jobs'][0]['config']))['materials_file']))
    eps={m['name']:m['eps_r']*d.matrix.spatial.m73.EPS0 for m in mat['materials']}
    volumes={m:np.zeros(geo.count) for m in eps};rc=[];cc=[];data=[];edgeparts={}
    for i,e in elements.items():
        nodes=e['nodes'];material=e['material']
        for k,n in enumerate(nodes):volumes[material][n]+=measure[i]['values'][vp[k]]*1e-12
        for k in range(3):
            u,v=sorted((nodes[(k+1)%3],nodes[(k+2)%3]));g=coeff[i]['values'][cp[k]]
            edgeparts.setdefault((u,v),[]).append(dict(cell=i,material=material,coefficient=g))
            value=eps[material]*g;rc.extend((u,u,v,v));cc.extend((u,v,u,v));data.extend((value,-value,-value,value))
    K=coo_matrix((data,(rc,cc)),shape=(geo.count,geo.count)).tocsr()
    info=dict(device=device,coordinate_max_um=coordinate,coefficient_permutation=list(cp),coefficient_mapping_max_absolute=ce,
        measure_permutation=list(vp),measure_mapping_max_um2=ve,cells=len(elements),
        signed_Si_volume_max_difference_um2=float(max(abs(volumes['Si']-geo.volumes['signed_si']))*1e12),
        legacy_K_relative_frobenius=float(np.linalg.norm((geo.matrices['legacy']-K).data)/np.linalg.norm(K.data)),
        region_local_K_relative_frobenius=float(np.linalg.norm((geo.matrices['region_local']-K).data)/np.linalg.norm(K.data)),
        cells_differing_from_raw_coefficients_at_1e_minus_10=sum(r['coefficient_max_absolute_difference']>1e-10 for r in modified),
        cells_differing_from_signed_measure_at_1e_minus_14=sum(r['measure_max_absolute_difference_um2']>1e-14 for r in modified))
    return geo,xy,elements,volumes,K,edgeparts,info,modified


def main():
    a.verify(OUT/'native_freeze.json');a.verify(OUT/'vela_freeze.json')
    allinfo=[];murows=[];edgeout=[];poisson=[];local=[];modifications=[];inputs=[Path(__file__).resolve(),mobility.__file__,box.__file__]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge'];Q=d.fixed.Q
    for device in ('n19','n23'):
        geo,xy,elements,volumes,K,edgeparts,info,modified=geometry(device);allinfo.append(info);modifications.extend(modified)
        src=LOCAL/'native_exports'/device;old=p.prev.LOCAL/f'native_exports/masetti/m65_{device}_vd_0p050000_endpoint/vg_020'
        nd=sc(src/'fields/DonorConcentration_region0.csv');na=sc(src/'fields/AcceptorConcentration_region0.csv')
        si=[i for i,e in elements.items() if e['material']=='Si'];assert si==list(range(len(si)))
        for carrier,pars in mobility.PARAMETERS.items():
            mu={i:mobility.formula(nd[i]+na[i],pars) for i in nd};native=np.array([float(x['component0']) for x in a.rows(src/f'fields/{carrier}Mobility_region0_cells.csv')])
            nodeplot=sc(old/f'fields/{carrier}Mobility_region0.csv')
            vertex=np.array([[mu[n] for n in elements[i]['nodes']] for i in si]);Ns=np.array([[nd[n]+na[n] for n in elements[i]['nodes']] for i in si])
            candidates={'arithmetic_local_vertex_mu':vertex.mean(axis=1),'harmonic_local_vertex_mu':3/(1/vertex).sum(axis=1),
                'geometric_local_vertex_mu':np.exp(np.log(vertex).mean(axis=1)),
                'formula_at_cell_mean_doping':np.array([mobility.formula(x,pars) for x in Ns.mean(axis=1)]),
                'arithmetic_native_node_plot':np.array([sum(nodeplot[n] for n in elements[i]['nodes'])/3 for i in si])}
            for name,value in candidates.items():
                error=abs(value/native-1)
                murows.append(dict(device=device,carrier=carrier,candidate=name,cells=len(si),median_relative=float(np.median(error)),p95_relative=float(np.percentile(error,95)),max_relative=float(max(error)),
                    worst_cell=int(np.argmax(error)),pointwise_identity_at_1e_minus_8=bool(max(error)<=1e-8)))
        masks,_,_=d.matrix.spatial.old.m78.supports(device,geo,.05)
        cfg=next(c for c in a.read(OUT/'vela_contract.json')['cases'] if c['device']==device)
        velaedges={tuple(sorted((int(e['node0']),int(e['node1'])))):e for e in a.rows(LOCAL/'vela'/cfg['key']/'strict/edges.csv')}
        for edge,parts in edgeparts.items():
            mats=set(x['material'] for x in parts)
            if 'Si' not in mats or len(mats)==1:continue
            gn=sum(x['coefficient'] for x in parts if x['material']=='Si')
            gv=float(velaedges[edge]['couple_m'])/float(velaedges[edge]['length_m'])
            kn=-K[edge[0],edge[1]];kv=-geo.matrices['legacy'][edge[0],edge[1]]
            edgeout.append(dict(device=device,node0=edge[0],node1=edge[1],materials=';'.join(sorted(mats)),
                native_Si_transport_geometry=gn,vela_all_cell_transport_geometry=gv,
                native_poisson_coefficient_F_per_m=float(kn),vela_poisson_coefficient_F_per_m=float(kv),
                native_Si_transport_is_zero=gn==0,geometry_only=True))
        for c in a.read(OUT/'vela_contract.json')['cases']:
            if c['device']!=device:continue
            root=LOCAL/'vela'/c['key'];rawsrc=p.prev.LOCAL/'native_exports/masetti'/c['case']/f"vg_{c['index']:03d}"
            psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(rawsrc,geo)
            assert spread<=1e-12
            rawnd=sc(rawsrc/'fields/DonorConcentration_region0.csv');rawna=sc(rawsrc/'fields/AcceptorConcentration_region0.csv')
            doping=np.array([(rawnd.get(i,0)-rawna.get(i,0))*1e6 for i in range(geo.count)])
            rawcharge=n-h-doping
            mapped=d.array(d.ordered(Path(c['mapped']),geo.count),('psi','electrons_m3','holes_m3'))
            probe=d.array(d.ordered(root/'mapped/residual.csv',geo.count),('psi_residual',))[0]/factor
            nativeR=K@psi+Q*rawcharge*volumes['Si']
            dielectric=(geo.matrices['legacy']-K)@psi
            volume=Q*rawcharge*(geo.volumes['all_cell']-volumes['Si'])
            convention=Q*((mapped[1]-n)-(mapped[2]-h))*geo.volumes['all_cell']
            reconstructed=nativeR+dielectric+volume+convention
            remainder=probe-reconstructed
            den=np.linalg.norm(probe[geo.free]);charge_norm=np.linalg.norm((Q*rawcharge*volumes['Si'])[geo.free])
            poisson.append(dict(key=c['key'],native_reconstructed_residual_over_charge=float(np.linalg.norm(nativeR[geo.free])/max(charge_norm,1e-300)),
                native_reconstructed_residual_over_vela_mapped_residual=float(np.linalg.norm(nativeR[geo.free])/max(den,1e-300)),
                vela_reconstruction_remainder_relative=float(np.linalg.norm(remainder[geo.free])/max(den,1e-300)),
                interpretation='Reconstructed native Poisson expression, not exported native row residual. Native charge quadrature still requires equivalence check.'))
            for name,source in [('native_reconstructed_residual',nativeR),('dielectric_geometry',dielectric),('charge_volume',volume),('thermal_density_convention',convention),('vela_reconstruction_remainder',remainder)]:
                response=geo.solve('legacy',source)
                for node in (320,324,338):local.append(dict(key=c['key'],component=name,node=node,conditional_poisson_response_V=float(response[node])))
    for name,rows in [('native_geometry',allinfo),('native_cell_geometry',modifications),('element_mobility_candidates',murows),('interface_geometry',edgeout),('native_poisson_replay',poisson),('native_geometry_poisson_components',local)]:a.write_csv(OUT/(name+'.csv'),rows)
    paths=[Path(x) for x in inputs]+[LOCAL/'results.tgz',OUT/'native_checks.csv']
    a.write(OUT/'native_geometry_audit.json',dict(status='completed_diagnostic_audit',input_hashes={a.rel(x):a.sha(x) for x in paths},
        mobility='No averaging candidate is accepted unless its all-cell maximum relative error is <=1e-8. Failed candidates are retained; cell plots are not silently equated to internal edge mobility.',
        poisson='Native coefficients/measures are actual exported geometry. The reconstructed charge equation and derived responses are kept distinct from an exported native residual or a calibrated self-consistent correction.',
        fitted_physical_parameters=False,production_changes=False))
    print('Native geometry and mobility/Poisson replay recorded',flush=True)


if __name__=='__main__':main()
