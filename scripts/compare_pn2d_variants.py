#!/usr/bin/env python3
"""Fail-closed exact-point terminal and state audit for PN2D variants."""
from __future__ import annotations
import argparse
import csv
import math
import hashlib
import json
import re
from pathlib import Path
from pn2d_variants import (REPO, SPEC, DEFAULT_ROOT, read_json, write_json, sha, biases, tag,
                          run_cmd, parse_quoted_list, parse_values_block)
from compare_genius_bjt_transport_fields import read_vtk_point_data
from audit_genius_bjt_conservative_sections import section_current


def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as f:
        data=list(csv.DictReader(f))
    if any(None in r or any(v is None for v in r.values()) for r in data):
        raise ValueError(f'incomplete CSV row: {path}')
    return data


def exact_row(data, bias, key='bias_V', contact=None):
    found=[r for r in data if abs(float(r[key])-bias)<=1e-9 and (contact is None or r['contact']==contact)]
    if not found:
        raise ValueError(f'missing solved bias {bias} contact {contact}')
    # SDevice can repeat a Goal endpoint as the next segment's starting point.
    # Keep the final native record, never interpolate between records.
    if any(any(not math.isfinite(float(v)) for k,v in r.items() if k==key) for r in found):
        raise ValueError('nonfinite bias')
    return found[-1]


def plt_rows(path):
    text=path.read_text(errors='strict')
    fields=parse_quoted_list(text,'datasets')
    if not fields: raise ValueError('no PLT datasets')
    return [dict(zip(fields,v,strict=True)) for v in parse_values_block(text,len(fields))]


def compare_current(ref, candidate, contract):
    if not math.isfinite(ref) or not math.isfinite(candidate):
        raise ValueError('nonfinite current')
    floor=contract['absolute_floor_A_per_um']
    absolute=abs(candidate-ref)
    relative=absolute/max(abs(ref),floor)
    active=abs(ref)>floor
    same_sign=ref*candidate>0
    log_error=abs(math.log10(abs(candidate/ref))) if active and same_sign else None
    passed=absolute<=floor or (active and same_sign and relative<=contract['relative_max'] and log_error<=contract['log_max_dex'])
    return {'reference_A_per_um':ref,'candidate_A_per_um':candidate,'absolute_error_A_per_um':absolute,
            'relative_error':relative,'log_error_dex':log_error,'pass':passed}


def continuity_satisfied(state, contract, ports=None, source_floor=1e-14):
    flag=str(state.get('global_continuity_closure_satisfied','')).lower() in ('1','true')
    if not flag: return False
    for carrier in ['electron','hole']:
        ratio=float(state[f'global_{carrier}_continuity_closure_ratio'])
        source=float(state[f'global_{carrier}_integrated_source'])
        if not math.isfinite(ratio) or ratio<0 or not math.isfinite(source): return False
        if abs(source)>=source_floor:
            if ratio>contract['conservation']['continuity_relative_max']: return False
        else:
            # The solver ratio divides a NET flux mismatch by net flux/source.
            # With SRH off, roundoff / roundoff can equal one. Verify actual
            # component conservation against through-current and the frozen
            # absolute floor instead of treating this undefined ratio as error.
            if ports is None: return False
            currents=[float(p[f'current_{carrier}_A_per_um']) for p in ports.values()]
            if not all(math.isfinite(v) for v in currents): return False
            tolerance=max(contract['conservation']['absolute_floor_A_per_um'],
                          max(abs(v) for v in currents)*contract['conservation']['continuity_relative_max'])
            if abs(sum(currents))>tolerance: return False
    return True


def terminals(d, contract):
    spec=read_json(SPEC); output=[]; numerical=[]
    for branch in spec['branches']:
        config=read_json(d/f'vela_{branch}.json')
        reference=plt_rows(d/f'{branch}.plt')
        candidate=rows(d/f'{branch}_terminal_balance.csv')
        curve=rows(d/f'{branch}.csv')
        for bias in biases(spec['branches'][branch]):
            r=exact_row(reference,bias,'Anode OuterVoltage')
            c={name:exact_row(candidate,bias,contact=name) for name in ['Anode','Cathode']}
            state=exact_row(curve,bias)
            for name in c:
                for component,refkey in [('total','TotalCurrent'),('electron','eCurrent'),('hole','hCurrent')]:
                    # Vela stores q times particle inflow: conventional hole current
                    # is minus its hole column; total = electron - hole.
                    factor=-1.0 if component=='hole' else 1.0
                    metric=compare_current(float(r[f'{name} {refkey}']),factor*float(c[name][f'current_{component}_A_per_um']),contract['current'])
                    output.append({'branch':branch,'bias_V':bias,'contact':name,'component':component,**metric})
            vkcl=abs(sum(float(c[n]['current_total_A_per_um']) for n in c))
            skcl=abs(sum(float(r[f'{n} TotalCurrent']) for n in c))
            scale=max(abs(float(r['Anode TotalCurrent'])),abs(float(c['Anode']['current_total_A_per_um'])))
            tolerance=max(contract['conservation']['absolute_floor_A_per_um'],contract['conservation']['kcl_relative_max']*scale)
            converged=str(state.get('converged','')).lower() in ('1','true')
            continuity=continuity_satisfied(state,contract,c,config['solver']['global_continuity_closure']['source_floor'])
            closure_details={'continuity_pass':continuity}
            for carrier in ['electron','hole']:
                factor=-1 if carrier=='hole' else 1
                currents=[factor*float(p[f'current_{carrier}_A_per_um']) for p in c.values()]
                closure_details[carrier+'_source_ratio_qualified']=abs(float(state[f'global_{carrier}_integrated_source']))>=config['solver']['global_continuity_closure']['source_floor']
                closure_details[carrier+'_port_net_A_per_um']=sum(currents)
                closure_details[carrier+'_port_net_over_through_current']=abs(sum(currents))/max(contract['conservation']['absolute_floor_A_per_um'],*(abs(v) for v in currents))
            numerical.append({'branch':branch,'bias_V':bias,'vela_kcl_A_per_um':vkcl,'sentaurus_kcl_A_per_um':skcl,
                              'electron_continuity_ratio':float(state['global_electron_continuity_closure_ratio']),
                              'hole_continuity_ratio':float(state['global_hole_continuity_closure_ratio']),
                              **closure_details,'converged':converged,'pass':converged and continuity and vkcl<=tolerance and skcl<=tolerance})
    logs=[r['log_error_dex'] for r in output if r['log_error_dex'] is not None]
    rmse=math.sqrt(sum(x*x for x in logs)/len(logs)) if logs else None
    return {'status':'pass' if all(r['pass'] for r in output) and rmse is not None and rmse<=contract['current']['log_rmse_max_dex'] else 'fail',
            'log_rmse_dex':rmse,'metrics':output,'numerical':numerical,
            'convergence_kcl_status':'pass' if all(r['pass'] for r in numerical) else 'fail'}


def percentile(a,f):
    a=sorted(a)
    if not a: raise ValueError('empty active population')
    return a[min(len(a)-1,round((len(a)-1)*f))]


def node_field(path, ids, components=1):
    data=rows(path); indexed={int(r['node_id']):r for r in data}
    if len(indexed)!=len(data) or set(indexed)!=set(ids): raise ValueError(f'node identity mismatch: {path}')
    result=[[float(indexed[i][f'component{k}']) for k in range(components)] for i in ids]
    if not all(math.isfinite(v) for r in result for v in r): raise ValueError(f'nonfinite field: {path}')
    return [r[0] for r in result] if components==1 else result


def validate_vtk_coordinates(path, expected, tolerance):
    with path.open(encoding='utf-8') as f:
        for line in f:
            if line.startswith('POINTS '):
                if int(line.split()[1])!=len(expected): raise ValueError('VTK point count mismatch')
                break
        else: raise ValueError('VTK POINTS missing')
        values=[]
        while len(values)<3*len(expected): values.extend(float(v) for v in next(f).split())
    for i,(x,y) in enumerate(expected):
        if max(abs(values[3*i]-x),abs(values[3*i+1]-y),abs(values[3*i+2]))>tolerance:
            raise ValueError('VTK point ordering or coordinates differ from the common mesh')


def spatial_regions(mesh, contract):
    """Use physical contact endpoints, including corners for full-edge contacts."""
    positions={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    anode=next(c for c in mesh['contacts'] if c['name']=='Anode')
    points=[positions[i] for i in anode['node_ids']]
    endpoints=[min(points,key=lambda p:p[1]),max(points,key=lambda p:p[1])]
    radius=contract['spatial']['contact_edge_region_radius_um']
    junction=read_json(SPEC)['geometry_um']['junction_x']
    return {'contact_edges':[min(math.hypot(n['x']-x,n['y']-y) for x,y in endpoints)<=radius for n in mesh['nodes']],
            'junction':[abs(n['x']-junction)<=0.1 for n in mesh['nodes']]}


def spatial(d, contract):
    mesh=read_json(d/'inputs/mesh.json'); ids=[n['id'] for n in mesh['nodes']]
    positions={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    areas={i:0.0 for i in ids}
    for t in mesh['triangles']:
        a,b,c=[positions[i] for i in t['node_ids']]
        area=abs((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))/2
        for i in t['node_ids']: areas[i]+=area/3
    weights=[areas[i] for i in ids]
    regions=spatial_regions(mesh,contract)
    output=[]
    for bias in contract['spatial']['biases_V']:
        branch='reverse' if bias<0 else 'forward'
        tdr=d/f'{branch}_{tag(bias)}_des.tdr'
        export=d/'fields'/f'{branch}_{tag(bias)}'
        export.mkdir(parents=True,exist_ok=True)
        export_identity={'tdr_sha256':sha(tdr),'importer_sha256':sha(REPO/'build-release/sentaurus_import.exe')}
        provenance=export/'export_provenance.json'
        if not provenance.exists() or read_json(provenance)!=export_identity:
            if run_cmd([REPO/'build-release/sentaurus_import.exe','--tdr',tdr,'--export-dir',export],REPO,export/'import.log'):
                raise ValueError(f'field TDR export failed: {tdr}')
            write_json(provenance,export_identity)
        exported_nodes=rows(export/'nodes.csv')
        if len(exported_nodes)!=len(ids) or len({int(n['id']) for n in exported_nodes})!=len(ids):
            raise ValueError('state TDR node cardinality mismatch')
        for n in exported_nodes:
            x,y=positions[int(n['id'])]
            if max(abs(float(n['x_um'])-x),abs(float(n['y_um'])-y))>contract['input']['coordinate_tolerance_um']:
                raise ValueError('state TDR coordinates differ from solver mesh')
        # Match by the bias in the VTK title/filename, not by an array ordinal.
        vtk_candidates=[]
        artifacts=read_json(d/'spatial_artifacts.json') if (d/'spatial_artifacts.json').exists() else {}
        if f'{bias:g}' in artifacts:
            record=artifacts[f'{bias:g}']
            if sha(d/record['state_file'])!=record['state_sha256']: raise ValueError('postprocessed state fingerprint mismatch')
            vtk_candidates=[d/record['vtk']]
        else:
            for p in (d/branch).glob('*.vtk'):
                match=re.search(r'_(-?[\d.eE+]+)V\.vtk$',p.name)
                if match and abs(float(match[1])-bias)<=1e-9: vtk_candidates.append(p)
        if len(vtk_candidates)!=1: raise ValueError(f'expected one VTK with native bias {bias}, found {len(vtk_candidates)}')
        vtk=vtk_candidates[0]; count,scalars,vectors=read_vtk_point_data(vtk)
        if count!=len(ids): raise ValueError('VTK node count mismatch')
        validate_vtk_coordinates(vtk,[positions[i] for i in ids],contract['input']['coordinate_tolerance_um'])
        mapping=[('ElectrostaticPotential','Potential','potential'),('eQuasiFermiPotential','ElectronQuasiFermi','qf'),
                 ('hQuasiFermiPotential','HoleQuasiFermi','qf'),('eDensity','Electrons','carrier'),('hDensity','Holes','carrier'),
                 ('EffectiveIntrinsicDensity','EffectiveIntrinsicDensity','carrier')]
        metrics=[]
        for sent,vela,kind in mapping:
            r=node_field(export/'fields'/f'{sent}_region0.csv',ids)
            if vela not in scalars: raise ValueError(f'missing VTK scalar {vela}; available {list(scalars)}')
            c=scalars[vela]
            active=[i for i in range(count) if kind!='carrier' or r[i]>=contract['spatial']['carrier_active_reference_min_cm3']]
            errors=[math.log10(max(c[i],1e-300)/r[i]) if kind=='carrier' else c[i]-r[i] for i in active]
            rmse=math.sqrt(sum(e*e for e in errors)/len(errors)); p95=percentile([abs(e) for e in errors],0.95)
            lim=contract['spatial']
            passed=(rmse<=lim['carrier_log_rmse_dex'] and p95<=lim['carrier_log_p95_dex']) if kind=='carrier' else (p95<=lim['qf_p95_V'] if kind=='qf' else rmse<=lim['potential_rmse_V'] and p95<=lim['potential_p95_V'])
            metrics.append({'field':sent,'active_nodes':len(active),'rmse':rmse,'p95':p95,'max_diagnostic':max(abs(e) for e in errors),'pass':passed})
        transport=[]
        for sent,vela,limit,vector in [('ElectricField','NodeElectricField_AreaAverageVector','field_relative_l2',True),
                                     ('eMobility','ElectronMobilityCm2PerVs','mobility_relative_l2',False),
                                     ('hMobility','HoleMobilityCm2PerVs','mobility_relative_l2',False),
                                     ('srhRecombination','SRHRecombinationCm3PerS','srh_relative_l2',False)]:
            r=node_field(export/'fields'/f'{sent}_region0.csv',ids,2 if vector else 1)
            c=vectors[vela] if vector else scalars[vela]
            differences=[sum((c[i][j]-r[i][j])**2 for j in range(2)) if vector else (c[i]-r[i])**2 for i in range(count)]
            magnitudes=[sum(t*t for t in v) if vector else v*v for v in r]
            numerator=sum(w*e for w,e in zip(weights,differences)); denominator=sum(w*e for w,e in zip(weights,magnitudes))
            relative=math.sqrt(numerator/max(denominator,1e-300))
            item={'field':sent,'relative_l2':relative,'pass':relative<=contract['spatial'][limit]}
            # Endpoint spikes get area-weighted region statistics; node maxima are not gates.
            for label,mask in regions.items():
                selected=[i for i in range(count) if mask[i]]
                area=sum(weights[i] for i in selected)
                region={'nodes':len(selected),'area_um2':area,'weighted_error_l2':math.sqrt(sum(weights[i]*differences[i] for i in selected))}
                for tool,field in [('reference',r),('candidate',c)]:
                    magnitude=[math.sqrt(sum(v*v for v in field[i][:2])) if vector else abs(field[i]) for i in selected]
                    region[tool+'_rms']=math.sqrt(sum(weights[i]*v*v for i,v in zip(selected,magnitude))/area)
                    region[tool+'_p95']=percentile(magnitude,.95)
                    region[tool+'_max_diagnostic']=max(magnitude)
                    if sent=='srhRecombination': region[tool+'_integral_A_per_um']=sum(weights[i]*field[i] for i in selected)*1e-12*1.602176634e-19
                item[label]=region
            if sent=='srhRecombination':
                ri=sum(w*v for w,v in zip(weights,r))*1e-12*1.602176634e-19
                ci=sum(w*v for w,v in zip(weights,c))*1e-12*1.602176634e-19
                floor=contract['spatial']['srh_integral_current_floor_A_per_um']
                item.update(reference_integral_A_per_um=ri,candidate_integral_A_per_um=ci,
                            integral_pass=abs(ci-ri)<=max(floor,abs(ri)*contract['spatial']['srh_integral_relative']))
                # Near equilibrium the signed source and its L2 denominator both vanish.
                if abs(ri)<=floor and abs(ci)<=floor: item['pass']=True; item['low_source_absolute_gate']=True
                item['pass']=item['pass'] and item['integral_pass']
            transport.append(item)
        token=('m' if bias<0 else '')+f'{abs(bias):.6f}'.replace('.','p')
        probe=read_json(d/f'vela_{branch}.json'); probe.pop('sweep')
        probe.update(simulation_type='sg_edge_flux_probe',state_file=f'{branch}/accepted_bias_{token}.csv',output_csv=f'{branch}/edges_{token}.csv')
        for c in probe['contacts']: c['bias']=bias if c['name']=='Anode' else 0.0
        write_json(d/'probe.json',probe)
        edge_path=d/probe['output_csv']; edge_provenance=edge_path.with_suffix('.provenance.json')
        identity={'state_sha256':sha(d/probe['state_file']),'config_sha256':sha(d/'probe.json'),
                  'runner_sha256':sha(REPO/'build-release/vela_example_runner.exe')}
        prior=read_json(edge_provenance) if edge_provenance.exists() else {}
        cached=(prior.get('inputs')==identity and prior.get('returncode')==0 and edge_path.exists()
                and prior.get('edge_sha256')==sha(edge_path))
        rc=0 if cached else run_cmd([REPO/'build-release/vela_example_runner.exe','--config',d/'probe.json'],REPO,d/'probe.log')
        if rc: raise ValueError('own-state SG flux diagnostic failed; stale output cannot be accepted')
        if not cached: write_json(edge_provenance,{'inputs':identity,'returncode':rc,'edge_sha256':sha(edge_path)})
        edges=rows(edge_path)
        if not edges: raise ValueError('empty edge flux diagnostics')
        sections=[{'x_um':x,**section_current(edges,axis='x',cut_um=x)} for x in contract['conservation']['sections_x_um']]
        terminal=exact_row(rows(d/f'{branch}_terminal_balance.csv'),bias,contact='Cathode')
        # SDevice and Vela terminals report conventional current supplied by the
        # external circuit. A +x section exits the right contact, so negate it.
        target=-float(terminal['current_total_A_per_um'])
        tolerance=max(contract['conservation']['absolute_floor_A_per_um'],abs(target)*contract['conservation']['section_relative_max'])
        for sec in sections: sec['pass']=abs(sec['total_A_per_um']-target)<=tolerance
        source_closure=[]
        source=next(t for t in transport if t['field']=='srhRecombination')
        sr=exact_row(plt_rows(d/f'{branch}.plt'),bias,'Anode OuterVoltage')
        vc=[exact_row(rows(d/f'{branch}_terminal_balance.csv'),bias,contact=name) for name in ['Anode','Cathode']]
        for tool in ['reference','candidate']:
            integral=source[tool+'_integral_A_per_um']
            for component,key,sign in [('electron','eCurrent',1),('hole','hCurrent',-1)]:
                currents=([float(sr[f'{name} {key}']) for name in ['Anode','Cathode']] if tool=='reference'
                          else [float(c[f'current_{component}_A_per_um'])*sign for c in vc])
                residual=abs(sum(currents)+sign*integral)
                scale=max(abs(integral),*(abs(v) for v in currents),contract['conservation']['absolute_floor_A_per_um'])
                passed=residual<=max(contract['conservation']['absolute_floor_A_per_um'],scale*contract['conservation']['continuity_relative_max'])
                source_closure.append({'tool':tool,'component':component,'residual_A_per_um':residual,'relative':residual/scale,'pass':passed})
        current_vectors=[]
        for sent,vela in [('eCurrentDensity','SentaurusElectronCurrentDensityVector'),('hCurrentDensity','SentaurusHoleCurrentDensityVector'),('TotalCurrentDensity','SentaurusTotalCurrentDensityVector')]:
            reference=node_field(export/'fields'/f'{sent}_region0.csv',ids,2)
            candidate=vectors[vela]
            norm=sum(weights[i]*sum(v*v for v in reference[i]) for i in range(count))
            diff=sum(weights[i]*sum((candidate[i][j]-reference[i][j])**2 for j in range(2)) for i in range(count))
            lateral=[sum(weights[i]*field[i][1]**2 for i in range(count)) for field in [reference,candidate]]
            longitudinal=[sum(weights[i]*field[i][0]**2 for i in range(count)) for field in [reference,candidate]]
            current_vectors.append({'field':sent,'relative_l2':math.sqrt(diff/max(norm,1e-300)),
                                    'reference_Jy_l2':math.sqrt(lateral[0]),'candidate_Jy_l2':math.sqrt(lateral[1]),
                                    'reference_Jx_l2':math.sqrt(longitudinal[0]),'candidate_Jx_l2':math.sqrt(longitudinal[1]),
                                    'reference_Jy_over_Jx_l2':math.sqrt(lateral[0]/max(longitudinal[0],1e-300)),
                                    'candidate_Jy_over_Jx_l2':math.sqrt(lateral[1]/max(longitudinal[1],1e-300)),
                                    'status':'diagnostic_reconstruction; conservative SG sections are authoritative'})
        height=read_json(SPEC)['geometry_um']['height']
        transverse=[{'y_um':height*f/8,**section_current(edges,axis='y',cut_um=height*f/8)} for f in range(1,8)]
        output.append({'bias_V':bias,'vtk':vtk.relative_to(d).as_posix(),'metrics':metrics,'transport':transport,'sections':sections,
                       'transverse_sg_sections':transverse,'source_closure':source_closure,'current_vectors':current_vectors,'diagnostic_probe_exitcode':rc})
    return {'status':'pass' if all(m['pass'] for p in output for m in p['metrics']) else 'fail','points':output,
            'transport_source_and_section_status':'pass' if all(m['pass'] for p in output for m in p['transport']+p['sections']+p['source_closure']) else 'fail'}


def audit_imported_mesh(d, entry, root):
    for name,digest in entry['input_sha256'].items():
        if sha(d/name)!=digest: raise ValueError(f'input modified after manifest: {name}')
    mesh=read_json(d/'inputs/mesh.json'); positions={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    imported=read_json(d/'input_audit.json')
    for file,key in [('inputs/mesh.json','mesh_sha256'),('inputs/doping.csv','doping_sha256'),('pn2d_msh.tdr','tdr_sha256')]:
        if sha(d/file)!=imported[key]: raise ValueError(f'imported input changed after independent solve: {file}')
    if len(positions)!=len(mesh['nodes']): raise ValueError('duplicate mesh node IDs')
    if not all(abs(v-t)<1e-10 for v,t in zip([min(x for x,y in positions.values()),max(x for x,y in positions.values()),min(y for x,y in positions.values()),max(y for x,y in positions.values())],[0,2,0,0.5])):
        raise ValueError('geometry mismatch')
    contact_results=[]
    for c in mesh['contacts']:
        x=0 if c['name']=='Anode' else 2
        length=entry['parameters']['anode_fraction']*0.5 if c['name']=='Anode' else 0.5
        points=[positions[i] for i in c['node_ids']]
        if not points or any(abs(px-x)>1e-10 for px,py in points): raise ValueError('contact on wrong boundary')
        low,high=min(y for _,y in points),max(y for _,y in points)
        if max(abs(low-(0.5-length)/2),abs(high-(0.5+length)/2))>1e-10: raise ValueError('contact extent mismatch')
        contact_results.append({'name':c['name'],'length_um':length,'nodes':len(points)})
    if {c['name'] for c in contact_results}!={'Anode','Cathode'}: raise ValueError('unexpected or missing electrodes')
    doping=rows(d/'inputs/doping.csv')
    if len(doping)!=len(positions) or {int(r['node_id']) for r in doping}!=set(positions): raise ValueError('doping node coverage mismatch')
    for r in doping:
        x,y=positions[int(r['node_id'])]
        donors=float(r['donors_cm3']); acceptors=float(r['acceptors_cm3'])
        # Inclusive source windows overlap on the junction plane. Preserve and
        # verify both reported species there, including asymmetric variants.
        expected=((entry['parameters']['ND_cm3'],entry['parameters']['NA_cm3']) if abs(x-1.0)<1e-10 else
                  ((0,entry['parameters']['NA_cm3']) if x<1 else (entry['parameters']['ND_cm3'],0)))
        if any(abs(a-b)>max(1.0,abs(b)*1e-10) for a,b in zip([donors,acceptors],expected)): raise ValueError('nodal doping plateau mismatch')
    same_mesh=None
    baseline_path=root/'P0'/entry['mesh']/'inputs/mesh.json'
    same_topology=None
    if baseline_path.exists():
        baseline=read_json(baseline_path)
        same_topology=mesh['nodes']==baseline['nodes'] and mesh['triangles']==baseline['triangles']
    if entry['id'] in ['P1','P2','P3']:
        same_mesh=sha(d/'pn2d_msh.tdr')==sha(root/'P0'/entry['mesh']/'pn2d_msh.tdr')
        if not same_mesh: raise ValueError('stage A mesh changed')
    return {'status':'pass','nodes':len(positions),'triangles':len(mesh['triangles']),'contacts':contact_results,
            'same_baseline_tdr':same_mesh,'same_baseline_coordinates_and_connectivity':same_topology,
            'junction_plane_doping_verified':True,'mesh_sha256':sha(d/'inputs/mesh.json'),'node_doping_sha256':sha(d/'inputs/doping.csv')}


def audit_input(d, entry, root):
    run=read_json(d/'sentaurus_run.json')
    if run.get('status')=='running': raise FileNotFoundError('Sentaurus run is in progress')
    if run['returncode']!=0 or run['fetch_returncode']!=0: raise ValueError('Sentaurus run or retrieval failed')
    if run['version']!=read_json(SPEC)['sentaurus_version']: raise ValueError('mixed Sentaurus versions')
    mesh_result=audit_imported_mesh(d,entry,root)
    for branch in ['forward','reverse']:
        cfg=read_json(d/f'vela_{branch}.json')
        if 'initial_state_file' in cfg['sweep']: raise ValueError('external initialization is forbidden')
        record=read_json(d/f'vela_{branch}_run.json')
        if record.get('status')=='running': raise FileNotFoundError('Vela run is in progress')
        if 'config_sha256' in record and record['config_sha256']!=sha(d/f'vela_{branch}.json'):
            raise ValueError('Vela configuration changed after execution')
        if 'runtime_config' in record:
            runtime=read_json(d/record['runtime_config'])
            if sha(d/record['runtime_config'])!=record['runtime_config_sha256']: raise ValueError('runtime configuration changed')
            expected=json.loads(json.dumps(cfg))
            if record.get('compact_spatial'): expected['sweep']['write_vtk']=False
            if runtime!=expected: raise ValueError('runtime override changes more than spatial output')
    return mesh_result


def mesh_gate(identifier, result_map, contract, pair=None):
    a,b=[result_map.get((identifier,m)) for m in (pair or contract['mesh']['required_pair'])]
    if not a or not b or 'metrics' not in a['terminal'] or 'metrics' not in b['terminal']:
        return {'status':'not_run'}
    lookup={(r['branch'],r['bias_V'],r['contact'],r['component']):r for r in b['terminal']['metrics']}
    metrics=[]
    for r in a['terminal']['metrics']:
        fine=lookup[(r['branch'],r['bias_V'],r['contact'],r['component'])]
        for tool in ['reference','candidate']:
            coarse_value=r[tool+'_A_per_um']; fine_value=fine[tool+'_A_per_um']
            absolute=abs(coarse_value-fine_value); relative=absolute/max(abs(fine_value),contract['mesh']['current_absolute_floor_A_per_um'])
            passed=absolute<=contract['mesh']['current_absolute_floor_A_per_um'] or relative<=contract['mesh']['current_relative_max']
            metrics.append({'branch':r['branch'],'bias_V':r['bias_V'],'contact':r['contact'],'component':r['component'],'tool':tool,'relative_error':relative,'absolute_error_A_per_um':absolute,'pass':passed})
    source=[]; regional=[]
    if 'points' not in a['spatial'] or 'points' not in b['spatial']: return {'status':'not_run','terminal_metrics':metrics}
    for p,q in zip(a['spatial']['points'],b['spatial']['points'],strict=True):
        for coarse,fine in zip(p['transport'],q['transport'],strict=True):
            for region in ['contact_edges','junction']:
                for tool in ['reference','candidate']:
                    for metric in ['rms','p95','max_diagnostic','integral_A_per_um']:
                        key=tool+'_'+metric
                        if key not in coarse.get(region,{}) or key not in fine.get(region,{}): continue
                        x,y=coarse[region][key],fine[region][key]
                        regional.append({'bias_V':p['bias_V'],'field':coarse['field'],'region':region,'tool':tool,
                                         'metric':metric,'coarse':x,'fine':y,'absolute_change':abs(x-y),
                                         'relative_change_diagnostic':abs(x-y)/max(abs(y),1e-300)})
        r=next(t for t in p['transport'] if t['field']=='srhRecombination')
        s=next(t for t in q['transport'] if t['field']=='srhRecombination')
        for tool in ['reference','candidate']:
            x,y=r[tool+'_integral_A_per_um'],s[tool+'_integral_A_per_um']
            passed=abs(x-y)<=max(contract['mesh']['current_absolute_floor_A_per_um'],abs(y)*contract['mesh']['integral_relative_max'])
            source.append({'bias_V':p['bias_V'],'tool':tool,'coarse_A_per_um':x,'fine_A_per_um':y,'pass':passed})
    return {'status':'pass' if all(r['pass'] for r in metrics+source) else 'fail','terminal_metrics':metrics,'srh_integrals':source,
            'regional_mesh_diagnostics':regional}


def model_effect(identifier, mesh, result_map, contract):
    base=result_map.get(('P0',mesh)); variant=result_map.get((identifier,mesh))
    if not base or not variant or 'metrics' not in base['terminal'] or 'metrics' not in variant['terminal']: return {'status':'not_run'}
    baseline={(r['branch'],r['bias_V'],r['contact'],r['component']):r for r in base['terminal']['metrics']}
    result=[]
    for r in variant['terminal']['metrics']:
        b=baseline[(r['branch'],r['bias_V'],r['contact'],r['component'])]
        vals=[abs(v[k+'_A_per_um']) for v in [r,b] for k in ['reference','candidate']]
        if min(vals)<=contract['model_effect']['absolute_floor_A_per_um']: continue
        sr=math.log10(abs(r['reference_A_per_um']/b['reference_A_per_um']))
        vr=math.log10(abs(r['candidate_A_per_um']/b['candidate_A_per_um']))
        result.append({'branch':r['branch'],'bias_V':r['bias_V'],'contact':r['contact'],'component':r['component'],
                       'reference_delta_dex':sr,'candidate_delta_dex':vr,'error_dex':abs(sr-vr),'pass':abs(sr-vr)<=contract['model_effect']['delta_log_current_max_dex']})
    return {'status':'pass' if result and all(r['pass'] for r in result) else 'fail','metrics':result}


def compare(root):
    contract=read_json(root/'contract.json'); campaign=read_json(root/'campaign.json'); results=[]
    tool_identity={name:sha(path) if path.exists() else None for name,path in
                   [('diagnostic_runner',REPO/'build-release/vela_example_runner.exe'),('importer',REPO/'build-release/sentaurus_import.exe')]}
    if sha(root/'contract.json')!=campaign['contract_sha256']: raise ValueError('contract fingerprint mismatch')
    for e in campaign['entries']:
        d=root/e['directory']; result={'id':e['id'],'mesh':e['mesh'],'status':'not_run'}
        evidence={}
        names=['sentaurus_run.json','vela_forward_run.json','vela_reverse_run.json','pn2d_msh.tdr','forward.plt','reverse.plt',
               'forward.csv','reverse.csv','forward_terminal_balance.csv','reverse_terminal_balance.csv','spatial_artifacts.json',
               'input_audit.json','inputs/mesh.json','inputs/doping.csv','manifest.json','vela_forward_runtime.json','vela_reverse_runtime.json',*e['input_sha256']]
        paths=[d/name for name in names]
        if (d/'spatial_artifacts.json').exists():
            for record in read_json(d/'spatial_artifacts.json').values(): paths.extend([d/record['vtk'],d/record['state_file']])
        for bias in contract['spatial']['biases_V']:
            branch='reverse' if bias<0 else 'forward'
            paths.append(d/f'{branch}_{tag(bias)}_des.tdr')
            paths.extend((d/branch).glob(f'*_{bias:g}V.vtk'))
        for path in paths:
            evidence[path.relative_to(d).as_posix()]=sha(path) if path.exists() else None
        fingerprint=hashlib.sha256(json.dumps({'inputs':e['input_sha256'],'evidence':evidence,'contract':sha(root/'contract.json'),
                                               'comparator':sha(Path(__file__)),
                                               **tool_identity},sort_keys=True).encode()).hexdigest()
        previous=d/'comparison.json'
        if previous.exists():
            cached=read_json(previous)
            if cached.get('evidence_fingerprint')==fingerprint:
                results.append(cached)
                continue
        result['evidence_fingerprint']=fingerprint
        result['evidence_sha256']=evidence
        try: result['input']=audit_input(d,e,root)
        except OSError as ex: result['input']={'status':'not_run','reason':str(ex)}
        except (ValueError,KeyError) as ex: result['input']={'status':'fail','reason':str(ex)}
        for name,fn in [('terminal',terminals),('spatial',spatial)]:
            if 'in progress' in result['input'].get('reason',''):
                result[name]={'status':'not_run','reason':'run in progress; incomplete output is not comparison evidence'}
                continue
            try: result[name]=fn(d,contract)
            except (OSError,ValueError,KeyError,StopIteration) as ex: result[name]={'status':'not_run','reason':str(ex)}
        result['status']='pass' if all(result[n]['status']=='pass' for n in ['input','terminal','spatial']) and result['terminal'].get('convergence_kcl_status')=='pass' and result['spatial'].get('transport_source_and_section_status')=='pass' else 'incomplete'
        if any(result[n]['status']=='fail' for n in ['input','terminal','spatial']) or result['spatial'].get('transport_source_and_section_status')=='fail' or result['terminal'].get('convergence_kcl_status')=='fail': result['status']='fail'
        for branch in ['forward','reverse']:
            run=d/f'vela_{branch}_run.json'
            if run.exists() and read_json(run)['returncode']:
                result['status']='fail'
                result[branch+'_failure']='Vela did not complete the independent scan; see run record and nonlinear trace'
        write_json(d/'comparison.json',result); results.append(result)
    result_map={(r['id'],r['mesh']):r for r in results}
    mesh={e['id']:mesh_gate(e['id'],result_map,contract) for e in read_json(SPEC)['experiments']}
    refinement=read_json(SPEC.parent/'mesh_refinement.json')
    if refinement['frozen_contract_sha256']!=sha(root/'contract.json'):
        raise ValueError('mesh continuation must use the original frozen thresholds')
    pair=refinement['additional_pair']
    refined_mesh={e['id']:mesh_gate(e['id'],result_map,contract,pair) for e in read_json(SPEC)['experiments']}
    historical_mesh={ '/'.join(p):{e['id']:mesh_gate(e['id'],result_map,contract,p) for e in read_json(SPEC)['experiments']} for p in refinement.get('previous_additional_pairs',[])}
    contact=refinement['contact_edge_control']
    contact_mesh={i:mesh_gate(i,result_map,contract,contact['pair']) for i in contact['ids']}
    grids=[m for m in read_json(SPEC)['meshes'] if m!='original']
    effects={i:{m:model_effect(i,m,result_map,contract) for m in grids} for i in ['P1','P2','P3']}
    stage_a=all(result_map.get((i,m),{}).get('status')=='pass' for i in ['P0','P1','P2','P3'] for m in pair) and all(refined_mesh[i]['status']=='pass' for i in ['P0','P1','P2','P3']) and all(effects[i][m]['status']=='pass' for i in effects for m in pair)
    required={('P0','original')}|{(e['id'],m) for e in read_json(SPEC)['experiments'] for m in pair}
    required|={(i,m) for i in contact['ids'] for m in contact['pair']}
    for r in results: r['required_for_final_qualification']=(r['id'],r['mesh']) in required
    report={'schema':'vela.pn2d.variants.results.v1','contract_sha256':sha(root/'contract.json'),'results':results,
            'overall_status':'pass_on_refined_mesh' if all(result_map[k]['status']=='pass' for k in required) and all(m['status']=='pass' for m in list(refined_mesh.values())+list(contact_mesh.values())) and stage_a else 'incomplete',
            'mesh':mesh,'historical_refinement_mesh':historical_mesh,'contact_edge_mesh':contact_mesh,'refined_mesh':refined_mesh,'refinement':refinement,'model_effects':effects,'stage_A_accepted':stage_a}
    write_json(root/'results.json',report)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    args=parser.parse_args(); compare(args.root.resolve())
