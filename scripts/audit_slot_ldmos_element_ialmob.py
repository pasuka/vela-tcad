"""Replay all SLOT native Si elements at a frozen Sentaurus state in SI units."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
try:
    from prepare_slot_ldmos_ialmob_probe import nearest_segment, read_field
except ModuleNotFoundError:
    from scripts.prepare_slot_ldmos_ialmob_probe import nearest_segment, read_field


def validate_export_mesh(mesh: dict, export: Path) -> None:
    """Require the same coordinate frame and topology before comparing fields."""
    with (export/'nodes.csv').open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    nodes={int(n['id']):(float(n['x_um']),float(n['y_um'])) for n in rows}
    expected={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    if len(nodes)!=len(rows) or nodes!=expected:
        raise ValueError('Native and local node coordinates differ')
    with (export/'elements.csv').open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    cells={int(c['id']):tuple(int(c['node'+str(k)]) for k in range(3)) for c in rows}
    expected={c['id']:tuple(c['node_ids']) for c in mesh['triangles']}
    if len(cells)!=len(rows) or cells!=expected:
        raise ValueError('Native and local element topology differs')


def audit(mesh_path: Path, geometry_path: Path, parameters_path: Path,
          export: Path, probe: Path, output: Path) -> dict:
    if output.exists(): raise ValueError("Refusing to overwrite an element audit")
    mesh = json.loads(mesh_path.read_text())
    validate_export_mesh(mesh,export)
    geometry = json.loads(geometry_path.read_text())
    cells_by_id={c['id']:c['node_ids'] for c in mesh['triangles']}
    if any(c['node_ids']!=cells_by_id.get(c['cell_id']) for c in geometry['cells']):
        raise ValueError('Mobility geometry topology differs')
    params = json.loads(parameters_path.read_text())
    metadata=json.loads((export/'metadata.json').read_text())
    units={(f['name'],f['region']):f['unit'] for f in metadata['fields']}
    for name,unit in {'ElectrostaticPotential':'V','eQuasiFermiPotential':'V','hQuasiFermiPotential':'V',
                      'DonorConcentration':'cm^-3','AcceptorConcentration':'cm^-3',
                      'eDensity':'cm^-3','hDensity':'cm^-3','eMobility':'cm^2*V^-1*s^-1',
                      'hMobility':'cm^2*V^-1*s^-1','LatticeTemperature':'K'}.items():
        if units.get((name,0))!=unit: raise ValueError('Unexpected native unit: '+name)
    if any(abs(t[0]-300.)>1e-10 for t in read_field(export,'LatticeTemperature').values()):
        raise ValueError('This replay requires uniform 300 K')
    xy = {n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    cells = {c['id']:c for c in mesh['triangles']}
    materials = {r['id']:r['material'] for r in mesh['regions']}
    adjacent = defaultdict(list)
    for c in cells.values():
        ns = c['node_ids']
        for k in range(3): adjacent[tuple(sorted((ns[k],ns[(k+1)%3])))].append(c['id'])
    contacts, effective = set(), set()
    for c in mesh['contacts']:
        nodes = set(c['node_ids'])
        if c['name'] in ('source','drain','substrate'): effective.update(nodes)
        contacts.update(e for e,v in adjacent.items() if len(v)==1 and set(e)<=nodes)
    segments, tangents, priority = [], {}, {}
    for edge, ids in adjacent.items():
        si = [i for i in ids if materials[cells[i]['region_id']]=='Si']
        other = [i for i in ids if materials[cells[i]['region_id']]!='Si']
        if not si: continue
        a,b = [xy[n] for n in edge]
        if any(materials[cells[i]['region_id']]=='SiO2' for i in other): segments.append((*a,*b))
        if other or (len(ids)==1 and edge not in contacts):
            length = math.dist(a,b);t=((b[0]-a[0])/length,(b[1]-a[1])/length)
            for cid in si:
                material = bool(other)
                if cid in tangents:
                    if priority[cid] and not material: continue
                    if priority[cid]==material:
                        old=tangents[cid]
                        if abs(t[0]*old[1]-t[1]*old[0])>1e-12: raise ValueError('Ambiguous boundary tangent')
                        continue
                tangents[cid],priority[cid]=t,material
    fields = {name:read_field(export,name) for name in (
        'ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential',
        'DonorConcentration','AcceptorConcentration','eDensity','hDensity')}
    si_nodes=set(fields['eDensity'])
    if any(set(v)!=si_nodes for v in fields.values()): raise ValueError('Si state supports differ')
    distance={n:nearest_segment(xy[n],segments)[0]*1e-6 for n in si_nodes}
    elements=[]
    for c in geometry['cells']:
        cid, ns=c['cell_id'],c['node_ids']
        vertices=[]
        for n in ns:
            vertices.append(dict(orientation_family=100,potential_V=fields['ElectrostaticPotential'][n][0],
                electron_qf_V=fields['eQuasiFermiPotential'][n][0],hole_qf_V=fields['hQuasiFermiPotential'][n][0],
                donors_m3=fields['DonorConcentration'][n][0]*1e6,acceptors_m3=fields['AcceptorConcentration'][n][0]*1e6,
                electrons_m3=fields['eDensity'][n][0]*1e6,holes_m3=fields['hDensity'][n][0]*1e6,
                electron_response_m3_per_V=fields['eDensity'][n][0]*1e6/.0258519997864355,
                hole_response_m3_per_V=fields['hDensity'][n][0]*1e6/.0258519997864355))
        elements.append(dict(cell_id=cid,coordinates_m=[[v*1e-6 for v in xy[n]] for n in ns],
            interface_distance_m=[distance[n] for n in ns],vertex_measure_m2=c['vertex_measure_m2'],
            partial_boundary_layer=cid in tangents,boundary_tangent=tangents.get(cid,[1.,0.]),
            touches_effective_electrode=bool(set(ns)&effective),vertices=vertices))
    output.mkdir(parents=True)
    inp=output/'input.json';result=output/'result.json'
    # Original SLOT has no RefDens interpolation; D4's explicit 1e12 cm^-3
    # belongs to another physical configuration.
    inp.write_text(json.dumps(dict(temperature_K=300,reference_density_m3=0.,**params,elements_SI=elements)))
    with result.open('w') as stream:
        subprocess.run([str(probe.resolve()),str(inp.resolve())],stdout=stream,check=True)
    actual={r['cell_id']:r for r in json.loads(result.read_text())['elements']}
    report={}
    for prefix,carrier in [('e','electron'),('h','hole')]:
        with (export/'fields'/f'{prefix}Mobility_region0_cells.csv').open(newline='') as stream:
            native={int(r['cell_id']):float(r['component0']) for r in csv.DictReader(stream)}
        if actual.keys()!=native.keys(): raise ValueError('Element mobility supports differ')
        rows=[]
        for cid in actual:
            value=actual[cid][carrier+'_m2_per_Vs']*1e4
            if not math.isfinite(value) or value<=0 or not math.isfinite(native[cid]) or native[cid]<=0:
                raise ValueError('Nonfinite or nonpositive element mobility')
            rows.append(dict(cell_id=cid,native_cm2_V_s=native[cid],local_cm2_V_s=value,
                             relative_error=value/native[cid]-1))
        errors=sorted(abs(r['relative_error']) for r in rows)
        report[carrier]=dict(count=len(rows),median=errors[len(errors)//2],p95=errors[int(.95*(len(errors)-1))],maximum=errors[-1],
            largest=sorted(rows,key=lambda r:abs(r['relative_error']),reverse=True)[:10])
        with (output/f'{carrier}.csv').open('w',newline='') as stream:
            w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report['sources']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (mesh_path,geometry_path,parameters_path,probe,export/'metadata.json')}
    report['scope']='Same-state element mobility comparison; no self-consistent qualification'
    (output/'summary.json').write_text(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','geometry','parameters','export','probe','output'): p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(audit(a.mesh,a.geometry,a.parameters,a.export,a.probe,a.output),indent=2))
