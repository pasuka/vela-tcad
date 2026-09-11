"""Prepare the explicit D4 geometry/config bundle from qualified local evidence."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from audit_templates_ldmos_averagebox_node4492 import (
    parse_debug_block,measure_in_tdr_vertex_order,coefficient_edge)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1];out=args.output.resolve()
    if out.exists():raise ValueError('Refusing to overwrite D4 preparation')
    out.mkdir(parents=True)
    read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    base=root/'reference_staging/templates_ldmos_followup_20260911'
    bundle=read(root/'reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_inputs.json')
    config=read(root/bundle['template'])
    mesh_path=Path(config['mesh_file'].replace('@workspace',str(root)));mesh=read(mesh_path)
    debug=base/'ialmob_difference_location/native/geometry/MeasureCoefficients.debug'
    measure=parse_debug_block(debug,'Measure');coeff=parse_debug_block(debug,'Coefficients')
    material={r['id']:r['material'] for r in mesh['regions']}
    xy={n['id']:(n['x']*1e-6,n['y']*1e-6) for n in mesh['nodes']}
    cells=[];edge_sums=defaultdict(float)
    for c in mesh['triangles']:
        if material[c['region_id']]!='Si':continue
        by_edge={tuple(sorted(coefficient_edge(c['node_ids'],k)[:2])):v for k,v in enumerate(coeff[c['id']])}
        weights=[]
        for k in range(3):
            a,b=c['node_ids'][k],c['node_ids'][(k+1)%3];key=tuple(sorted((a,b)))
            value=by_edge[key]
            if not math.isfinite(value) or value<0:raise ValueError('Invalid coefficient')
            edge_sums[key]+=value*math.dist(xy[a],xy[b]);weights.append(value)
        cells.append(dict(cell_id=c['id'],node_ids=c['node_ids'],vertex_measure_m2=[v*1e-12 for v in measure_in_tdr_vertex_order(measure[c['id']])],edge_coefficients=weights))
    couples=Path(config['mesh_geometry']['external_averagebox_couples_file'].replace('@workspace',str(root)))
    with couples.open(newline='') as f:
        expected={tuple(sorted((int(r['node0']),int(r['node1'])))):float(r['couple_m']) for r in csv.DictReader(f)}
    mismatch=[(key,v,expected.get(key)) for key,v in edge_sums.items() if key not in expected or not math.isclose(v,expected[key],rel_tol=1e-10,abs_tol=0.)]
    if mismatch:raise ValueError(mismatch[:5])
    geometry=out/'geometry.json';geometry.write_text(json.dumps(dict(schema='vela.ialmob.transport_geometry.v1',node_count=len(mesh['nodes']),cell_count=len(mesh['triangles']),cells=cells)))
    ial=dict(geometry_file='@workspace/'+geometry.relative_to(root).as_posix(),effective_electrodes=['source','drain','substrate'],crystal_x=[1,0,0],crystal_y=[0,1,0],high_field=True,reference_density_m3=1e18)
    params_sources=[]
    for prefix,carrier in [('e','electron'),('h','hole')]:
        groups={}
        for family in ('100','110'):
            p=base/f'ialmob_closure_r2/raw_distance_gradient_replay/{prefix}_absolute_{family}_input.json'
            groups[family]=read(p)['parameters_cm'];params_sources.append(p)
        ial[carrier+'_parameters_cm']=groups
    config['_comment']='Experimental coupled D4 profile; qualification is tracked separately from D5.'
    config['solver']['mobility'].update(model='ialmob',ialmob=ial,contact_electric_field_fallback=False)
    profile=out/'config.json';profile.write_text(json.dumps(config,indent=2))
    refs={str(g):str(Path(v.replace('D5-no-IALMob','D4-classical')).as_posix()) for g,v in bundle['references'].items()}
    dependencies={p:h for p,h in bundle['files'].items() if p not in bundle['references'].values() and p!=bundle['template']}
    for p in [geometry,profile,debug,root/'scripts/prepare_templates_ldmos_ialmob_transport.py',*params_sources,*[root/v for v in refs.values()]]:
        dependencies[p.relative_to(root).as_posix()]=sha(p)
    bundle.update(schema='vela.templates_ldmos.linked_d4_inputs.v1',template=profile.relative_to(root).as_posix(),references=refs,files=dependencies,
        note='D5 zero states are initial guesses only; D4 initialize must reclose them under live IALMob before any drain advance. Original numerical and Id-Vd scoring gates retained.')
    (out/'inputs.json').write_text(json.dumps(bundle,indent=2))
    (out/'geometry_audit.json').write_text(json.dumps(dict(silicon_cells=len(cells),transport_edges=len(edge_sums),edge_coupling_mismatches=0,source_debug_sha256=sha(debug),geometry_sha256=sha(geometry)),indent=2))
    print(out/'inputs.json')


if __name__=='__main__':main()
