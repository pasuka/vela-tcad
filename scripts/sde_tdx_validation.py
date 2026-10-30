"""Validation-only mapping of regionalized TDX fields back to the input mesh."""
import math
from collections import Counter,defaultdict
from compare_sentaurus_tdr_tdx import parse_grd
from vela_sdevice.mesh_import import _edge_nodes
from audit_sde_native_corpus import read_native


def validate_pair(bundle,path,tolerance=1e-6):
    native=read_native(path);grid=parse_grd(path.with_suffix('.grd'))
    points=bundle['nodes'];returned=native['vertices']
    if len(points)!=len(returned):raise ValueError('TDX changed node count')
    bins=defaultdict(list);scale=1e-10
    for p in points:bins[tuple(round(p[k]/scale) for k in ('x','y'))].append(p['vela_id'])
    mapping={}
    for i,p in enumerate(returned):
        if math.dist(p,(points[i]['x'],points[i]['y']))<=1e-12:
            mapping[i]=i;continue
        key=tuple(round(q/scale) for q in p);matches=[]
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                matches += [n for n in bins[key[0]+dx,key[1]+dy]
                            if math.dist(p,(points[n]['x'],points[n]['y']))<=1e-12]
        if len(matches)!=1:raise ValueError('TDX changed or ambiguously mapped node coordinates')
        mapping[i]=matches[0]
    if len(set(mapping.values()))!=len(points):raise ValueError('TDX node mapping is not one-to-one')
    names={r['id']:r['name'] for r in bundle['regions']}
    original=Counter((names[t['region_id']],tuple(sorted(t['node_ids']))) for t in bundle['triangles'])
    actual=Counter((name,tuple(sorted(mapping[n] for n in ids))) for name,ids in native['cells'])
    if actual!=original:raise ValueError('TDX changed region topology')
    aliases=bundle['provenance']['material_aliases']
    material_names={m['id']:aliases.get(m['name'],m['name']) for m in bundle['materials']}
    for region in bundle['regions']:
        expected={material_names[t['material_id']] for t in bundle['triangles'] if t['region_id']==region['id']}
        raw=native['regions'][region['name']]['material']
        if expected!={aliases.get(raw,raw)}:raise ValueError('TDX changed material identity')
    for kind in ('contacts','interfaces'):
        observed={}
        for name,record in grid['region_blocks'].items():
            if record['material']!=('Contact' if kind=='contacts' else 'Interface'):continue
            edges=[]
            for i in grid['region_elements'][name]:
                row=grid['elements'][i]
                if row[0]!=1:raise ValueError('TDX boundary became non-linear')
                ids=row[1:] if len(row)==3 else _edge_nodes(row[1],grid['edges'])
                edges.append(tuple(sorted(mapping[int(n)] for n in ids)))
            observed[name]=Counter(edges)
        expected={g['name']:Counter(tuple(sorted(line['node_ids'])) for line in bundle['lines']
                    if g['physical_tag'] in line['physical_curve_tags']) for g in bundle[kind]}
        if observed!=expected:raise ValueError('TDX changed '+kind+' edges or names')
    expected={k:v['values'] for k,v in bundle['fields'].items() if k!='species'}
    expected.update({k:v['values'] for k,v in bundle['fields'].get('species',{}).items()})
    seen={k:set() for k in expected};maximum={k:0. for k in expected};entries=Counter()
    for region,fields in native['fields'].items():
        for field,values in fields.items():
            if field not in expected:continue
            for node,value in values.items():
                target=mapping[node];reference=expected[field][target]
                error=abs(value-reference)/max(abs(reference),1.)
                maximum[field]=max(maximum[field],error);seen[field].add(target);entries[field]+=1
                if not math.isfinite(value) or error>tolerance:
                    raise ValueError('TDX changed stored '+field+' beyond the frozen tolerance')
    implicit={}
    for field,values in expected.items():
        missing=set(range(len(points)))-seen[field]
        if any(values[n]!=0. for n in missing):raise ValueError('TDX lost nonzero stored '+field)
        implicit[field]=len(missing)
    return {'status':'passed','nodes':len(points),'triangles':sum(actual.values()),
            'maximum_stored_field_relative_errors':maximum,'implicit_zero_nodes':implicit,
            'regional_field_entries':dict(entries),'validation_only_regional_mapping':len(native['regions'])>1,
            'net_policy':'compare stored net directly; never recompute it from rounded carrier text'}
