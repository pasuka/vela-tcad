"""Local y=0 mesh reflection. No shell or vendor program is executed."""
import copy
import hashlib
import json


def reflect(half,operation):
    from sentaurus_gmsh_backend import _groups, validate_bundle
    if operation['axis']!='y' or operation['coordinate_um']!=0:
        raise ValueError('only y=0 reflection is supported')
    if min(n['y'] for n in half['nodes']) < -1e-8:
        raise ValueError('reflection input crosses its axis')
    contact_names={c['name'] for c in half['contacts']}
    if 'drain' not in contact_names or 'source' in contact_names:
        raise ValueError('drain-to-source reflection requires a half-device drain and no existing source')
    b=copy.deepcopy(half)
    original=half['nodes'];node_map={};b['nodes']=copy.deepcopy(original)
    for n in original:
        i=n['vela_id']
        if abs(n['y'])<=1e-8:
            node_map[i]=i;b['nodes'][i]['y']=0.0
        else:
            node_map[i]=len(b['nodes']);q=copy.deepcopy(n);q['y']=-q['y'];b['nodes'].append(q)
    region_name={r['id']:r['name'] for r in half['regions']}
    joined=set()
    for t in half['triangles']:
        if sum(abs(original[i]['y'])<=1e-8 for i in t['node_ids'])>=2:joined.add(t['region_id'])
    names={r['id']:r['name'] if r['id'] in joined else r['name']+'_mirrored' for r in half['regions']}
    b['regions']=_groups(list(region_name.values())+list(names.values()),'region',1000)
    regions={r['name']:r for r in b['regions']}
    contact_map={c['name']:operation['rename'].get(c['name'],c['name']) for c in half['contacts']}
    b['contacts']=_groups([c['name'] for c in half['contacts']]+list(contact_map.values()),'contact',3000)
    contact_tags={c['name']:c['physical_tag'] for c in b['contacts']}
    rename_region={region_name[i]:name for i,name in names.items()}
    interface_map={c['name']:'+'.join(rename_region[n] for n in c['name'].split('+')) for c in half['interfaces']}
    b['interfaces']=_groups([c['name'] for c in half['interfaces']]+list(interface_map.values()),'interface',4000)
    interface_tags={c['name']:c['physical_tag'] for c in b['interfaces']}
    b['triangles']=[]
    material={m['id']:m for m in b['materials']}
    for mirrored in (False,True):
        for t in half['triangles']:
            q=copy.deepcopy(t);name=names[t['region_id']] if mirrored else region_name[t['region_id']]
            q['region_id']=regions[name]['id']
            if mirrored:q['node_ids']=[node_map[t['node_ids'][i]] for i in (0,2,1)]
            q['physical_surface_tags']=sorted([regions[name]['physical_tag'],material[q['material_id']]['physical_tag']])
            b['triangles'].append(q)
    b['lines']=[];seen=set()
    old_names={c['physical_tag']:c['name'] for k in ('contacts','interfaces') for c in half[k]}
    for mirrored in (False,True):
        for line in half['lines']:
            q=copy.deepcopy(line);name=old_names[q['physical_curve_tags'][0]]
            if mirrored:
                q['node_ids']=[node_map[i] for i in reversed(line['node_ids'])]
                name=contact_map[name] if line['kind']=='contact' else interface_map[name]
            tag=contact_tags[name] if line['kind']=='contact' else interface_tags[name]
            q['physical_curve_tags']=[tag]
            key=(tag,tuple(sorted(q['node_ids'])))
            if key not in seen:b['lines'].append(q);seen.add(key)
    for field in list(b['fields'].values())+list(b['fields']['species'].values()):
        if 'values' not in field:continue
        source=list(field['values'])
        field['values'].extend(source[n['vela_id']] for n in original if node_map[n['vela_id']]!=n['vela_id'])
    owners=copy.deepcopy(half['provenance']['contact_owners'])
    for c in half['contacts']:
        name=contact_map[c['name']];mirrored=[rename_region[n] for n in owners[c['name']]]
        owners[name]=sorted(set(owners.get(name,[])+mirrored))
    b['provenance']['contact_owners']=owners
    b['provenance']['mirror']={'axis':'y','coordinate_um':0.0,'half_nodes':len(original),'half_triangles':len(half['triangles']),
        'merged_axis_nodes':sum(i==j for i,j in node_map.items()),'half_bundle_sha256':hashlib.sha256(json.dumps(half,sort_keys=True).encode()).hexdigest(),
        'region_renames':rename_region,'contact_renames':contact_map,'source_location':operation['location']}
    b['provenance']['half_geometry_audit']=b['provenance'].pop('geometry_audit',None)
    b['provenance']['mirror']['refinement_scope']='authored_half_before_reflection'
    full_areas={r['name']:0.0 for r in b['regions']}
    region_names={r['id']:r['name'] for r in b['regions']}
    for t in b['triangles']:
        a,c,d=[b['nodes'][i] for i in t['node_ids']]
        full_areas[region_names[t['region_id']]]+=((c['x']-a['x'])*(d['y']-a['y'])-(c['y']-a['y'])*(d['x']-a['x']))/2
    b['provenance']['geometry_audit']={'algorithm':'shared-axis mesh reflection','region_areas_um2':full_areas}
    b['id_mapping']={}
    for key,kind,offset in [('nodes','node',0),('triangles','triangle',0),('lines','line',len(b['triangles']))]:
        for i,item in enumerate(b[key]):item.update(vela_id=i,dfise_id=i,gmsh_tag=i+1+offset)
        b['id_mapping'][f'gmsh_{kind}_to_vela']={str(item['gmsh_tag']):item['vela_id'] for item in b[key]}
        b['id_mapping'][f'vela_{kind}_to_dfise']={str(item['vela_id']):item['dfise_id'] for item in b[key]}
    validate_bundle(b)
    return b


def load_discrete(gmsh,bundle):
    """Rebuild a discrete Gmsh export model from the authoritative full bundle."""
    gmsh.clear();gmsh.model.add('vela_reflected_mesh')
    for region in bundle['regions']:gmsh.model.addDiscreteEntity(2,region['id']+1)
    # All nodes are retained in a physical entity; cells/curves share their IDs.
    gmsh.model.mesh.addNodes(2,bundle['regions'][0]['id']+1,[n['gmsh_tag'] for n in bundle['nodes']],
                            [v for n in bundle['nodes'] for v in (n['x'],n['y'],n['z'])])
    for region in bundle['regions']:
        cells=[t for t in bundle['triangles'] if t['region_id']==region['id']]
        gmsh.model.mesh.addElementsByType(region['id']+1,2,[t['gmsh_tag'] for t in cells],[i+1 for t in cells for i in t['node_ids']])
        gmsh.model.addPhysicalGroup(2,[region['id']+1],region['physical_tag'])
        gmsh.model.setPhysicalName(2,region['physical_tag'],region['physical_name'])
    for material in bundle['materials']:
        entities=sorted({t['region_id']+1 for t in bundle['triangles'] if t['material_id']==material['id']})
        gmsh.model.addPhysicalGroup(2,entities,material['physical_tag']);gmsh.model.setPhysicalName(2,material['physical_tag'],material['physical_name'])
    for curve,group in enumerate(bundle['contacts']+bundle['interfaces'],1):
        gmsh.model.addDiscreteEntity(1,curve)
        lines=[line for line in bundle['lines'] if group['physical_tag'] in line['physical_curve_tags']]
        gmsh.model.mesh.addElementsByType(curve,1,[line['gmsh_tag'] for line in lines],[i+1 for line in lines for i in line['node_ids']])
        gmsh.model.addPhysicalGroup(1,[curve],group['physical_tag']);gmsh.model.setPhysicalName(1,group['physical_tag'],group['physical_name'])
