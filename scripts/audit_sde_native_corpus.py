"""Compare the sealed native SDE corpus with independently exported Vela meshes."""
from __future__ import annotations
import argparse
from collections import defaultdict,Counter
import hashlib,json,math,re
from pathlib import Path

from compare_sentaurus_tdr_tdx import parse_grd
from vela_sdevice.mesh_import import _edge_nodes
from sentaurus_sde_cases import load_manifest,DEFAULT_MANIFEST
from sentaurus_bjt_contract import concentration,check as check_bjt

REPO=Path(__file__).resolve().parents[1]
BJT_PROFILES=json.loads((REPO/'reference_tcad/genius_bjt_sentaurus2022/contracts/geometry_contract.json').read_text())['profiles']


def read_native(path):
    grid=parse_grd(path.with_suffix('.grd'));vertices=grid['vertices'];regions={};cells=[];contacts={}
    for name,record in grid['region_blocks'].items():
        ids=grid['region_elements'][name];material=record['material']
        if material in {'Contact','Interface'}:
            if material=='Contact':
                contacts[name]=[list(map(int,grid['elements'][i][1:])) if len(grid['elements'][i])==3 else
                                list(_edge_nodes(grid['elements'][i][1],grid['edges'])) for i in ids]
            continue
        regions[name]={'material':material,'vertices':set(),'cells':[]}
        for index in ids:
            row=grid['elements'][index]
            if row[0] not in {2,3}:raise ValueError('native volume element is not Tri3/Quad4')
            nodes=sorted({n for edge in row[1:1+(3 if row[0]==2 else 4)] for n in _edge_nodes(edge,grid['edges'])})
            cx=sum(vertices[i][0] for i in nodes)/len(nodes);cy=sum(vertices[i][1] for i in nodes)/len(nodes)
            nodes.sort(key=lambda i:math.atan2(vertices[i][1]-cy,vertices[i][0]-cx))
            regions[name]['vertices'].update(nodes);regions[name]['cells'].append(nodes);cells.append((name,nodes))
    fields=defaultdict(dict)
    pattern=r'Dataset\s*\("([^"]+)"\)\s*\{(.*?)Values\s*\((\d+)\)\s*\{(.*?)\n\s*\}'
    for match in re.finditer(pattern,path.with_suffix('.dat').read_text(),re.S):
        validity=re.search(r'validity\s*=\s*\[([^]]*)',match[2])[1]
        names=re.findall(r'"([^"]+)"',validity);values=[float(v) for v in match[4].split()]
        if len(names)!=1:raise ValueError('native auditor requires explicit region-local validity')
        region=names[0];nodes=sorted(regions[region]['vertices'])
        if len(values)!=len(nodes):raise ValueError('native region-local dataset length mismatch')
        fields[region][match[1]]=dict(zip(nodes,values))
    return {'vertices':vertices,'regions':regions,'cells':cells,'contacts':contacts,'fields':fields}


def own_model(path):
    bundle=json.loads(path.read_text());vertices=[(n['x'],n['y']) for n in bundle['nodes']]
    names={r['id']:r['name'] for r in bundle['regions']};regions={n:{'vertices':set(),'cells':[]} for n in names.values()}
    cells=[]
    for item in bundle['triangles']:
        name=names[item['region_id']];ids=item['node_ids'];cells.append((name,ids))
        regions[name]['vertices'].update(ids);regions[name]['cells'].append(ids)
    contacts={c['name']:[line['node_ids'] for line in bundle['lines'] if c['physical_tag'] in line['physical_curve_tags']]
              for c in bundle['contacts']}
    return {'vertices':vertices,'regions':regions,'cells':cells,'contacts':contacts,'bundle':bundle}


def area(points):
    return abs(sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1])))/2


def species(model,region,node):
    if 'bundle' in model:
        fields=model['bundle']['fields']
        return tuple(fields[n]['values'][node] for n in ('DonorConcentration','AcceptorConcentration','DopingConcentration'))
    fields=model['fields'][region]
    def value(aggregate,names):
        if aggregate in fields:return fields[aggregate][node]
        return sum(fields[n][node] for n in names if n in fields)
    donor=value('DonorConcentration',('ArsenicActiveConcentration','PhosphorusActiveConcentration','AntimonyActiveConcentration'))
    acceptor=value('AcceptorConcentration',('BoronActiveConcentration','IndiumActiveConcentration'))
    return donor,acceptor,fields.get('DopingConcentration',{}).get(node,donor-acceptor)


def expected(case,region,x,y):
    contract=case['contract'];kind=contract['doping']['kind'];material=contract['regions'][region]
    if material not in {'Si','PolySilicon'}:return 0.,0.
    if kind=='pn':
        rule=contract['doping'];c=rule['concentration_cm3'];j=rule['junction_x_um']
        return (c if x>=j-1e-12 else 0.),(c if x<=j+1e-12 else 0.)
    if kind=='donor':return (contract['doping']['concentration_cm3'] if material=='Si' else 0.),0.
    if kind=='bjt':
        return concentration(BJT_PROFILES,x,y)
    if kind=='mos':
        if region=='R.Polygate':return 1e20,0.
        if region!='R.Substrate':return 0.,0.
        donor=0.
        for peak,depth,factor,start in [(6e20,.05,.4,.085),(2e20,.0125,.25,.025)]:
            donor+=peak*math.exp(-math.log(peak/1.5e18)*((x/depth)**2+(max(start-abs(y),0)/(factor*depth))**2))
        return donor,1.5e18
    if kind=='bvmethods':
        if region!='R.Substrate':return 0.,0.
        donor=0.
        for peak,depth,factor,start in [(5e20,.1,.4,.165),(1e19,.05,.8,.065)]:
            sigma=depth/math.sqrt(2*math.log(peak/1e18));scale=factor*sigma
            donor+=peak*math.exp(-.5*(y/sigma)**2)*.5*(math.erf((abs(x)-start)/scale)-math.erf((abs(x)-1.13)/scale))
        return donor,1e18
    raise ValueError('unknown frozen case doping kind')


def edges(model):
    result=defaultdict(set)
    for region,ids in model['cells']:
        for a,b in zip(ids,ids[1:]+ids[:1]):result[tuple(sorted((a,b)))].add(region)
    return result


def curve_error(model,kind):
    maximum=0.;count=0
    if kind not in {'mos','bvmethods'}:return maximum,count
    counts=Counter(tuple(sorted((a,b))) for _,ids in model['cells'] for a,b in zip(ids,ids[1:]+ids[:1]))
    for (a,b),owners in edges(model).items():
        if counts[a,b]!=1:continue
        # Only exterior edges: a single volume owner can also contain two cells.
        # The caller's curve square is on the known exterior fillet.
        p,q=model['vertices'][a],model['vertices'][b]
        if kind=='mos':
            if not all(-.1012-1e-10<=v[0]<=-.0512+1e-10 and .035-1e-10<=abs(v[1])<=.085+1e-10 for v in (p,q)):continue
            if p[1]*q[1]<0:continue
            transform=lambda v:(v[0],abs(v[1]));center=(-.0512,.035);radius=.05
        else:
            if not all(.085-1e-10<=abs(v[0])<=.165+1e-10 and -.202-1e-10<=v[1]<=-.122+1e-10 for v in (p,q)):continue
            if p[0]*q[0]<0:continue
            transform=lambda v:(abs(v[0]),v[1]);center=(.085,-.122);radius=.08
        p,q=transform(p),transform(q)
        for v in (p,q,((p[0]+q[0])/2,(p[1]+q[1])/2)):
            maximum=max(maximum,abs(radius-math.hypot(v[0]-center[0],v[1]-center[1])))
        count+=1
    return maximum,count


def audit_case(case,native,own):
    contract=case['contract'];kind=contract['doping']['kind'];failures=[]
    if set(native['regions'])!=set(contract['regions']):failures.append('native volume region identity')
    aliases=own['bundle']['provenance']['material_aliases']
    for name,record in native['regions'].items():
        if aliases.get(record['material'],record['material'])!=contract['regions'].get(name):
            failures.append('native material identity '+name)
    points=native['vertices']
    bounds=[min(p[0] for p in points),min(p[1] for p in points),max(p[0] for p in points),max(p[1] for p in points)]
    if max(abs(a-b) for a,b in zip(bounds,contract['bounds_um']))>1e-8:failures.append('native bounds')
    point_errors=[0.,0.,0.];worst=[None,None,None]
    for region,data in native['regions'].items():
        if region not in contract['regions']:continue
        for node in data['vertices']:
            x,y=points[node];nd,na=expected(case,region,x,y)
            for i,(actual,reference) in enumerate(zip(species(native,region,node),(nd,na,nd-na))):
                error=abs(actual-reference)/max(abs(reference),1.)
                if error>point_errors[i]:point_errors[i]=error;worst[i]={'xy_um':[x,y],'expected':reference,'actual':actual}
    if max(point_errors)>1e-6:failures.append('native pointwise doping')
    doses={};areas={};curves={}
    for label,model in [('native',native),('own',own)]:
        dose=[0.,0.];region_area=defaultdict(float)
        for region,ids in model['cells']:
            volume=area([model['vertices'][i] for i in ids]);region_area[region]+=volume
            if contract['regions'].get(region) not in {'Si','PolySilicon'}:continue
            for k in (0,1):dose[k]+=volume*sum(species(model,region,n)[k] for n in ids)/len(ids)
        doses[label]=dose;areas[label]=dict(region_area);curves[label]=curve_error(model,kind)
    integral=[abs(a-b)/max(abs(b),1.) for a,b in zip(doses['own'],doses['native'])]
    if max(integral)>.01:failures.append('comparative doping integral')
    for region,reference in contract['region_area_um2'].items():
        budget=1e-8
        if region.startswith('R.Spacer') and kind in {'mos','bvmethods'}:
            radius=.05 if kind=='mos' else .08
            budget+=math.pi*radius/2*curves['native'][0]
        if abs(areas['native'].get(region,0)-reference)>budget:failures.append('native region area '+region)
    contact_errors={};contact_geometry={}
    if set(native['contacts'])!=set(contract['contacts']):failures.append('native contact identity')
    for name in contract['contacts']:
        lengths=[]
        for model in (own,native):
            lengths.append(sum(math.dist(model['vertices'][a],model['vertices'][b]) for a,b in model['contacts'].get(name,[])))
        contact_errors[name]=abs(lengths[0]-lengths[1])
        if contact_errors[name]>1e-8:failures.append('native contact length '+name)
        rule=contract['contacts'][name]
        nodes={n for edge in native['contacts'].get(name,[]) for n in edge}
        if 'segments_um' in rule:
            def distance(point,segment):
                a,b=segment;dx,dy=b[0]-a[0],b[1]-a[1]
                t=max(0.,min(1.,((point[0]-a[0])*dx+(point[1]-a[1])*dy)/(dx*dx+dy*dy)))
                return math.hypot(point[0]-a[0]-t*dx,point[1]-a[1]-t*dy)
            samples=[points[n] for n in nodes]+[tuple((points[a][k]+points[b][k])/2 for k in (0,1)) for a,b in native['contacts'].get(name,[])]
            normal=max((min(distance(q,seg) for seg in rule['segments_um']) for q in samples),default=math.inf)
            span_error=max(min((math.dist(endpoint,points[n]) for n in nodes),default=math.inf)
                           for seg in rule['segments_um'] for endpoint in seg)
        else:
            axis='xy'.index(rule['axis'])
            normal=max((abs(points[n][axis]-rule['coordinate_um']) for n in nodes),default=math.inf)
            span=sorted(points[n][1-axis] for n in nodes)
            span_error=max(abs(a-b) for a,b in zip((span[0],span[-1]),rule['span_um'])) if span and 'span_um' in rule else 0.
        contact_geometry[name]={'normal_error_um':normal,'span_error_um':span_error}
        if normal>1e-8 or span_error>1e-8:failures.append('native contact position '+name)
    interface={}
    for label,model in [('native',native),('own',own)]:
        values=defaultdict(float)
        for (a,b),owners in edges(model).items():
            if len(owners)==2:values['|'.join(sorted(owners))]+=math.dist(model['vertices'][a],model['vertices'][b])
        interface[label]=dict(values)
    if set(interface['native'])!=set(interface['own']):failures.append('native interface identity')
    for name,length in interface['own'].items():
        if abs(length-interface['native'].get(name,0))>1e-8:failures.append('native interface length '+name)
    from sentaurus_device_ir import parse_sde_device_ir
    from sentaurus_ir_v2 import migrate_device_ir_v1
    from sde_native_metrics import cuts,sizing
    ir=migrate_device_ir_v1(parse_sde_device_ir(REPO/case['source'],parameters=case.get('parameters',{})))
    cut_report=cuts(case,ir,native,own,species)
    if any(not c['passed'] for c in cut_report):failures.append('comparative cut junction')
    local_size={label:sizing(ir,model,edges(model)) for label,model in [('native',native),('own',own)]}
    for label,rows in local_size.items():
        for row in rows:
            if not row['passed']:failures.append(label+' local spacing '+row['name'])
    return {'id':case['id'],'status':'failed' if failures else 'passed','failures':failures,
            'cuts':cut_report,'local_spacing':local_size,
            'native_pointwise_relative_error':point_errors,'worst_native_points':worst,'doses_cm3_um2':doses,
            'integral_relative_error':integral,'region_areas_um2':areas,'curve_error_um':curves,
            'native_bounds_um':bounds,'contact_position_error_um':contact_geometry,'contact_length_error_um':contact_errors,'interface_lengths_um':interface}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-root',type=Path,required=True);parser.add_argument('--own-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--run-index',default='1');args=parser.parse_args();rows=[]
    for case in load_manifest(DEFAULT_MANIFEST)['cases']:
        try:rows.append(audit_case(case,read_native(args.native_root/case['id']/'native'),own_model(args.own_root/case['id']/args.run_index/'mesh_field_bundle.json')))
        except Exception as error:rows.append({'id':case['id'],'status':'failed','error':str(error)})
        print(case['id'],rows[-1]['status'],rows[-1].get('failures',rows[-1].get('error')),flush=True)
    report={'schema':'vela.sde_native_corpus_audit.v1','status':'passed' if all(r['status']=='passed' for r in rows) else 'failed','run_index':args.run_index,'cases':rows}
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return 0 if report['status']=='passed' else 2


if __name__=='__main__':raise SystemExit(main())
