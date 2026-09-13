"""Prepare explicit SI geometry and frozen native state for the experimental D0 operator.
This is a state audit, not a self-consistent curve or acceptance runner.
"""
import argparse, csv, hashlib, json, math
from pathlib import Path


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def fields(export, name, count, required=True, required_nodes=()):
    values={}
    for path in sorted((export/'fields').glob(name+'_region*.csv')):
        for row in csv.DictReader(path.open(encoding='utf-8')):
            i=int(row['node_id']);v=float(row['component0'])
            if not 0<=i<count or not math.isfinite(v):raise ValueError(f'Invalid field {name} node/value')
            if i in values and not math.isclose(v,values[i],rel_tol=1e-9,abs_tol=1e-9):
                raise ValueError(f'Inconsistent shared node {i} for {name}')
            values[i]=v
    if not set(required_nodes)<=set(values):raise ValueError(f'Incomplete semiconductor field {name}')
    if required and len(values)!=count:
        raise ValueError(f'Incomplete field {name}: {len(values)}/{count}')
    return [values.get(i,0.) for i in range(count)]


def contact_boundary_lengths(mesh, contact, coordinate_to_metres):
    """Lumped contact measure; reject missing, duplicate or interior edges."""
    if not math.isfinite(coordinate_to_metres) or coordinate_to_metres<=0.:
        raise ValueError('Positive coordinate conversion required')
    nodes={n['id']:(n['x'],n['y']) for n in mesh['nodes']}
    incidence={}
    for cell in mesh['triangles']:
        ns=cell['node_ids']
        for i,j in zip(ns,ns[1:]+ns[:1]):
            incidence.setdefault(tuple(sorted((i,j))),[]).append(cell['region_id'])
    lengths={i:0. for i in contact['node_ids']};seen=set()
    for a,b in contact.get('edge_node_ids',[]):
        key=tuple(sorted((a,b)))
        if a not in lengths or b not in lengths or key in seen or incidence.get(key)!=[0]:
            raise ValueError('Finite contact requires unique exterior silicon edges')
        seen.add(key);half=.5*math.dist(nodes[a],nodes[b])*coordinate_to_metres
        lengths[a]+=half;lengths[b]+=half
    if not lengths or any(v<=0. for v in lengths.values()):
        raise ValueError('Finite contact node has no positive boundary measure')
    return lengths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('thermal-input','export','profile','couples','output'):
        p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--gate',type=float,required=True);p.add_argument('--drain',type=float,required=True)
    p.add_argument('--state-result',type=Path)
    p.add_argument('--isothermal-state',type=Path)
    p.add_argument('--freeze-temperature',type=float)
    p.add_argument('--initial-temperature',type=float)
    p.add_argument('--potential-origin',type=float,default=0.)
    p.add_argument('--diagnostic-newton-iterations',type=int,default=0)
    p.add_argument('--hole-recombination-velocity-cm-s',type=float,
                   help='Explicit constant source/drain hRecVelocity; original D0 is 1.93e6 cm/s')
    p.add_argument('--auger-with-generation',action='store_true',help='Explicit signed Auger control; not enabled by original D0')
    p.add_argument('--native-poisson-debug',type=Path,help='Explicit audited AverageBox edge/charge-volume pair')
    p.add_argument('--native-poisson-export',type=Path,help='Native nodes.csv/elements.csv for strict geometry identity')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    thermal=read(a.thermal_input);mesh=read(Path(thermal['mesh_file']));profile=read(a.profile)
    nodes=mesh['nodes'];count=len(nodes);assert [n['id'] for n in nodes]==list(range(count))
    xy=[(n['x'],n['y']) for n in nodes];length=thermal['coordinate_to_metres']
    eps0=8.8541878128e-12;area=[0.]*count;global_area=[0.]*count;edge={};edge_eps={}
    for cell in mesh['triangles']:
        ns=cell['node_ids'];assert len(ns)==3
        x0,y0=xy[ns[0]];x1,y1=xy[ns[1]];x2,y2=xy[ns[2]]
        double_area=abs((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))
        if double_area<=0.:raise ValueError('Degenerate triangle')
        silicon=cell['region_id']==0;eps=eps0*(11.7 if silicon else 3.9)
        for i in ns:global_area[i]+=double_area*length**2/6.
        if silicon:
            for i in ns:area[i]+=double_area*length**2/6.
        for j in range(3):
            i,k,opposite=ns[j],ns[(j+1)%3],ns[(j+2)%3]
            u=(xy[i][0]-xy[opposite][0],xy[i][1]-xy[opposite][1])
            v=(xy[k][0]-xy[opposite][0],xy[k][1]-xy[opposite][1])
            weight=(u[0]*v[0]+u[1]*v[1])/(2.*double_area)
            # Match the qualified BoxGeometryBuilder negative-cotangent fallback.
            if weight<0.:weight=double_area/(6.*math.dist(xy[i],xy[k])**2)
            key=tuple(sorted((i,k)));edge[key]=edge.get(key,0.)+weight
            edge_eps.setdefault(key,[]).append(eps)
    couples={tuple(sorted((int(r['node0']),int(r['node1'])))):float(r['couple_m']) for r in csv.DictReader(a.couples.open())}
    if not set(couples)<=set(edge):raise ValueError('Unknown AverageBox edge')
    geometry=[]
    for (i,j),weight in sorted(edge.items()):
        eps=weight*sum(edge_eps[(i,j)])/len(edge_eps[(i,j)])
        distance=math.dist(xy[i],xy[j])*length
        geometry.append(dict(nodes=[i,j],poisson_F_per_m=eps,transport_weight=couples.get((i,j),0.)/distance))
    names=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','LatticeTemperature')
    data=[fields(a.export,name,count) for name in names]
    state=[(0. if area[i]==0. and k in (1,2) else data[k][i]) for i in range(count) for k in range(4)]
    donors=[v*1e6 for v in fields(a.export,'DonorConcentration',count,False,[i for i in range(count) if area[i]>0.])]
    acceptors=[v*1e6 for v in fields(a.export,'AcceptorConcentration',count,False,[i for i in range(count) if area[i]>0.])]
    mobility=profile['solver']['mobility']
    # The supplied linked profile uses TCAD cm/s; this new operator is SI.
    for name in ('electron_saturation_velocity_m_s','hole_saturation_velocity_m_s'):
        mobility[name]*=.01
    mobility['ialmob']['geometry_file']=str(Path(mobility['ialmob']['geometry_file'].replace('@workspace/',str(Path.cwd())+'/')).resolve())
    boundaries=[]
    for contact in mesh['contacts']:
        name=contact['name']
        if name=='th_lat':continue
        if name=='gate':kind='psi';value=a.gate-next(c['flatband_voltage'] for c in profile['contacts'] if c['name']=='gate')
        elif name in ('source','drain','substrate'):kind='neutral_contact';value=a.drain if name=='drain' else 0.
        else:raise ValueError(f'Unknown contact {name}')
        lengths=None
        if name in ('source','drain') and a.hole_recombination_velocity_cm_s is not None:
            if not math.isfinite(a.hole_recombination_velocity_cm_s) or a.hole_recombination_velocity_cm_s<0.:
                raise ValueError('Invalid hole recombination velocity')
            lengths=contact_boundary_lengths(mesh,contact,length)
        for i in contact['node_ids']:
            boundary=dict(node=i,kind=kind,value=value)
            if lengths is not None:
                boundary.update(hole_recombination_velocity_m_per_s=a.hole_recombination_velocity_cm_s*.01,boundary_length_m=lengths[i])
            boundaries.append(boundary)
    if a.freeze_temperature is not None:
        boundaries.extend(dict(node=i,kind='temperature',value=a.freeze_temperature) for i in range(count))
    cfg={k:thermal[k] for k in ('mesh_file','coordinate_to_metres','region_conductivity','thermodes')}
    retained=None
    if a.state_result:
        previous=read(a.state_result);state=previous['state_interleaved']
        if 'referenced_state_interleaved' in previous:
            retained=[[],[],[],[]]
            for i in range(count):
                for k,key in enumerate(('electron_qf_reference_V','hole_qf_reference_V')):
                    ref=previous[key][i];inc=previous['referenced_state_interleaved'][4*i+1+k];offset=previous.get('potential_origin_V',0.)
                    absolute=ref+offset;retained[k].append(absolute);retained[k+2].append(math.fsum((inc,ref,offset,-absolute)))
        for i in range(count):
            for k in range(3):
                if area[i]>0. or k==0:state[4*i+k]+=previous.get('potential_origin_V',0.)
    if a.isothermal_state:
        rows=list(csv.DictReader(a.isothermal_state.open()))
        assert [int(r['node_id']) for r in rows]==list(range(count))
        source_node=next(c for c in mesh['contacts'] if c['name']=='source')['node_ids'][0]
        offset=-float(rows[source_node]['phin'])
        if 'electron_qf_increment_V' in rows[0]:
            retained=[[],[],[],[]]
            for row in rows:
                for k,carrier in enumerate(('electron','hole')):
                    ref=float(row[carrier+'_qf_reference_V']);inc=float(row[carrier+'_qf_increment_V']);absolute=ref+offset
                    retained[k].append(absolute);retained[k+2].append(math.fsum((inc,ref,offset,-absolute)))
        for i,row in enumerate(rows):
            for k,name in enumerate(('psi','phin','phip')):state[4*i+k]=float(row[name])+offset if area[i]>0 or k==0 else 0.
    initial_temperature=a.freeze_temperature if a.freeze_temperature is not None else a.initial_temperature
    if initial_temperature is not None:
        for i in range(count):state[4*i+3]=initial_temperature
    for i in range(count):
        for k in range(3):
            if area[i]>0. or k==0:state[4*i+k]-=a.potential_origin
    for b in boundaries:
        if b['kind']!='temperature':b['value']-=a.potential_origin
    if retained:
        raw=list(state);refs=[[],[]]
        for i in range(count):
            for k in range(2):
                if area[i]>0.:
                    old_ref=retained[k][i];new_ref=old_ref-a.potential_origin
                    raw[4*i+1+k]=math.fsum((retained[k+2][i],old_ref,-a.potential_origin,-new_ref));refs[k].append(new_ref)
                else:raw[4*i+1+k]=0.;refs[k].append(0.)
        cfg.update(referenced_state_interleaved=raw,electron_qf_reference_V=refs[0],hole_qf_reference_V=refs[1])
    cfg.update(recombination_area_m2=[global_area[i] if area[i]>0. else 0. for i in range(count)],potential_origin_V=a.potential_origin,electrical_current_scale_A_per_m=1.602176634e-19*max(abs(d-n) for d,n in zip(donors,acceptors))*.1417*(1.380649e-23*300/1.602176634e-19),electrical_gate_solver={k:profile['solver'][k] for k in ('carrier_row_convergence','continuity_row_scaling','block_absolute_convergence')},diagnostic_newton_max_iterations=a.diagnostic_newton_iterations,silicon_area_m2=area,fixed_charge_C_per_m=[0.]*count,edge_geometry=geometry,donors_m3=donors,acceptors_m3=acceptors,state_interleaved=state,mobility_SI=mobility,boundaries=boundaries,auger_with_generation=a.auger_with_generation)
    if bool(a.native_poisson_debug)!=bool(a.native_poisson_export):raise ValueError('Both native Poisson debug and geometry export are required')
    if a.native_poisson_debug:
        from templates_ldmos_native_poisson import apply_native_poisson
        cfg=apply_native_poisson(cfg,mesh,a.native_poisson_debug,a.native_poisson_export)
    output=a.output/'input.json';output.write_text(json.dumps(cfg),encoding='utf-8')
    sources=[Path(__file__),a.thermal_input,a.profile,a.couples,Path(thermal['mesh_file']),Path(mobility['ialmob']['geometry_file']),*sorted((a.export/'fields').glob('*.csv'))]
    if a.state_result:sources.append(a.state_result)
    if a.isothermal_state:sources.append(a.isothermal_state)
    if a.native_poisson_debug:
        sources.extend(Path(p) for p in cfg['poisson_geometry_provenance']['sources_sha256'])
    manifest=dict(scope='Explicit experimental SI four-equation input; provenance lists the actual native or Vela seed; not a qualification report',gate_V=a.gate,drain_V=a.drain,nodes=count,cells=len(mesh['triangles']),transport_edges=sum(e['transport_weight']>0 for e in geometry),input_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in manifest.items() if k!='sources_sha256'}))

if __name__=='__main__':main()
