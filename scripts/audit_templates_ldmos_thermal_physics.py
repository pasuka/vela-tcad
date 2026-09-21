"""Compare local temperature-dependent silicon physics with native node fields."""
import argparse,csv,hashlib,json,math,os,re,subprocess
from pathlib import Path
from audit_templates_ldmos_g3_idvg_shift_kcl import scalar_field

def auger_parameters(path):
    text=path.read_text(encoding='utf-8')
    match=re.search(r'^Auger\b[^\n]*\n\{(.*?)\n\}',text,re.M|re.S)
    if not match:raise ValueError('Missing audited Auger parameter block')
    result={}
    for key in ('A','B','C','H','N0'):
        row=re.search(r'^\s*'+key+r'\s*=\s*([^#\n]+)',match[1],re.M)
        if not row:raise ValueError('Missing Auger '+key)
        pair=[float(v.strip()) for v in row[1].split(',')]
        if len(pair)!=2 or not all(math.isfinite(v) for v in pair):raise ValueError('Invalid Auger pair')
        factor=1e-12 if key in ('A','B','C') else 1e6 if key=='N0' else 1.
        result[key]=[v*factor for v in pair]
    if min(result['N0'])<=0 or min(result['H'])<0:raise ValueError('Invalid Auger scales')
    return result

def auger_at_native_density(n,p,t,fn,fp,parameters,with_generation=False):
    """Isolate the documented Auger law using native carriers and QF splitting.

    SI densities/rates; the Fermi equilibrium product is obtained from the
    native QF splitting. Negative splitting uses logarithms to avoid overflow.
    """
    if min(n,p,t)<=0:raise ValueError('Positive native densities and temperature required')
    split=(fp-fn)/(1.380649e-23/1.602176634e-19*t)
    excess=n*p*(-math.expm1(-split)) if split>=0 else math.exp(math.log(n)+math.log(p)-split)*math.expm1(split)
    ratio=t/300.
    coefficients=[(parameters['A'][i]+ratio*(parameters['B'][i]+ratio*parameters['C'][i]))*
                  (1+parameters['H'][i]*math.exp(-density/parameters['N0'][i])) for i,density in enumerate((n,p))]
    rate=(coefficients[0]*n+coefficients[1]*p)*excess
    return rate if with_generation else max(0.,rate)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('export','probe','output'):parser.add_argument('--'+name,required=True,type=Path)
    parser.add_argument('--parameters',type=Path,help='Audited Siliconc100.par for native-density Auger isolation')
    parser.add_argument('--require-auger',action='store_true')
    parser.add_argument('--vela-state',type=Path,help='Compare a reclosed Vela state instead of evaluating the native state')
    parser.add_argument('--auger-with-generation',action='store_true',help='Explicit signed-law diagnostic; original D0 leaves this off')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    names=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','LatticeTemperature','DonorConcentration','AcceptorConcentration','eDensity','hDensity','BandGap','BandgapNarrowing','ElectronAffinity','ConductionBandEnergy','ValenceBandEnergy','srhRecombination')
    manifest=json.loads((args.export/'field_manifest.json').read_text(encoding='utf-8'))
    auger_names={f['name'] for f in manifest['fields'] if f['name'].casefold()=='augerrecombination'}
    if len(auger_names)>1 or (args.require_auger and not auger_names):raise ValueError('Missing/ambiguous native Auger field')
    if auger_names:names=(*names,next(iter(auger_names)))
    fields={name:scalar_field(args.export,name) for name in names};nodes=sorted(fields['eDensity'])
    assert all(set(values)==set(nodes) for values in fields.values())
    states=[dict(id=n,potential_V=fields['ElectrostaticPotential'][n],electron_qf_V=fields['eQuasiFermiPotential'][n],hole_qf_V=fields['hQuasiFermiPotential'][n],temperature_K=fields['LatticeTemperature'][n],donors_m3=fields['DonorConcentration'][n]*1e6,acceptors_m3=fields['AcceptorConcentration'][n]*1e6) for n in nodes]
    if args.vela_state:
        solved=json.loads(args.vela_state.read_text(encoding='utf-8'));origin=solved.get('potential_origin_V',0.)
        x=solved['referenced_state_interleaved']
        for state in states:
            n=state['id'];state.update(potential_V=x[4*n]+origin,electron_qf_V=x[4*n+1],hole_qf_V=x[4*n+2],temperature_K=x[4*n+3],
                electron_qf_reference_V=solved['electron_qf_reference_V'][n]+origin,hole_qf_reference_V=solved['hole_qf_reference_V'][n]+origin)
    else:
        for state in states:
            n=state['id'];state.update(reference_conduction_band_eV=fields['ConductionBandEnergy'][n],reference_valence_band_eV=fields['ValenceBandEnergy'][n])
    source=args.output/'input.json';source.write_text(json.dumps(dict(states_SI=states,auger_with_generation=args.auger_with_generation)),encoding='utf-8')
    env=dict(os.environ)
    if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;'+env['PATH']
    output=args.output/'output.json'
    with output.open('w',encoding='utf-8') as out:subprocess.run([str(args.probe.resolve()),str(source.resolve())],stdout=out,env=env,check=True)
    data=json.loads(output.read_text(encoding='utf-8'));rows=data['results'];assert [r['id'] for r in rows]==nodes
    mapping={'bandgap_eV':('BandGap',1.),'affinity_eV':('ElectronAffinity',1.),'conduction_band_eV':('ConductionBandEnergy',1.),'valence_band_eV':('ValenceBandEnergy',1.),'electrons_m3':('eDensity',1e6),'holes_m3':('hDensity',1e6),'srh_m3_per_s':('srhRecombination',1e6)}
    mapping['bgn_eV']=('BandgapNarrowing',1.)
    if not args.vela_state:
        mapping.update(electrons_at_native_band_m3=('eDensity',1e6),holes_at_native_band_m3=('hDensity',1e6))
    if auger_names:mapping['auger_m3_per_s']=(next(iter(auger_names)),1e6)
    if args.parameters:
        if not auger_names:raise ValueError('Auger isolation requires a native Auger field')
        parameters=auger_parameters(args.parameters)
        for row,n in zip(rows,nodes):
            row['auger_native_density_m3_per_s']=auger_at_native_density(fields['eDensity'][n]*1e6,fields['hDensity'][n]*1e6,
                fields['LatticeTemperature'][n],fields['eQuasiFermiPotential'][n],fields['hQuasiFermiPotential'][n],parameters,args.auger_with_generation)
        mapping['auger_native_density_m3_per_s']=(next(iter(auger_names)),1e6)
    summary={};records=[]
    for key,(field,scale) in mapping.items():
        errors=[];relative=[];roundoff_excluded=0
        reference=[fields[field][n]*scale for n in nodes];floor=max(abs(v) for v in reference)*1e-12
        for row,ref in zip(rows,reference):
            value=row[key]['value'] if isinstance(row[key],dict) else row[key]
            error=value-ref;errors.append(error)
            n=row['id'];qf_scale=max(1.,abs(fields['ElectrostaticPotential'][n]),abs(fields['eQuasiFermiPotential'][n]),abs(fields['hQuasiFermiPotential'][n]))
            qf_resolved=abs(fields['hQuasiFermiPotential'][n]-fields['eQuasiFermiPotential'][n])>32*math.ulp(qf_scale)
            rate=key.startswith(('srh_','auger_'))
            if rate and not qf_resolved:roundoff_excluded+=1
            elif abs(ref)>max(1e-100,floor):relative.append(abs(error/ref))
            records.append(dict(node_id=row['id'],quantity=key,native=ref,vela=value,difference=error))
        relative.sort();summary[key]=dict(max_abs=max(map(abs,errors)),rms_abs=math.sqrt(sum(e*e for e in errors)/len(errors)),resolved_nodes=len(relative),relative_floor=floor,relative_median=relative[len(relative)//2] if relative else None,relative_max=max(relative) if relative else None,
            relative_qf_roundoff_excluded=roundoff_excluded,relative_qf_resolution='Rates require |fp-fn| > 32 ULP of max(1 V, |psi|, |fn|, |fp|); absolute errors retain every node')
    with (args.output/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    reference_potential=[fields['ConductionBandEnergy'][n]+fields['ElectronAffinity'][n]+fields['ElectrostaticPotential'][n] for n in nodes]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sources=[args.probe,Path(__file__),source,output,args.export/'field_manifest.json',*list((args.export/'fields').glob('*.csv'))]
    sources += [p for p in (args.parameters,args.vela_state) if p is not None]
    result=dict(scope='Reclosed Vela local state versus native fields; diagnostic field errors, not new acceptance gates' if args.vela_state else 'Frozen native local temperature state; no coupled solve or electrical gates',nodes=len(nodes),reference_potential_V=dict(vela=data['reference_potential_V'],native_min=min(reference_potential),native_max=max(reference_potential)),results=summary,sha256={str(p):sha(p) for p in sources})
    result['auger_with_generation']=args.auger_with_generation
    (args.output/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result['results'],indent=2))

if __name__=='__main__':main()
