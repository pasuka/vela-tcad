"""Compare local temperature-dependent silicon physics with native node fields."""
import argparse,csv,hashlib,json,math,os,subprocess
from pathlib import Path
from audit_templates_ldmos_g3_idvg_shift_kcl import scalar_field

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('export','probe','output'):parser.add_argument('--'+name,required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    names=('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential','LatticeTemperature','DonorConcentration','AcceptorConcentration','eDensity','hDensity','BandGap','BandgapNarrowing','ElectronAffinity','ConductionBandEnergy','ValenceBandEnergy','srhRecombination')
    fields={name:scalar_field(args.export,name) for name in names};nodes=sorted(fields['eDensity'])
    assert all(set(values)==set(nodes) for values in fields.values())
    states=[dict(id=n,potential_V=fields['ElectrostaticPotential'][n],electron_qf_V=fields['eQuasiFermiPotential'][n],hole_qf_V=fields['hQuasiFermiPotential'][n],temperature_K=fields['LatticeTemperature'][n],donors_m3=fields['DonorConcentration'][n]*1e6,acceptors_m3=fields['AcceptorConcentration'][n]*1e6) for n in nodes]
    source=args.output/'input.json';source.write_text(json.dumps(dict(states_SI=states)),encoding='utf-8')
    env=dict(os.environ)
    if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;'+env['PATH']
    output=args.output/'output.json'
    with output.open('w',encoding='utf-8') as out:subprocess.run([str(args.probe.resolve()),str(source.resolve())],stdout=out,env=env,check=True)
    data=json.loads(output.read_text(encoding='utf-8'));rows=data['results'];assert [r['id'] for r in rows]==nodes
    mapping={'bandgap_eV':('BandGap',1.),'affinity_eV':('ElectronAffinity',1.),'conduction_band_eV':('ConductionBandEnergy',1.),'valence_band_eV':('ValenceBandEnergy',1.),'electrons_m3':('eDensity',1e6),'holes_m3':('hDensity',1e6),'srh_m3_per_s':('srhRecombination',1e6)}
    summary={};records=[]
    for key,(field,scale) in mapping.items():
        errors=[];relative=[]
        reference=[fields[field][n]*scale for n in nodes];floor=max(abs(v) for v in reference)*1e-12
        for row,ref in zip(rows,reference):
            error=row[key]['value']-ref;errors.append(error)
            if abs(ref)>max(1e-100,floor):relative.append(abs(error/ref))
            records.append(dict(node_id=row['id'],quantity=key,native=ref,vela=row[key]['value'],difference=error))
        relative.sort();summary[key]=dict(max_abs=max(map(abs,errors)),rms_abs=math.sqrt(sum(e*e for e in errors)/len(errors)),resolved_nodes=len(relative),relative_floor=floor,relative_median=relative[len(relative)//2] if relative else None,relative_max=max(relative) if relative else None)
    with (args.output/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
    reference_potential=[fields['ConductionBandEnergy'][n]+fields['ElectronAffinity'][n]+fields['ElectrostaticPotential'][n] for n in nodes]
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    result=dict(scope='Frozen native local temperature state; no coupled solve or electrical gates',nodes=len(nodes),reference_potential_V=dict(vela=data['reference_potential_V'],native_min=min(reference_potential),native_max=max(reference_potential)),results=summary,sha256={str(p):sha(p) for p in [args.probe,Path(__file__),args.export/'field_manifest.json',*list((args.export/'fields').glob('*.csv'))]})
    (args.output/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result['results'],indent=2))

if __name__=='__main__':main()
