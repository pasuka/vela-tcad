"""Derive SLOT intrinsic references from native bands and prepare isolated controls."""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from audit_slot_ldmos_element_ialmob import validate_export_mesh
from run_slot_ldmos_noavalanche_continuation import redirect_diagnostics


def derive_reference(export: Path, thermal_voltage: float) -> dict:
    if not math.isfinite(thermal_voltage) or thermal_voltage<=0:
        raise ValueError('Native thermal voltage must be positive finite')
    metadata=json.loads((export/'metadata.json').read_text())
    units={(f['name'],f['region']):f['unit'] for f in metadata['fields']}
    expected={'eDensity':'cm^-3','hDensity':'cm^-3','ElectrostaticPotential':'V',
        'eQuasiFermiPotential':'V','hQuasiFermiPotential':'V','ConductionBandEnergy':'eV',
        'ValenceBandEnergy':'eV','ElectronAffinity':'eV','BandGap':'eV','BandgapNarrowing':'eV'}
    result={}
    for name,region in [('Si',0),('PolySilicon',3)]:
        data={}
        for field,unit in expected.items():
            if units.get((field,region))!=unit:raise ValueError('Native band unit mismatch: '+field)
            with (export/'fields'/f'{field}_region{region}.csv').open(newline='') as stream:
                rows=list(csv.DictReader(stream))
            data[field]={int(r['node_id']):float(r['component0']) for r in rows}
            if len(rows)!=len(data[field]) or not rows:raise ValueError('Invalid field support')
        ids=set(data['eDensity'])
        if any(set(values)!=ids for values in data.values()):raise ValueError('Native field support mismatch')
        nc=[];nv=[];origin=[];chi=[];ni=[]
        for i in sorted(ids):
            n,p=data['eDensity'][i],data['hDensity'][i]
            if not (n>0 and p>0):raise ValueError('Nonpositive native carrier density')
            nc.append(n*math.exp((data['ConductionBandEnergy'][i]+data['eQuasiFermiPotential'][i])/thermal_voltage))
            nv.append(p*math.exp((-data['ValenceBandEnergy'][i]-data['hQuasiFermiPotential'][i])/thermal_voltage))
            origin.append(data['ConductionBandEnergy'][i]+data['ElectrostaticPotential'][i]+data['ElectronAffinity'][i])
            chi.append(data['ElectronAffinity'][i]-.5*data['BandgapNarrowing'][i])
            ni.append(math.sqrt(n*p)*math.exp(-(data['BandgapNarrowing'][i]+data['hQuasiFermiPotential'][i]-data['eQuasiFermiPotential'][i])/(2*thermal_voltage)))
        def constant(values,label):
            mean=statistics.median(values)
            if not all(math.isfinite(v) for v in values) or max(values)-min(values)>1e-10*max(abs(mean),1.):
                raise ValueError('Native field is not material-constant: '+label)
            return mean
        values={'Nc_cm3':constant(nc,'Nc'),'Nv_cm3':constant(nv,'Nv'),'ni_cm3':constant(ni,'ni'),
            'vacuum_reference_eV':constant(origin,'origin'),'affinity_eV':constant(chi,'affinity'),
            'bandgap_eV':constant(list(data['BandGap'].values()),'Eg')}
        values['intrinsic_reference_V']=(values['vacuum_reference_eV']-values['affinity_eV']
            -.5*values['bandgap_eV']-.5*thermal_voltage*math.log(values['Nc_cm3']/values['Nv_cm3']))
        result[name]=values
    if abs(result['Si']['vacuum_reference_eV']-result['PolySilicon']['vacuum_reference_eV'])>1e-10:
        raise ValueError('Native materials do not share the same vacuum energy origin')
    # A common origin cancels; the device reference remains the calibrated Si one.
    result['PolySilicon']['relative_to_Si_V']=result['PolySilicon']['intrinsic_reference_V']-result['Si']['intrinsic_reference_V']
    return result


def prepare(base_path: Path,state: Path,export: Path,output: Path,thermal_voltage: float) -> dict:
    base_path,state,export,output=[p.resolve() for p in (base_path,state,export,output)]
    if output.exists():raise ValueError('Refusing to overwrite energy-reference controls')
    base=json.loads(base_path.read_text());solver=base['solver']
    if base.get('scaling',{}).get('mode')!='unit_scaling' or solver.get('temperature_K',300)!=300:
        raise ValueError('This SLOT control requires unit_scaling and 300 K')
    for key in ('mesh_file','node_doping_file','materials_file'):
        base[key]=str((base_path.parent/base[key]).resolve())
    mesh=json.loads(Path(base['mesh_file']).read_text());validate_export_mesh(mesh,export)
    native=derive_reference(export,thermal_voltage)
    materials=json.loads(Path(base['materials_file']).read_text())
    if solver.get('recombination_materials')!=['Si']:raise ValueError('Expected Si-only recombination')
    if not all(c['type']=='ohmic' for c in base['contacts']):raise ValueError('Expected restored ohmic contacts')
    output.mkdir(parents=True);cases={}
    for name in ('reference_only','aligned_plain','aligned_full','aligned_material_sources'):
        folder=output/name;folder.mkdir()
        cfg=copy.deepcopy(base);mats=copy.deepcopy(materials)
        mats['_comment']='Restored semiconductor PolySilicon with native doping, ohmic contacts and a band-derived intrinsic reference.'
        poly=next(m for m in mats['materials'] if m['name']=='PolySilicon')
        poly['intrinsic_reference_potential_V']=native['PolySilicon']['relative_to_Si_V']
        if name!='reference_only':
            for key in ('Si','PolySilicon'):
                mat=next(m for m in mats['materials'] if m['name']==key);source=native[key]
                mat.update(bandgap_eV=source['bandgap_eV'],electron_affinity_eV=source['affinity_eV'],
                    Nc_m3=source['Nc_cm3'],Nv_m3=source['Nv_cm3'])
        if name in ('aligned_full','aligned_material_sources'):
            cfg['solver']['auger_density_dependence']={'enabled':True,
                'electron':{'enhancement':3.46667,'reference_density_m3':1e18},
                'hole':{'enhancement':8.25688,'reference_density_m3':1e18}}
        if name=='aligned_material_sources':
            cfg['solver']['recombination_volume_policy']='material_local_mixed_voronoi'
        material=folder/'materials.json';material.write_text(json.dumps(mats,indent=2))
        cfg['materials_file']=str(material);cfg['sweep']['initial_state_file']=str(state)
        cfg['output_csv']=str(folder/'iv.csv');cfg['sweep']['write_state_file']=str(folder/'final_state.csv')
        cfg['sweep']['write_state_every_point_prefix']=str(folder/'states/state')
        redirect_diagnostics(cfg,folder)
        config=folder/'config.json';config.write_text(json.dumps(cfg,indent=2));cases[name]=str(config)
    paths={base_path,state,Path(base['materials_file']),Path(base['node_doping_file']),Path(base['mesh_file']),export/'metadata.json',export/'nodes.csv',export/'elements.csv'}
    for region in (0,3):paths.update((export/'fields').glob(f'*_region{region}.csv'))
    manifest={'native':native,'native_thermal_voltage_V':thermal_voltage,'cases':cases,
        'initial_state':str(state),'sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2));return manifest


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('base-config','initial-state','native-export','output'):
        parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--native-thermal-voltage',type=float,required=True)
    args=parser.parse_args()
    print(json.dumps(prepare(args.base_config,args.initial_state,args.native_export,args.output,args.native_thermal_voltage),indent=2))
