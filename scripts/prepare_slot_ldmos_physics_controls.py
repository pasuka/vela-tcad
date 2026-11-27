"""Prepare independent SLOT Auger, band and poly controls; never edit a baseline."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

from run_slot_ldmos_noavalanche_continuation import redirect_diagnostics
from audit_slot_ldmos_element_ialmob import validate_export_mesh


def native_values(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    parts = re.split(r'Material\s*=\s*"([^"]+)"\s*\{', text)
    materials = dict(zip(parts[1::2], parts[2::2]))
    def scalar(body: str, key: str) -> float:
        match = re.search(r'^\s*' + re.escape(key) + r'\s*=\s*([-+\d.eE]+)', body, re.M)
        if not match:
            raise ValueError(f"Missing native parameter {key}")
        return float(match[1])
    si = materials["Silicon"]
    auger = re.search(r'^Auger\s[^\n]*\n\{(.*?)^\}', si, re.M | re.S)
    if not auger:
        raise ValueError("Missing native Auger block")
    result = {"Si": {"bandgap_eV": scalar(si, "Eg0"), "electron_affinity_eV": scalar(si, "Chi0")},
              "SiO2": {"electron_affinity_eV": scalar(materials["Oxide"], "Chi0")}}
    pairs = {}
    for key in ("A", "B", "C", "H", "N0"):
        match = re.search(r'^\s*' + key + r'\s*=\s*([-+\d.eE]+)\s*,\s*([-+\d.eE]+)', auger[1], re.M)
        if not match:
            raise ValueError(f"Missing native Auger {key}")
        pairs[key] = [float(match[1]), float(match[2])]
    result["auger"] = pairs
    return result


def prepare(base_path: Path, state: Path, parameter_path: Path, output: Path) -> dict:
    base_path, state, parameter_path, output = [p.resolve() for p in (base_path, state, parameter_path, output)]
    if output.exists():
        raise ValueError("Refusing to overwrite a control bundle")
    base = json.loads(base_path.read_text(encoding="utf-8"))
    for key in ('materials_file','mesh_file','node_doping_file'):
        if key in base:
            p=Path(base[key])
            base[key]=str(p.resolve() if p.is_absolute() else (base_path.parent/p).resolve())
    if base.get("scaling", {}).get("mode") != "unit_scaling":
        raise ValueError("SLOT controls require explicit cm-based unit_scaling")
    if base["solver"].get("temperature_K", 300.) != 300.:
        raise ValueError("Native A+B+C controls are restricted to 300 K")
    if base["solver"].get("auger_with_generation") is not False:
        raise ValueError("Expected the qualified plain-Auger baseline")
    if base["solver"].get("auger_density_dependence", {}).get("enabled", False):
        raise ValueError("Expected a baseline without Auger density enhancement")
    if not state.is_file():
        raise ValueError("Missing qualified seed")
    material_path = Path(base["materials_file"])
    if not material_path.is_absolute():
        material_path = base_path.parent / material_path
    materials = json.loads(material_path.read_text(encoding="utf-8"))
    native = native_values(parameter_path)
    for i, key in enumerate(("auger_cn_m6_per_s", "auger_cp_m6_per_s")):
        expected = sum(native["auger"][name][i] for name in ("A", "B", "C"))
        if abs(base["solver"][key] / expected - 1.) > 1e-12:
            raise ValueError("Baseline Auger coefficient does not match native 300 K")
    output.mkdir(parents=True)
    cases = {}
    for name in ("baseline", "auger_density", "si_bandgap", "si_affinity", "oxide_affinity", "band_parameters"):
        folder = output / name
        folder.mkdir()
        cfg, mats = copy.deepcopy(base), copy.deepcopy(materials)
        cfg["_comment"] = "Independent SLOT physics control: " + name
        cfg["sweep"]["initial_state_file"] = str(state)
        cfg["output_csv"] = str(folder / "iv.csv")
        cfg["sweep"]["write_state_file"] = str(folder / "final_state.csv")
        cfg["sweep"]["write_state_every_point_prefix"] = str(folder / "states/state")
        if "boundary_control" in cfg["sweep"]:
            cfg["sweep"].pop("boundary_control")
        redirect_diagnostics(cfg, folder)
        if name == "auger_density":
            cfg["solver"]["auger_density_dependence"] = {"enabled": True}
            for i, carrier in enumerate(("electron", "hole")):
                cfg["solver"]["auger_density_dependence"][carrier] = {
                    "enhancement": native["auger"]["H"][i],
                    "reference_density_m3": native["auger"]["N0"][i]}
        changes = {}
        if name in ("si_bandgap", "band_parameters"):
            changes.setdefault("Si", {})["bandgap_eV"] = native["Si"]["bandgap_eV"]
        if name in ("si_affinity", "band_parameters"):
            changes.setdefault("Si", {})["electron_affinity_eV"] = native["Si"]["electron_affinity_eV"]
        if name in ("oxide_affinity", "band_parameters"):
            changes["SiO2"] = native["SiO2"]
        if changes:
            for material, patch in changes.items():
                row = next((r for r in mats["materials"] if r["name"] == material), None)
                if row is None:
                    row = {"name": material}
                    mats["materials"].append(row)
                row.update(patch)
            cfg["materials_file"] = str(folder / "materials.json")
            (folder / "materials.json").write_text(json.dumps(mats, indent=2), encoding="utf-8")
        config_path = folder / "config.json"
        config_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
        cases[name] = str(config_path)
    manifest = {"native_parameters": native, "cases": cases, "initial_state": str(state),
                "sources": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in (base_path, state, parameter_path, material_path)},
                "note": "Single-factor controls; ni, mesh, contacts, solver gates and baseline files preserved."}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def prepare_poly(base_path: Path, state: Path, fields: Path, output: Path) -> dict:
    """Restore an isolated poly region without applying Silicon-only sources.

    Intrinsic density is inferred from the native Boltzmann identity, with the
    fitted native thermal voltage recorded separately from Vela's constant.
    This is a controlled material restoration, not a claim that intrinsic
    energy references or all material-local volumes are already equivalent.
    """
    base_path,state,output=[p.resolve() for p in (base_path,state,output)]
    if output.exists(): raise ValueError("Refusing to overwrite poly controls")
    base=json.loads(base_path.read_text(encoding="utf-8"))
    for key in ('materials_file','mesh_file','node_doping_file'):
        p=Path(base[key]);base[key]=str(p if p.is_absolute() else (base_path.parent/p).resolve())
    if base.get("scaling",{}).get("mode")!="unit_scaling":
        raise ValueError("Poly control requires unit_scaling")
    if base["solver"].get("temperature_K",300.)!=300.:
        raise ValueError("Native poly controls are restricted to 300 K")
    impact=base["solver"].get("impact_ionization","none")
    if (impact.get('model','none') if isinstance(impact,dict) else impact)!="none":
        raise ValueError("This independent poly control starts without avalanche")
    fields=fields.resolve()
    metadata=fields.parent/'metadata.json'
    units={(f['name'],f['region']):f['unit'] for f in json.loads(metadata.read_text())['fields']}
    required_units={**{k:'cm^-3' for k in ('eDensity','hDensity','DonorConcentration','AcceptorConcentration')},
        **{k:'V' for k in ('eQuasiFermiPotential','hQuasiFermiPotential')},
        **{k:'eV' for k in ('BandgapNarrowing','BandGap','ElectronAffinity')},
        **{k:'cm^2*V^-1*s^-1' for k in ('eMobility','hMobility')}}
    for name,unit in required_units.items():
        if units.get((name,3))!=unit:raise ValueError('Unexpected native poly unit: '+name)
    dependencies={metadata,fields.parent/'nodes.csv',fields.parent/'elements.csv'}
    def read(name):
        path=fields/f'{name}_region3.csv'
        if not path.exists():path=fields/f'{name}_region3_cells.csv'
        with path.open(encoding="utf-8",newline='') as stream:
            rows=list(csv.DictReader(stream))
        result={int(r.get('node_id',r.get('cell_id'))):float(r['component0']) for r in rows}
        if not result or len(result)!=len(rows) or not all(math.isfinite(v) for v in result.values()):
            raise ValueError('Invalid native poly field: '+name)
        dependencies.add(path)
        return result
    values={name:read(name) for name in ('eDensity','hDensity','eQuasiFermiPotential','hQuasiFermiPotential',
        'DonorConcentration','AcceptorConcentration','BandgapNarrowing','BandGap','ElectronAffinity')}
    ids=set(values['eDensity'])
    if any(set(v)!=ids for v in values.values()):raise ValueError("Native poly field node mismatch")
    mesh=json.loads(Path(base['mesh_file']).read_text())
    validate_export_mesh(mesh,fields.parent)
    regions={r['id']:r['material'] for r in mesh['regions']}
    actual={i for c in mesh['triangles'] if regions[c['region_id']]=='PolySilicon' for i in c['node_ids']}
    if ids!=actual:raise ValueError("Poly field support does not match the local mesh")
    x=[values['BandgapNarrowing'][i]+values['hQuasiFermiPotential'][i]-values['eQuasiFermiPotential'][i] for i in ids]
    y=[.5*math.log(values['eDensity'][i]*values['hDensity'][i]) for i in ids]
    xm,ym=statistics.mean(x),statistics.mean(y)
    slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x)
    ni=math.exp(ym-slope*xm)
    fit_error=max(abs(b-(ym+slope*(a-xm))) for a,b in zip(x,y))
    if fit_error>1e-8:raise ValueError("Native poly state does not fit the Boltzmann/OldSlotboom contract")
    mobilities={k:statistics.mean(read(k).values()) for k in ('eMobility','hMobility')}
    for k,mean in mobilities.items():
        if any(abs(v/mean-1)>1e-12 for v in read(k).values()):raise ValueError("Expected native constant poly mobility")
    mats=json.loads(Path(base['materials_file']).read_text())
    poly=next(r for r in mats['materials'] if r['name']=='PolySilicon')
    poly.update(ni=ni,mun=mobilities['eMobility'],mup=mobilities['hMobility'],
        bandgap_eV=statistics.mean(values['BandGap'].values()),
        electron_affinity_eV=statistics.mean(values['ElectronAffinity'][i]-.5*values['BandgapNarrowing'][i] for i in ids),
        Nc_m3=2.8e19,Nv_m3=1.04e19,temperature_K=300.)
    with Path(base['node_doping_file']).open(encoding="utf-8",newline='') as stream:
        doping=list(csv.DictReader(stream))
    with state.open(encoding="utf-8",newline='') as stream:
        states=list(csv.DictReader(stream))
    if not ids.issubset({int(r['node_id']) for r in doping}):raise ValueError("Incomplete local doping support")
    output.mkdir(parents=True)
    cases={}
    vt=1.380649e-23*300/1.602176634e-19
    for name in ('poly_intrinsic','poly_doped'):
        folder=output/name;folder.mkdir()
        cfg=copy.deepcopy(base);cfg['_comment']='Independent no-avalanche poly restoration: '+name
        cfg['solver']['recombination_materials']=['Si']
        cfg['solver']['mobility']['ialmob']['constant_material_mobility_cm2_per_Vs']={
            'PolySilicon':{'electron':mobilities['eMobility'],'hole':mobilities['hMobility']}}
        for contact in cfg['contacts']:
            if contact['name'] in ('gate','SLOT'):
                contact['type']='ohmic';contact.pop('flatband_voltage',None)
        rows=copy.deepcopy(doping);seed=copy.deepcopy(states)
        for row in rows:
            i=int(row['node_id'])
            if i in ids:
                row['donors_cm3']=values['DonorConcentration'][i] if name=='poly_doped' else 0.
                row['acceptors_cm3']=values['AcceptorConcentration'][i] if name=='poly_doped' else 0.
        for row in seed:
            i=int(row['node_id'])
            if i not in ids:continue
            for key in ('phin','phip','electron_qf_increment_V','hole_qf_increment_V','electron_qf_reference_V','hole_qf_reference_V'):
                if key in row:row[key]=0.
            nd,na=(values['DonorConcentration'][i],values['AcceptorConcentration'][i]) if name=='poly_doped' else (0.,0.)
            bgn=.009*(math.log((nd+na)/1e17)+math.sqrt(math.log((nd+na)/1e17)**2+.5)) if nd+na>0 else 0.
            effective=ni*math.exp(bgn/(2*vt))
            psi=vt*math.log(values['eDensity'][i]/effective) if name=='poly_doped' else 0.
            row['psi']=psi;row['electrons_m3']=effective*math.exp(psi/vt)*1e6;row['holes_m3']=effective*math.exp(-psi/vt)*1e6
        for filename,rows_ in [('doping.csv',rows),('seed.csv',seed)]:
            with (folder/filename).open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows_[0]));writer.writeheader();writer.writerows(rows_)
        (folder/'materials.json').write_text(json.dumps(mats,indent=2),encoding='utf-8')
        cfg['materials_file']=str(folder/'materials.json');cfg['node_doping_file']=str(folder/'doping.csv')
        cfg['sweep']['initial_state_file']=str(folder/'seed.csv');cfg['sweep'].pop('boundary_control',None)
        cfg['output_csv']=str(folder/'iv.csv');cfg['sweep']['write_state_file']=str(folder/'final_state.csv')
        cfg['sweep']['write_state_every_point_prefix']=str(folder/'states/state')
        redirect_diagnostics(cfg,folder)
        (folder/'config.json').write_text(json.dumps(cfg,indent=2),encoding='utf-8');cases[name]=str(folder/'config.json')
    manifest={'cases':cases,'poly_nodes':len(ids),'native_thermal_voltage_V':1/(2*slope),'ni_cm3':ni,
        'max_log_fit_error':fit_error,'mobility_cm2_per_Vs':mobilities,'source_fields':str(fields),
        'qualification':'no-avalanche material restoration; energy reference and volume parity remain explicit controls'}
    manifest['sources']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(dependencies|{base_path,state,Path(base['mesh_file']),Path(base['materials_file']),Path(base['node_doping_file'])})}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base-config", "initial-state", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--parameters',type=Path)
    source.add_argument('--poly-fields',type=Path)
    args = parser.parse_args()
    result=(prepare_poly(args.base_config,args.initial_state,args.poly_fields,args.output) if args.poly_fields
            else prepare(args.base_config,args.initial_state,args.parameters,args.output))
    print(json.dumps(result,indent=2))
