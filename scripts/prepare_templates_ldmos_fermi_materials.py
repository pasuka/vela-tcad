"""Audit a 300 K Silicon parameter source and prepare a consistent material profile.

No curve fitting; historical material contracts are immutable inputs. The
generated intrinsic density excludes BGN, which remains in the live BGN model.
"""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path


def block(source, name):
    source='\n'.join(line for line in source.splitlines() if not line.lstrip().startswith('*'))
    match=re.search(r'(?m)^'+re.escape(name)+r'\s*\{([^}]*)\}',source)
    if not match:raise ValueError(f'Missing {name} parameter block')
    result={}
    for line in match[1].splitlines():
        match=re.match(r'\s*([\w()]+)\s*=\s*([+-]?[\d.]+(?:[eE][+-]?\d+)?)',line)
        if match:result[match[1]]=float(match[2])
    return result


def material_values(source, constants):
    b=block(source,'Bandgap');e=block(source,'eDOSMass');h=block(source,'hDOSMass')
    if e['Formula']!=1 or h['Formula']!=1:raise ValueError('This audit requires DOS Formula 1')
    if b.get('dEg0(OldSlotboom)',0)!=0:raise ValueError('Explicit BGN offset mapping required')
    t=300.;eg=b['Eg0']+b['alpha']*b['Tpar']**2/(b['beta']+b['Tpar'])-b['alpha']*t*t/(b['beta']+t)
    egzero=b['Eg0']+b['alpha']*b['Tpar']**2/(b['beta']+b['Tpar'])
    me=((6*e['a']*egzero/eg)**2*e['ml'])**(1/3)+e['mm']
    numerator=sum(h[k]*t**i for i,k in enumerate(('a','b','c','d','e')))
    denominator=1+sum(h[k]*t**(i+1) for i,k in enumerate(('f','g','h','i')))
    mh=(numerator/denominator)**(2/3)+h['mm']
    prefactor=2*(2*math.pi*constants['m0']*constants['kb']*t/constants['h']**2)**1.5/1e6
    nc=prefactor*me**1.5;nv=prefactor*mh**1.5;vt=constants['kb']/constants['q']*t
    return dict(bandgap_eV=eg,electron_affinity_eV=b['Chi0'],
        conduction_band_density_of_states_cm3=nc,valence_band_density_of_states_cm3=nv,
        intrinsic_carrier_density_cm3=math.sqrt(nc*nv)*math.exp(-eg/(2*vt)),temperature_K=t)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-par',required=True,type=Path)
    parser.add_argument('--base',type=Path,default=Path('reference_tcad/templates_ldmos_sentaurus2022/contracts/materials.json'))
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();root=Path(__file__).resolve().parents[1]
    constants_file=root/'include/vela/core/PhysicalConstants.h'
    constants={k:float(re.search(r'constexpr double '+k+r'\s*=\s*([\d.eE+-]+)',constants_file.read_text())[1]) for k in ('kb','q','h','m0')}
    base=json.loads(args.base.read_text());si=next(m for m in base['materials'] if m['name']=='Si');old=dict(si)
    values=material_values(args.source_par.read_text(),constants);si.update(values)
    si['provenance']='Siliconc100.par DOS Formula 1 and Eg(T) at 300 K; current Vela SI constants; ni=sqrt(Nc*Nv)*exp(-Eg/(2*kT)); BGN applied separately, no curve fit'
    base['revision']+=1
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(base,indent=2)+'\n')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    audit=dict(temperature_K=300,constants=constants,revision=base['revision'],
        source_files={str(p):digest(p) for p in (args.source_par,args.base,constants_file)},
        changes={k:dict(old=old[k],new=v,relative_change=v/old[k]-1) for k,v in values.items()},
        output_sha256=digest(args.output),note='Historical input contract preserved; curve acceptance is a separate requirement.')
    args.output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit['changes'],indent=2))


if __name__=='__main__':main()
