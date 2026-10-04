"""Reconstruct the interface coefficients from geometry and source formulas."""
from pathlib import Path
import math
import localize_simplemos_remaining_poisson_20260907 as l

a=l.a;d=l.d;OUT=l.OUT

def main():
    a.verify(OUT/'ledger_evidence.json');cases=a.read(OUT/'ledger_contract.json')['cases'];rows=[];summary=[];files=[Path(__file__).resolve(),OUT/'ledger_evidence.json']
    for device in ('n19','n23'):
        geo,xy,elements,volumes,K,parts,info,_=l.g.geometry(device)
        c=next(c for c in cases if c['device']==device);cfg=l.v.config(c)
        materials=a.read(Path(cfg['materials_file']));eps={m['name']:m['eps_r']*d.matrix.spatial.m73.EPS0 for m in materials['materials']}
        for (i,j),part in parts.items():
            if {x['material'] for x in part}!={'Si','SiO2'}:continue
            assert len(part)==2
            raw=[]
            for x in part:
                e=elements[x['cell']];nodes=e['nodes'];k=next(k for k in range(3) if nodes[k] not in (i,j))
                coefficient=l.g.box.raw_coefficient([xy[n] for n in nodes],k)
                raw.append(dict(material=x['material'],coefficient=coefficient,native=x['coefficient']))
            gs={x['material']:x['coefficient'] for x in raw}
            old_expected=.5*(eps['Si']+eps['SiO2'])*(gs['Si']+gs['SiO2'])
            native_expected=eps['Si']*gs['Si']+eps['SiO2']*gs['SiO2']
            old=float(-geo.matrices['legacy'][i,j]);new=float(-K[i,j]);region=float(-geo.matrices['region_local'][i,j])
            # For two adjacent regions, the averaging defect has this exact factorization.
            defect=.5*(eps['Si']-eps['SiO2'])*(gs['SiO2']-gs['Si'])
            error_old=abs(old-old_expected)/abs(old);error_new=abs(new-native_expected)/abs(new)
            error_defect=abs((old-new)-defect)/max(abs(defect),1e-300)
            error_local=abs(region/new-1)
            ok=error_old<=1e-10 and error_new<=1e-10 and error_defect<=1e-8 and error_local<=1e-10
            rows.append(dict(device=device,node0=i,node1=j,x0_um=xy[i][0],y0_um=xy[i][1],x1_um=xy[j][0],y1_um=xy[j][1],
                Si_g=gs['Si'],SiO2_g=gs['SiO2'],old_K_F_per_m=old,native_geometry_K_F_per_m=new,
                old_over_native_minus_one=old/new-1,old_formula_relative=error_old,native_formula_relative=error_new,
                defect_identity_relative=error_defect,region_local_to_native_relative=error_local,qualified=ok))
        for c in (c for c in cases if c['device']==device):
            spatial=[r for r in a.rows(OUT/'dielectric_spatial.csv') if r['key']==c['key'] and r['state']=='native_mu']
            total=math.fsum(float(r['signed_prediction_A_per_um']) for r in spatial)
            interface=float(next(r['signed_prediction_A_per_um'] for r in spatial if r['materials']=='Si;SiO2'))
            summary.append(dict(key=c['key'],interface_signed_prediction_A_per_um=interface,total_signed_prediction_A_per_um=total,
                interface_over_total=interface/total,other_regions_sum_A_per_um=total-interface,interpretation='Signed contributions can cancel; fraction may exceed one. Same-direction DC calibration is reported separately.'))
        files += [Path(cfg['materials_file']),l.g.LOCAL/'native_exports'/device/'nodes.csv',l.g.LOCAL/'native_exports'/device/'elements.csv',l.g.LOCAL/'native_raw/bundle'/device/'MeasureCoefficients.debug']
    a.write_csv(OUT/'interface_formula.csv',rows);a.write_csv(OUT/'interface_share.csv',summary)
    files += [OUT/'interface_formula.csv',OUT/'interface_share.csv',l.p.REPO/'include/vela/equation/AssemblerUtils.h',l.p.REPO/'src/equation/CoupledDDAssembler.cpp',Path(l.g.__file__).resolve(),Path(l.g.box.__file__).resolve()]
    d.matrix.freeze(OUT/'interface_formula_evidence.json',files)
    print('Interface averaging formula identities',sum(r['qualified'] for r in rows),'/',len(rows),flush=True)
    print('Legacy coefficient excess:',min(r['old_over_native_minus_one'] for r in rows),max(r['old_over_native_minus_one'] for r in rows),flush=True)

if __name__=='__main__':main()
