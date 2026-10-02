"""Independent manual-constant replay and exported element support checks.

No physical constant is fitted to residuals. The original replay is immutable;
this addendum tests the constants printed in the 2022 manual's built-in PMI
example and retains both the original and replayed residual.
"""
from pathlib import Path
import math
import numpy as np
import audit_simplemos_masetti_native_geometry_20260907 as g

p=g.p;a=p.a;d=p.d;OUT=p.OUT;LOCAL=p.LOCAL
NATIVE_Q=1.602192e-19
NATIVE_EPS0=8.8542e-12


def main():
    a.verify(OUT/'native_geometry_audit.json');a.verify(OUT/'analysis_scope.json')
    manual=p.MANUAL.read_text(encoding='utf-8')
    assert 'const double eps0 = 8.8542e-14;' in manual and 'const double e0 = 1.602192e-19;' in manual
    results=[];element_checks=[];projection=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for device in ('n19','n23'):
        geo,xy,elements,volumes,K,parts,info,modified=g.geometry(device)
        nativeK=K*(NATIVE_EPS0/d.matrix.spatial.m73.EPS0)
        # BM_ElementVolume is an independently named element field, so compare
        # every exported cell with its coordinate triangle, not only array size.
        src=LOCAL/'native_exports'/device
        for region,material in [(0,'Si'),(1,'SiO2'),(2,'Nitride'),(3,'Nitride')]:
            rows=a.rows(src/f'fields/BM_ElementVolume_region{region}_cells.csv')
            errors=[]
            for row in rows:
                cell=int(row['cell_id']);e=elements[cell];assert e['material']==material
                area=d.matrix.spatial.m73.triangle_area([xy[n] for n in e['nodes']])
                errors.append(abs(float(row['component0'])/area-1))
            assert max(errors)<=1e-10
            element_checks.append(dict(device=device,region=region,cells=len(rows),volume_to_coordinate_area_max_relative=max(errors),passed=True))
        for c in a.read(OUT/'vela_contract.json')['cases']:
            if c['device']!=device:continue
            raw=p.prev.LOCAL/'native_exports/masetti'/c['case']/f"vg_{c['index']:03d}"
            psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(raw,geo);assert spread<=1e-12
            nd=g.sc(raw/'fields/DonorConcentration_region0.csv');na=g.sc(raw/'fields/AcceptorConcentration_region0.csv')
            charge=n-h-np.array([(nd.get(i,0)-na.get(i,0))*1e6 for i in range(geo.count)])
            nativeR=nativeK@psi+NATIVE_Q*charge*volumes['Si']
            charge_scale=np.linalg.norm((NATIVE_Q*charge*volumes['Si'])[geo.free])
            relative=np.linalg.norm(nativeR[geo.free])/charge_scale
            # Residual replay is an independent numerical check, not an exported
            # native residual and not a self-consistent candidate A/B gate.
            results.append(dict(key=c['key'],native_constant_replay_over_charge=float(relative),
                replay_below_1e_minus_10=bool(relative<=1e-10),native_constants_source='2022 built-in Schottky resistance PMI example; no fitted value'))
            mapped=d.array(d.ordered(Path(c['mapped']),geo.count),('psi','electrons_m3','holes_m3'))
            root=LOCAL/'vela'/c['key'];mappedR=d.array(d.ordered(root/'mapped/residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            strictR=d.array(d.ordered(root/'strict/residual.csv',geo.count),('psi_residual','phin_residual','phip_residual'))
            weights=d.array(d.ordered(root/'adjoint.csv',geo.count),('lambda_poisson','lambda_electron','lambda_hole'))
            sources={
                'native_constant_replay':nativeR,
                'epsilon0_convention':(K-nativeK)@psi,
                'charge_constant_convention':(d.fixed.Q-NATIVE_Q)*charge*volumes['Si'],
                'dielectric_geometry':(geo.matrices['legacy']-K)@psi,
                'charge_volume':d.fixed.Q*charge*(geo.volumes['all_cell']-volumes['Si']),
                'thermal_density_convention':d.fixed.Q*((mapped[1]-n)-(mapped[2]-h))*geo.volumes['all_cell']}
            sources['replay_remainder']=mappedR[0]/factor-sum(sources.values())
            for name,source in sources.items():
                value=-math.fsum(float(x*y) for x,y in zip(weights[0],source*factor))
                response=geo.solve('legacy',source)
                projection.append(dict(key=c['key'],component=name,screening_current_A_per_um=value,
                    conditional_interface_338_psi_V=float(response[338]),screening_only=True,
                    continuity_channel_calibration_does_not_qualify_this_source=True))
            for block,name in [(1,'mapped_electron_residual'),(2,'mapped_hole_residual')]:
                value=-math.fsum(float(x*y) for x,y in zip(weights[block],mappedR[block]))
                projection.append(dict(key=c['key'],component=name,screening_current_A_per_um=value,
                    conditional_interface_338_psi_V='',screening_only=True,continuity_channel_calibration_does_not_qualify_this_source=True))
            projection.append(dict(key=c['key'],component='minus_strict_residual',screening_current_A_per_um=-d.project(weights,strictR),
                conditional_interface_338_psi_V='',screening_only=True,continuity_channel_calibration_does_not_qualify_this_source=True))
            total=math.fsum(x['screening_current_A_per_um'] for x in projection if x['key']==c['key'])
            expected=d.project(weights,mappedR-strictR)
            assert abs(total-expected)<=max(1e-12*abs(expected),1e-24)
    a.write_csv(OUT/'native_poisson_constant_replay.csv',results)
    a.write_csv(OUT/'native_element_support_check.csv',element_checks)
    a.write_csv(OUT/'poisson_transport_screening.csv',projection)
    a.write(OUT/'poisson_constants_addendum.json',dict(status='completed_read_only_addendum',
        native_charge_C=NATIVE_Q,native_epsilon0_F_per_m=NATIVE_EPS0,constants_fitted=False,
        interpretation='The prescribed manual constants are tested by eight residual replays. Success supports this mapped Poisson expression; it does not prove a self-consistent current correction.',
        original_all_acute_geometry_identity_attempt='Failed before any geometry output; the retained geometry audit uses the independently established 40-cell Si/SiO2 numbering witness and records all nonmatching cells.',
        input_hashes={a.rel(f):a.sha(f) for f in [Path(__file__).resolve(),p.MANUAL,OUT/'native_geometry_audit.json',OUT/'analysis_scope.json']}))
    print('Native constant Poisson replays',sum(r['replay_below_1e_minus_10'] for r in results),'/ 8; element support checks',len(element_checks),flush=True)


if __name__=='__main__':main()
