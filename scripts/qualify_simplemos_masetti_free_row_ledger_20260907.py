"""Separate physical free Poisson rows from replaced Dirichlet rows.

The earlier raw screening sum was algebraically correct, but its separate
physical-contact terms canceled and were not valid equation contributions.
Keep that file unchanged and supersede its component interpretation here.
"""
from pathlib import Path
import math
import numpy as np
import audit_simplemos_masetti_poisson_constants_20260907 as s

p=s.p;a=p.a;d=p.d;OUT=p.OUT;LOCAL=p.LOCAL


def main():
    a.verify(OUT/'poisson_constants_addendum.json');rows=[];checks=[]
    factor=a.read(d.matrix.p.CONTRACT)['frozen_poisson_residual_per_physical_charge']
    for device in ('n19','n23'):
        geo,xy,elements,volumes,K,parts,info,modified=s.g.geometry(device)
        nativeK=K*(s.NATIVE_EPS0/d.matrix.spatial.m73.EPS0)
        for c in a.read(OUT/'vela_contract.json')['cases']:
            if c['device']!=device:continue
            src=p.prev.LOCAL/'native_exports/masetti'/c['case']/f"vg_{c['index']:03d}"
            psi,n,h,spread=d.matrix.spatial.m73.sentaurus_state(src,geo)
            nd=s.g.sc(src/'fields/DonorConcentration_region0.csv');na=s.g.sc(src/'fields/AcceptorConcentration_region0.csv')
            charge=n-h-np.array([(nd.get(i,0)-na.get(i,0))*1e6 for i in range(geo.count)])
            root=LOCAL/'vela'/c['key'];mapped=d.array(d.ordered(Path(c['mapped']),geo.count),('psi','electrons_m3','holes_m3'))
            residual={r:d.array(d.ordered(root/r/'residual.csv',geo.count),('psi_residual','phin_residual','phip_residual')) for r in ('mapped','strict')}
            weights=d.array(d.ordered(root/'adjoint.csv',geo.count),('lambda_poisson','lambda_electron','lambda_hole'))
            sources={
                'native_constant_replay':nativeK@psi+s.NATIVE_Q*charge*volumes['Si'],
                'epsilon0_convention':(K-nativeK)@psi,
                'charge_constant_convention':(d.fixed.Q-s.NATIVE_Q)*charge*volumes['Si'],
                'dielectric_geometry':(geo.matrices['legacy']-K)@psi,
                'charge_volume':d.fixed.Q*charge*(geo.volumes['all_cell']-volumes['Si']),
                'thermal_density_convention':d.fixed.Q*((mapped[1]-n)-(mapped[2]-h))*geo.volumes['all_cell']}
            sources['replay_remainder']=residual['mapped'][0]/factor-sum(sources.values())
            total=0.
            for name,source in sources.items():
                vector=np.zeros_like(residual['mapped']);vector[0,geo.free]=source[geo.free]*factor
                value=d.project(weights,vector);total+=value
                rows.append(dict(key=c['key'],component=name,screening_current_A_per_um=value,screening_only=True,source_support='free_poisson_rows'))
            for block,name in [(1,'mapped_electron_residual'),(2,'mapped_hole_residual')]:
                vector=np.zeros_like(residual['mapped']);vector[block]=residual['mapped'][block]
                value=d.project(weights,vector);total+=value
                rows.append(dict(key=c['key'],component=name,screening_current_A_per_um=value,screening_only=True,source_support='actual_mapped_residual_rows'))
            contact=np.zeros_like(residual['mapped']);contact[0,geo.contact_nodes]=residual['mapped'][0,geo.contact_nodes]
            for name,vector in [('mapped_poisson_dirichlet_rows',contact),('minus_strict_residual',-residual['strict'])]:
                value=d.project(weights,vector);total+=value
                rows.append(dict(key=c['key'],component=name,screening_current_A_per_um=value,screening_only=True,source_support='actual_replaced_equation_rows'))
            expected=d.project(weights,residual['mapped']-residual['strict'])
            error=abs(total-expected)/max(abs(expected),1e-300);assert error<=1e-10
            checks.append(dict(key=c['key'],projection_sum_relative_error=error,passed=True))
    a.write_csv(OUT/'poisson_transport_free_row_screening.csv',rows)
    a.write_csv(OUT/'free_row_screening_checks.csv',checks)
    a.write(OUT/'free_row_screening_scope.json',dict(status='completed_screening_with_explicit_dirichlet_rows',
        supersedes_component_interpretation_of='poisson_transport_screening.csv',
        prior_issue='Physical contact-node Poisson charges are not Dirichlet equation residuals. The raw total retained cancellation, but native/remainder component projections were not usable. The prior file remains diagnostic and is not used for attribution.',
        limitations='Channel electron source calibration does not qualify these distinct Poisson or distributed continuity sources, nor a finite model replacement.',
        input_hashes={a.rel(x):a.sha(x) for x in [Path(__file__).resolve(),OUT/'poisson_constants_addendum.json',OUT/'poisson_transport_screening.csv']}))
    print('Eight free-row screening identities passed; replaced contact rows recorded separately',flush=True)


if __name__=='__main__':main()
