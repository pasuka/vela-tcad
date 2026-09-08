"""Numerical evidence gates for the independent column and full-field follow-up."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

REPO=Path(__file__).resolve().parents[2]
ROOT=REPO/'reference_tcad/simplemos_sentaurus2022'


def rows(path):
    with path.open() as f:return list(csv.DictReader(f))


class FollowupEvidenceTests(unittest.TestCase):
    def test_all_target_biases_have_a_strict_state(self):
        points=rows(ROOT/'fullfield_validation_20260906/points.csv')
        expected={(device,vd,round(i*.05,10)) for device in ('n19','n23') for vd in (.05,1.) for i in range(21)}
        self.assertEqual({(r['device'],float(r['vd']),float(r['vg'])) for r in points},expected)
        self.assertEqual(len(points),84)
        for r in points:
            self.assertEqual(r['qualified'],'True')
            self.assertEqual(int(r['local_violations']),0)
            self.assertLessEqual(float(r['vela_kcl_over_Id']),1e-8)
            # Passing the floor-scaled global gate is not source-scale qualification.
            self.assertEqual(r['global_electron_source_qualified'],'False')
            self.assertEqual(r['global_hole_source_qualified'],'False')

    def test_failed_initial_attempts_are_retained(self):
        root=ROOT/'fullfield_validation_20260906'
        old=rows(root/'vela_ascending.csv')+rows(root/'vela_recovery.csv')
        self.assertTrue(any(r['converged']=='True' and r['qualified']=='False' for r in old))
        tighter=rows(root/'vela_reclosure.csv')
        self.assertEqual(len(tighter),2)
        self.assertTrue(all(r['qualified']=='True' for r in tighter))

    def test_srh_unit_reconstruction_matches_production(self):
        rates=rows(ROOT/'fullfield_validation_20260906/srh.csv')
        self.assertEqual(len(rates),84*6)
        for r in rates:self.assertLess(float(r['production_rate_reconstruction_L1']),1e-7)

    def test_jacobian_diagnostics_preserve_the_physical_residual(self):
        root=REPO/'build-release/simplemos_jacobian_columns_20260906'
        digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        for device in ('n19','n23'):
            for vd in ('0p050000','1p000000'):
                case=f'm65_{device}_vd_{vd}_endpoint'
                original=digest(root/case/'baseline_off/columns.csv.residual.csv')
                self.assertEqual(original,digest(root/case/'baseline_on/columns.csv.residual.csv'))
                if vd=='1p000000':
                    for step in ('1e-07','1e-08'):
                        self.assertEqual(original,digest(root/'surface_step'/case/step/'columns.csv.residual.csv'))

    def test_surface_step_defects_converge_under_independent_difference(self):
        controls=rows(ROOT/'jacobian_columns_20260906/surface_step/defect_controls.csv')
        self.assertEqual(len(controls),8)
        for r in controls:
            self.assertGreater(float(r['baseline_relative_error']),1e-3)
            self.assertEqual(r['fd_stable'],'True')
            self.assertLess(float(r['relative_error']),1e-3)
            self.assertLess(float(r['fine_relative_error']),2e-6)
            self.assertEqual(float(r['outside_pattern_relative']),0.)

    def test_unresolved_difference_signals_are_not_promoted(self):
        columns=rows(ROOT/'jacobian_columns_20260906/columns.csv')
        baseline=[r for r in columns if r['variant']=='baseline_on']
        self.assertTrue(any(r['active']=='True' and r['fd_stable']=='False' for r in baseline))
        for r in baseline:
            if r['active']=='False' or r['fd_stable']=='False':self.assertEqual(r['passed'],'False')

    def test_calibrated_poisson_response_closes_at_both_amplitudes(self):
        checks=rows(ROOT/'channel_poisson_coupling_20260906/checks.csv')
        perturbations=[r for r in checks if r['comparison'] in ('full','half')]
        self.assertEqual(len(perturbations),4)
        for r in perturbations:
            self.assertLess(float(r['charge_only_relative']),1e-3)
            self.assertLess(float(r['full_half_psi_relative']),1e-3)
        # A mapped native state is not a solution of the Vela Poisson operator.
        for r in checks:
            if r['comparison']=='absolute_strict_minus_mapped':
                self.assertEqual(r['charge_response_qualified'],'False')
                self.assertGreater(float(r['charge_only_relative']),1.)

    def test_field_visibility_retains_unqualified_carrier_nodes(self):
        fields=rows(ROOT/'fullfield_validation_20260906/native_precision/comparison/visibility/fields.csv')
        self.assertEqual(len(fields),84*2*3)
        self.assertTrue(any(float(r['ignored_squared_error_fraction'])>.9 for r in fields if r['carrier']=='hole'))
        for r in fields:
            self.assertEqual(int(r['selected_nodes']),int(r['qualified_nodes'])+int(r['ignored_nodes']))

    def test_refined_native_currents_meet_the_independent_kcl_gate(self):
        root=ROOT/'fullfield_validation_20260906/native_precision'
        native=rows(root/'native_points.csv')
        self.assertEqual(len(native),42)
        for r in native:
            self.assertEqual(r['native_kcl_qualified'],'True')
            self.assertLessEqual(float(r['kcl_over_Id']),1e-8)
        pairs=rows(root/'comparison/points.csv')
        self.assertEqual(len(pairs),84)
        for r in pairs:
            self.assertEqual(r['qualified'],'True')
            self.assertLessEqual(float(r['native_kcl_over_Id']),1e-8)
            self.assertLessEqual(float(r['vela_kcl_over_Id']),1e-8)

    def test_native_precision_control_keeps_the_physical_field(self):
        stability=rows(ROOT/'fullfield_validation_20260906/native_precision/native_field_stability.csv')
        self.assertEqual(len(stability),42*6)
        for r in stability:
            if 'Potential' in r['field']:self.assertLess(float(r['max_abs_change']),1e-10)


if __name__=='__main__':unittest.main()
