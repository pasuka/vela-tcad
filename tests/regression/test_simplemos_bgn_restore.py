"""Guard the isolation and evidence semantics of the BGN restoration runner."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
import simplemos_bgn_restore_native_20260908 as n
import simplemos_bgn_restore_vela_20260908 as v
import simplemos_bgn_intrinsic_control_20260908 as ni


class BgnRestorationContract(unittest.TestCase):
    def test_native_only_model_keyword_changes(self):
        header = ('Physics { EffectiveIntrinsicDensity(NoBandGapNarrowing) }\n'
                  'Physics(Material="Silicon") { Mobility(DopingDependence) '
                  'Recombination(SRH(DopingDependence)) }\n' + n.prior.prior.MATH + '\n')
        off, on = [n.deck(header, model) for model in n.MODELS]
        self.assertEqual(off.replace('EffectiveIntrinsicDensity(NoBandGapNarrowing)',
                         'EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom))'), on)
        self.assertEqual(on.count('Load(FilePrefix='), 4)
        self.assertNotIn('Quasistationary', on)
        self.assertIn('vg_000', on)
        self.assertIn('vg_050', on)

    def test_vela_preserves_all_other_configuration(self):
        cfg = dict(materials_file='frozen.json', mesh_geometry=dict(cell_box_policy='delaunay_transfer'),
                   solver=dict(mobility=dict(model='masetti', doping_concentration_basis='total_impurity',
                                           edge_averaging='element_box'),
                               bandgap_narrowing=dict(model='none', fermi_statistics_correction=False,
                                                    equal_ni_flux_evaluation='compensated_log_expm1'),
                               carrier_row_convergence=dict(eps_row=1e-6), linear_refinement_iterations=4))
        original = copy.deepcopy(cfg)
        on = v.model_config(cfg, 'old_slotboom')
        self.assertEqual(cfg, original)
        on['solver']['bandgap_narrowing']['model'] = 'none'
        self.assertEqual(on, cfg)
        self.assertEqual(v.model_config(cfg, 'no_bgn'), cfg)

    def test_invalid_baseline_cannot_be_silently_relabelled(self):
        cfg = dict(solver=dict(mobility=dict(model='phumob', doping_concentration_basis='total_impurity',
                                            edge_averaging='element_box'), bandgap_narrowing=dict(model='none')))
        with self.assertRaises(AssertionError):
            v.model_config(cfg, 'old_slotboom')
        with self.assertRaises(AssertionError):
            n.deck('Physics { EffectiveIntrinsicDensity(OldSlotboom) }', 'no_bgn')

    def test_failed_initialization_is_not_qualified_by_small_current_error(self):
        delta = dict(psi_max_V=1e-10, phin_max_V=1e-10, phip_max_V=1e-10, density_max_relative=1e-9)
        self.assertFalse(v.v.dual_qualified(True, False, delta, 0.))
        delta['phip_max_V'] = 1.01e-6
        self.assertFalse(v.v.dual_qualified(True, True, delta, 0.))

    def test_intrinsic_control_rejects_inconsistent_exported_convention(self):
        rows = [dict(model='old_slotboom', source_identity_max_V=1e-15,
                     bgn_max_abs_error_eV=1e-17, inferred_base_ni_min_cm3=1.46e10,
                     inferred_base_ni_max_cm3=1.46e10) for _ in range(16)]
        self.assertEqual(ni.native_intrinsic(rows), 1.46e10)
        rows[3]['inferred_base_ni_max_cm3'] *= 1.0001
        with self.assertRaises(AssertionError):
            ni.native_intrinsic(rows)


if __name__ == '__main__':
    unittest.main()
