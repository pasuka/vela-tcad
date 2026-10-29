import sys
import json
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from analyze_bvmethods_mobility_loadline import jvp_block_metrics, endpoint_probe_eligible, audit_accepted_point, audit_saved_state

class MobilityLoadlineAuditContracts(unittest.TestCase):
    def test_tiny_display_error_with_zero_fd_is_not_a_derivative_pass(self):
        row = {'analytic_phip_norm': '4e-50', 'finite_difference_phip_norm': '0', 'phip_relative_error': '4e-50'}
        result = jvp_block_metrics(row, 'phip')
        self.assertEqual(result['relative_error'], 1)
        self.assertFalse(result['finite_difference_signal_resolved'])

    def test_small_resolved_signal_uses_its_own_norm(self):
        row = {'analytic_phin_norm': '2e-9', 'finite_difference_phin_norm': '2.01e-9', 'phin_relative_error': '1e-11'}
        self.assertAlmostEqual(jvp_block_metrics(row, 'phin')['relative_error'], 1 / 201)

    def test_nonfinite_derivative_metric_is_not_accepted(self):
        with self.assertRaises(ValueError):
            jvp_block_metrics({'analytic_phin_norm': 'nan', 'finite_difference_phin_norm': '1',
                               'phin_relative_error': '0'}, 'phin')

    def test_boltzmann_endpoint_ledger_is_ineligible_for_fermi_physics(self):
        self.assertFalse(endpoint_probe_eligible({'solver': {'carrier_statistics': 'fermi_dirac'}}))
        self.assertTrue(endpoint_probe_eligible({'solver': {'carrier_statistics': 'boltzmann'}}))

    def test_rejected_csv_current_placeholders_are_refused(self):
        with self.assertRaises(ValueError):
            audit_accepted_point({'converged': '0'})

    def test_accepted_csv_is_independently_checked_against_resistor(self):
        row = {'converged': '1', 'inner_voltage_V': '4e-6', 'outer_voltage_V': '1e-5',
               'current_total_A_per_um': '6e-13', 'carrier_row_max_ratio': '0.0009',
               'global_electron_continuity_closure_ratio': '0.001', 'global_hole_continuity_closure_ratio': '0.001',
               'contact_drain_total_A_per_um': '6e-13', 'contact_source_total_A_per_um': '-6e-13',
               'contact_gate_total_A_per_um': '0', 'contact_substrate_total_A_per_um': '0'}
        self.assertTrue(audit_accepted_point(row)['pass'])
        row['outer_voltage_V'] = '.01'
        self.assertFalse(audit_accepted_point(row)['pass'])

    def test_reused_seed_cannot_claim_fresh_carrier_and_continuity_gates(self):
        row = {'converged': '1', 'solver_method': 'boundary_seed_state',
               'inner_voltage_V': '4e-6', 'outer_voltage_V': '1e-5',
               'current_total_A_per_um': '6e-13',
               'contact_drain_total_A_per_um': '6e-13', 'contact_source_total_A_per_um': '-6e-13',
               'contact_gate_total_A_per_um': '0', 'contact_substrate_total_A_per_um': '0'}
        result = audit_accepted_point(row)
        self.assertTrue(result['seed_gates_require_parent_evidence'])
        self.assertIsNone(result['checks']['carrier_rows'])
        self.assertFalse(result['pass'])

    def test_coupled_default_zero_diagnostics_are_not_physical_gate_evidence(self):
        row = {'converged': '1', 'solver_method': 'coupled_load_line_newton',
               'inner_voltage_V': '4e-6', 'outer_voltage_V': '1e-5', 'current_total_A_per_um': '6e-13',
               'carrier_row_max_ratio': '0', 'global_electron_continuity_closure_ratio': '0',
               'global_hole_continuity_closure_ratio': '0',
               'contact_drain_total_A_per_um': '6e-13', 'contact_source_total_A_per_um': '-6e-13',
               'contact_gate_total_A_per_um': '0', 'contact_substrate_total_A_per_um': '0'}
        result = audit_accepted_point(row)
        self.assertTrue(result['coupled_carrier_and_continuity_metrics_unmeasured'])
        self.assertIsNone(result['checks']['global_continuity'])
        self.assertFalse(result['pass'])

    def test_checked_coupled_metrics_still_must_pass_actual_physical_gates(self):
        proof = {k: True for k in ('coupled_acceptance_checked', 'physical_residuals_finite',
            'actual_load_line_satisfied', 'carrier_row_enforced', 'carrier_row_satisfied',
            'global_continuity_enforced', 'global_continuity_satisfied', 'terminal_kcl_satisfied')}
        row = {'converged': '1', 'solver_method': 'coupled_load_line_newton',
               'validation_diagnostics': json.dumps(proof),
               'inner_voltage_V': '4e-6', 'outer_voltage_V': '1e-5', 'current_total_A_per_um': '6e-13',
               'carrier_row_max_ratio': '.0009', 'global_electron_continuity_closure_ratio': '.001',
               'global_hole_continuity_closure_ratio': '.001',
               'contact_drain_total_A_per_um': '6e-13', 'contact_source_total_A_per_um': '-6e-13',
               'contact_gate_total_A_per_um': '0', 'contact_substrate_total_A_per_um': '0'}
        self.assertTrue(audit_accepted_point(row)['pass'])
        row['global_hole_continuity_closure_ratio'] = '.1'
        self.assertFalse(audit_accepted_point(row)['pass'])

    def test_saved_state_gate_ignores_claimed_pass_and_checks_actual_currents(self):
        totals = {'drain': 6e-13, 'source': -6e-13, 'gate': 0, 'substrate': 0}
        terminal = {'inner_bias_V': 4e-6, 'state_sha256': 'test', 'kcl_pass': True,
                    'totals_A_per_um': totals,
                    'components': {c: {'contact_current_extractor_electron_A_per_um': i,
                                       'contact_current_extractor_hole_A_per_um': 0} for c, i in totals.items()}}
        self.assertTrue(audit_saved_state(terminal, 1e-5)['pass'])
        totals['source'] = 6e-13
        self.assertFalse(audit_saved_state(terminal, 1e-5)['checks']['terminal_kcl'])

    def test_raw_hole_flux_sign_is_preserved_in_component_identity(self):
        terminal = {'inner_bias_V': 0, 'state_sha256': 'test',
                    'totals_A_per_um': {c: 0 for c in ('drain', 'source', 'gate', 'substrate')},
                    'components': {c: {'contact_current_extractor_electron_A_per_um': 2e-13,
                                       'contact_current_extractor_hole_A_per_um': 2e-13}
                                   for c in ('drain', 'source', 'gate', 'substrate')}}
        self.assertTrue(audit_saved_state(terminal, 0)['pass'])
        terminal['components']['drain']['contact_current_extractor_hole_A_per_um'] = -2e-13
        self.assertFalse(audit_saved_state(terminal, 0)['checks']['component_identity'])

if __name__ == '__main__':
    unittest.main()
