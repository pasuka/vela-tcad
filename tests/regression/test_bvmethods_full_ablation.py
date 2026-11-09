"""Scientific-output qualification for the isolated full-Jacobian ablation."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[2] / 'scripts/diagnostics/bvmethods_full_ablation.py'
spec = importlib.util.spec_from_file_location('bvmethods_ablation', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class QualificationTests(unittest.TestCase):
    def test_signed_current_interpolation_does_not_extrapolate(self):
        rows = [{'current_total_A_per_um': str(i), 'inner_voltage_V': str(v)}
                for i, v in [(0., 6.), (-2e-4, 7.)]]
        result = module.interpolate(rows, 1e-4)
        self.assertEqual(result['inner_V'], 6.5)
        self.assertEqual(result['bracket_width_V'], 1.)
        self.assertIsNone(module.interpolate(rows, 1.3e-3))
        self.assertIsNone(module.interpolate(rows[:1], 1e-4))

    def test_small_update_cannot_replace_the_joint_device_gate(self):
        g = dict(physical_residuals_finite=True, terminal_kcl_satisfied=True,
                 carrier_row_enforced=False, carrier_row_satisfied=False,
                 global_continuity_enforced=False, global_continuity_satisfied=False,
                 actual_load_line_satisfied=True, density_update_mode='enforce_joint',
                 density_update_valid=True, density_update_rms=.1,
                 coupled_device_residual_inf=1.01157e-8, actual_load_line_residual_V=0.)
        row = {'outer_voltage_V': '.56', 'validation_diagnostics': json.dumps(g)}
        self.assertFalse(module.gate(row))
        g['coupled_device_residual_inf'] = 9e-9
        row['validation_diagnostics'] = json.dumps(g)
        self.assertTrue(module.gate(row))
        g['terminal_kcl_satisfied'] = False
        row['validation_diagnostics'] = json.dumps(g)
        self.assertFalse(module.gate(row))

    def test_atomic_status_roundtrip(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary)/'case/status.json'
            module.put(target, {'phase': 'P0', 'finished': False})
            self.assertEqual(module.read(target)['phase'], 'P0')
            self.assertFalse(target.with_suffix('.json.tmp').exists())


if __name__ == '__main__':
    unittest.main()
