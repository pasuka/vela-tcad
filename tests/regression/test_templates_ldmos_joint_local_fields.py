"""Physical-field gate weighting and resolved-density mask checks."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from audit_templates_ldmos_joint_local_fields import carrier_metric
from evidence_paths import candidate_path


class LocalFieldGateTests(unittest.TestCase):
    def test_copied_evidence_mapping_preserves_root_boundary(self):
        self.assertEqual(candidate_path('/remote/run/results/a.json', ('/remote/run', 'copy')),
                         Path('copy/results/a.json'))
        self.assertEqual(candidate_path('local/a.json'), Path('local/a.json'))
        for value in ('/remote/runner/a.json', '/remote/run/../outside.json'):
            with self.assertRaises(ValueError):
                candidate_path(value, ('/remote/run', 'copy'))

    def test_volume_weighted_relative_not_absolute_density(self):
        result = carrier_metric([1., 100.], [1.1, 120.], [3., 1.], [True, True], 0.)
        self.assertAlmostEqual(result['weighted_relative_rms']**2, (3*.1**2+.2**2)/4)

    def test_floor_is_global_peak_and_strict_for_interface(self):
        result = carrier_metric([100., 1., 2.], [100., 9., 2.02],
                                [1., 1., 1.], [False, True, True], .01)
        self.assertEqual(result['resolved_nodes'], 1)
        self.assertAlmostEqual(result['weighted_relative_rms'], .01)

    def test_zero_volume_does_not_change_metric(self):
        result = carrier_metric([1., 1.], [2., 1.], [0., 3.], [True, True], 0.)
        self.assertEqual(result['weighted_relative_rms'], 0.)

    def test_no_resolved_interface_is_error(self):
        with self.assertRaises(ValueError):
            carrier_metric([1.], [1.], [1.], [False], 0.)

    def test_invalid_inputs_are_errors(self):
        for reference, candidate, areas, floor in (
            ([1.], [float('nan')], [1.], 0.),
            ([0.], [1.], [1.], 0.),
            ([1.], [-1.], [1.], 0.),
            ([1.], [1.], [-1.], 0.),
            ([1.], [1., 2.], [1.], 0.),
            ([1.], [1.], [1.], 1.),
        ):
            with self.subTest(reference=reference, candidate=candidate, areas=areas, floor=floor):
                with self.assertRaises(ValueError):
                    carrier_metric(reference, candidate, areas, [True], floor)


if __name__ == '__main__':
    unittest.main()
