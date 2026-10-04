import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_linear_closeout import policy_input, trajectory


class CloseoutProtocolTest(unittest.TestCase):
    def test_default_removes_legacy_overrides_without_changing_physics(self):
        source = {'electrothermal_linear_solver': 'sparselu_colamd',
                  'reuse_linear_analysis': False, 'reuse_sparselu_symbolic': False,
                  'diagnostic_ialmob_screening_method': 'legacy',
                  'mobility_SI': {'ialmob': {'screening_method':'legacy','delta':2.4}},
                  'sweep': {'maximum_step_V': 1.3333333333333333},
                  'reuse_jacobian_structure': True, 'temperature_K':[300.,311.]}
        original = copy.deepcopy(source)
        default = policy_input(source, True)
        control = policy_input(source, False)
        self.assertEqual(source, original)
        self.assertNotIn('electrothermal_linear_solver', default)
        self.assertNotIn('reuse_sparselu_symbolic', default)
        self.assertNotIn('reuse_linear_analysis', default)
        self.assertNotIn('screening_method', default['mobility_SI']['ialmob'])
        self.assertEqual(default['mobility_SI']['ialmob']['delta'],2.4)
        self.assertEqual(default['temperature_K'], source['temperature_K'])
        self.assertEqual(default['sweep'], source['sweep'])
        self.assertEqual(control, dict(default, reuse_linear_analysis=False))

    def test_default_does_not_add_model_options(self):
        self.assertEqual(policy_input({'mesh_file':'mesh.json'},True),{'mesh_file':'mesh.json'})

    def test_trajectory_includes_failed_attempts(self):
        ledger = {'runs':[{'bias_V':1.,'gate':{'pass_gate':False}},
                          {'bias_V':.5,'gate':{'pass_gate':True}}]}
        self.assertEqual(trajectory(ledger),[(1.,False),(.5,True)])


if __name__=='__main__': unittest.main()
