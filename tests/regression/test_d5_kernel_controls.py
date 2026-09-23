import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_d5_kernel_controls import configure_candidate
from analyze_templates_ldmos_d5_kernel_controls import analyze


class KernelControlsTest(unittest.TestCase):
    def test_analysis_refuses_incomplete_missing_and_duplicate_evidence(self):
        with self.assertRaises(ValueError):
            analyze(dict(status='running'))
        with self.assertRaises(ValueError):
            analyze(dict(status='completed',cases=[],planned_curves=4))
        duplicate=dict(repeat=0,mode='structure',gate=4)
        with self.assertRaises(ValueError):
            analyze(dict(status='completed',cases=[duplicate,duplicate],planned_curves=2))

    def test_fermi_and_combined_are_explicit_separate_candidates(self):
        base=dict(solver=dict(temperature=300))
        f=configure_candidate(base,'fermi');both=configure_candidate(base,'combined')
        self.assertTrue(f['solver'].pop('diagnostic_fermi_node_cache'))
        self.assertFalse(f['solver'].pop('diagnostic_reuse_jacobian_structure'))
        self.assertTrue(both['solver'].pop('diagnostic_fermi_node_cache'))
        self.assertTrue(both['solver'].pop('diagnostic_reuse_jacobian_structure'))
        self.assertEqual(f,base);self.assertEqual(both,base)

    def test_structure_candidate_changes_only_explicit_cache_switch(self):
        base=dict(solver=dict(temperature=300,block_absolute_convergence=dict(psi=5e-8)),
                  sweep=dict(step=.2),mesh_file='frozen_mesh')
        off=configure_candidate(base,'baseline');on=configure_candidate(base,'structure')
        self.assertFalse(off['solver'].pop('diagnostic_reuse_jacobian_structure'))
        self.assertTrue(on['solver'].pop('diagnostic_reuse_jacobian_structure'))
        self.assertEqual(off,base);self.assertEqual(on,base)
        on['solver']['block_absolute_convergence']['psi']=0
        self.assertEqual(base['solver']['block_absolute_convergence']['psi'],5e-8)

    def test_unknown_candidate_rejected(self):
        with self.assertRaises(ValueError):configure_candidate({'solver':{}},'typo')


if __name__=='__main__':unittest.main()
