"""Cross-point study must preserve an explicit control and a fixed backend."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from run_templates_ldmos_analysis_reuse import validate_options


class StudyOptions(unittest.TestCase):
    def test_all_four_backends_share_worker_control(self):
        for backend in ('sparselu','umfpack','mumps','superlu_mt'):
            self.assertEqual(validate_options(backend,1,['worker','reuse']),'worker')
    def test_old_three_mode_control_is_preserved(self):
        self.assertEqual(validate_options('strumpack',2,['subprocess','worker','reuse']),'subprocess')
    def test_reject_ambiguous_or_incomparable_controls(self):
        for modes in (['reuse'],['worker','worker','reuse'],['subprocess','worker']):
            with self.assertRaises(ValueError):validate_options('strumpack',2,modes)
        for backend in ('sparselu','umfpack'):
            with self.assertRaises(ValueError):validate_options(backend,2,['worker','reuse'])


if __name__=='__main__':unittest.main()
