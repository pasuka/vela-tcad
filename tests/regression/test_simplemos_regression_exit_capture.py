import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
spec=importlib.util.spec_from_file_location('runner',Path(__file__).resolve().parents[2]/'scripts/run_simplemos_engineering_joined_regression_20260930.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class ExitCaptureTests(unittest.TestCase):
    def test_zero_is_recorded(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t);m.checked_step('zero',[sys.executable,'-c','raise SystemExit(0)'],p,p)
            self.assertEqual((p/'zero.exit').read_text(),'0\n')
    def test_nonzero_is_recorded_and_stops_pipeline(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)
            with self.assertRaisesRegex(RuntimeError,'exited 7'):m.checked_step('seven',[sys.executable,'-c','raise SystemExit(7)'],p,p)
            self.assertEqual((p/'seven.exit').read_text(),'7\n')

if __name__=='__main__':unittest.main()
