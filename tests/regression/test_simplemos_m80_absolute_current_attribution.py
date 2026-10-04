"""Off-state current must remain separate from actual cross-code error."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
import run_simplemos_m80_absolute_current_attribution as m


class ErrorAccounting(unittest.TestCase):
    def test_large_mapped_current_cancellation_does_not_become_actual_error(self):
        r=m.exact_split(1.1e-8,3e-5,1e-8)
        self.assertAlmostEqual(r["actual"],1e-9,delta=1e-23)
        self.assertAlmostEqual(r["state"]+r["target"],r["actual"],delta=1e-20)
        self.assertGreater(abs(r["state"]),1000*abs(r["actual"]))
        self.assertLess(r["state"],0)
        self.assertGreater(r["target"],0)


if __name__=="__main__":unittest.main()
