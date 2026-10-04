"""Validate higher-order derivative and physical-state perturbation accounting."""
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
import run_simplemos_m79c_current_difference_calibration as m


class DifferenceCalibration(unittest.TestCase):
    def test_richardson_cancels_cubic_current_truncation(self):
        f=lambda x: 3*x+7*x**3
        central=lambda h: (f(h)-f(-h))/(2*h)
        self.assertAlmostEqual(m.richardson(central(.08),central(.04)),3,places=14)

    def test_physical_and_reference_state_and_boltzmann_are_consistent(self):
        row={"psi":.4,"phin":.05,"phip":0,"electron_qf_increment_V":.001,
             "electron_qf_reference_V":.049,"holes_m3":5e18,"electrons_m3":1e17}
        result=m.perturb([row],[[.5,.2]],.004,.025)[0]
        self.assertEqual(row["psi"],.4)
        self.assertAlmostEqual(result["phin"],result["electron_qf_reference_V"]+result["electron_qf_increment_V"])
        self.assertAlmostEqual(math.log(result["electrons_m3"]/row["electrons_m3"]),(.002-.0008)/.025)
        self.assertAlmostEqual(math.log(result["holes_m3"]/row["holes_m3"]),-.002/.025)
        self.assertEqual(result["phip"],row["phip"])


if __name__ == "__main__":unittest.main()
