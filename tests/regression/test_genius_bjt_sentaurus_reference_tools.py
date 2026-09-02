from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_genius_bjt_sentaurus_reference.py"
SPEC = importlib.util.spec_from_file_location("audit_genius_bjt", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class GeniusBjtReferenceToolsTest(unittest.TestCase):
    def test_outside_distance_is_zero_inside_and_distance_outside(self) -> None:
        self.assertEqual(MODULE.outside_distance(2.0, 1.25, 4.75), 0.0)
        self.assertEqual(MODULE.outside_distance(1.0, 1.25, 4.75), 0.25)
        self.assertEqual(MODULE.outside_distance(5.0, 1.25, 4.75), 0.25)

    def test_expected_doping_matches_flat_profile_peaks(self) -> None:
        donor, acceptor = MODULE.expected_doping(3.5, 0.0)
        self.assertLess(abs(donor - (7.0e19 + 5.0e15)) / donor, 1.0e-12)
        self.assertLess(abs(acceptor - 4.6e18) / acceptor, 1.0e-12)

    def test_analytical_profile_has_genius_one_e_fold_distance(self) -> None:
        value = MODULE.analytical_profile(
            1.13, 0.2, 6.0e17, (1.25, 0.0, 4.75, 0.35), (0.12, 0.16)
        )
        self.assertLess(abs(value / 6.0e17 - 1.0 / MODULE.math.e), 1.0e-12)


if __name__ == "__main__":
    unittest.main()
