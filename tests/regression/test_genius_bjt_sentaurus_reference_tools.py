from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_genius_bjt_sentaurus_reference.py"
SPEC = importlib.util.spec_from_file_location("audit_genius_bjt", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONVERTER_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "convert_tcad_export.py"
CONVERTER_SPEC = importlib.util.spec_from_file_location("convert_tcad_export", CONVERTER_SCRIPT)
assert CONVERTER_SPEC is not None and CONVERTER_SPEC.loader is not None
CONVERTER = importlib.util.module_from_spec(CONVERTER_SPEC)
CONVERTER_SPEC.loader.exec_module(CONVERTER)

COMPARISON_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "compare_genius_bjt_sentaurus_vela.py"
)
COMPARISON_SPEC = importlib.util.spec_from_file_location(
    "compare_genius_bjt_sentaurus_vela", COMPARISON_SCRIPT
)
assert COMPARISON_SPEC is not None and COMPARISON_SPEC.loader is not None
COMPARISON = importlib.util.module_from_spec(COMPARISON_SPEC)
COMPARISON_SPEC.loader.exec_module(COMPARISON)


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

    def test_bjt_converter_selects_collector_for_voltage_sweep(self) -> None:
        contacts = [
            {"name": "base"},
            {"name": "emitter"},
            {"name": "collector"},
        ]
        self.assertEqual(CONVERTER.choose_sweep_contact(contacts, "bjt2d"), "collector")

    def test_comparison_selects_three_direct_terminal_currents(self) -> None:
        rows = [
            {
                "bias_V": "0.30000000000000004",
                "contact": contact,
                "current_total_A_per_um": current,
                "converged": "1",
            }
            for contact, current in (
                ("collector", "2e-6"),
                ("base", "1e-8"),
                ("emitter", "-2.01e-6"),
            )
        ]
        terminals, error = COMPARISON.select_vela_terminals(rows, 0.3)
        self.assertLess(error, 1.0e-12)
        self.assertEqual(set(terminals), {"collector", "base", "emitter"})

    def test_comparison_rejects_reconstructed_or_missing_terminal(self) -> None:
        rows = [
            {"bias_V": "0.3", "contact": "collector"},
            {"bias_V": "0.3", "contact": "emitter"},
        ]
        with self.assertRaisesRegex(ValueError, "missing contacts"):
            COMPARISON.select_vela_terminals(rows, 0.3)


if __name__ == "__main__":
    unittest.main()
