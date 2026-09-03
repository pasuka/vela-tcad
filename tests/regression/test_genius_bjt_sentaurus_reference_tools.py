from __future__ import annotations

import importlib.util
import tempfile
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

SPATIAL_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "compare_genius_bjt_spatial_fields.py"
)
SPATIAL_SPEC = importlib.util.spec_from_file_location(
    "compare_genius_bjt_spatial_fields", SPATIAL_SCRIPT
)
assert SPATIAL_SPEC is not None and SPATIAL_SPEC.loader is not None
SPATIAL = importlib.util.module_from_spec(SPATIAL_SPEC)
SPATIAL_SPEC.loader.exec_module(SPATIAL)

TRANSPORT_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "compare_genius_bjt_transport_fields.py"
)
TRANSPORT_SPEC = importlib.util.spec_from_file_location(
    "compare_genius_bjt_transport_fields", TRANSPORT_SCRIPT
)
assert TRANSPORT_SPEC is not None and TRANSPORT_SPEC.loader is not None
TRANSPORT = importlib.util.module_from_spec(TRANSPORT_SPEC)
TRANSPORT_SPEC.loader.exec_module(TRANSPORT)


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

    def test_comparison_expands_unique_accepted_terminal_row(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "terminals.csv"
            path.write_text(
                "index,vce_V,collector_A_per_um,base_A_per_um,emitter_A_per_um\n"
                "3,0.3,2e-6,1e-8,-2.01e-6\n",
                encoding="utf-8",
            )
            rows = COMPARISON.read_vela_terminal_rows(path)
        terminals, error = COMPARISON.select_vela_terminals(rows, 0.3)
        self.assertEqual(error, 0.0)
        self.assertEqual(set(terminals), {"collector", "base", "emitter"})

    def test_asserted_numerical_parity_checks_all_four_observables(self) -> None:
        rows = [
            {
                "Ic_absolute_log10_error": 0.01,
                "Ib_absolute_log10_error": 0.04,
                "Ie_absolute_log10_error": 0.02,
                "beta_absolute_log10_error": 0.03,
            },
            {
                "Ic_absolute_log10_error": 0.02,
                "Ib_absolute_log10_error": 0.06,
                "Ie_absolute_log10_error": 0.03,
                "beta_absolute_log10_error": 0.04,
            },
        ]
        contract = {
            "status": "asserted",
            "maximum_Ic_absolute_log10_error": 0.05,
            "maximum_Ib_absolute_log10_error": 0.05,
            "maximum_Ie_absolute_log10_error": 0.05,
            "maximum_beta_absolute_log10_error": 0.05,
            "reason": "test",
        }
        result = COMPARISON.evaluate_numerical_parity(rows, contract)
        self.assertFalse(result["pass"])
        self.assertEqual(
            result["observed_maximum_absolute_log10_error"]["Ib"], 0.06
        )

    def test_asserted_numerical_parity_fails_on_emitter_current(self) -> None:
        rows = [
            {
                "Ic_absolute_log10_error": 0.01,
                "Ib_absolute_log10_error": 0.01,
                "Ie_absolute_log10_error": 0.051,
                "beta_absolute_log10_error": 0.01,
            }
        ]
        contract = {
            "status": "asserted",
            "maximum_Ic_absolute_log10_error": 0.05,
            "maximum_Ib_absolute_log10_error": 0.05,
            "maximum_Ie_absolute_log10_error": 0.05,
            "maximum_beta_absolute_log10_error": 0.05,
            "reason": "test",
        }
        result = COMPARISON.evaluate_numerical_parity(rows, contract)
        self.assertFalse(result["pass"])
        self.assertEqual(
            result["observed_maximum_absolute_log10_error"]["Ie"], 0.051
        )

    def test_characterization_only_parity_is_not_asserted(self) -> None:
        result = COMPARISON.evaluate_numerical_parity(
            [], {"status": "characterization_only", "reason": "baseline"}
        )
        self.assertIsNone(result["pass"])

    def test_spatial_density_gate_uses_sentaurus_reference_floor(self) -> None:
        contract = {
            "units": "decade",
            "error": "log10_vela_minus_log10_sentaurus",
            "mask": {"type": "sentaurus_reference_min_cm3", "minimum_cm3": 1e10},
            "maximum_rmse": 0.05,
            "maximum_p95_absolute_error": 0.05,
            "maximum_absolute_error": 0.05,
        }
        result = SPATIAL.evaluate_gate([4.0, 0.04], [1.0, 1.0e12], contract)
        self.assertEqual(result["selected_node_count"], 1)
        self.assertEqual(result["excluded_node_count"], 1)
        self.assertTrue(result["pass"])

    def test_spatial_gate_checks_rmse_p95_and_maximum(self) -> None:
        contract = {
            "units": "V",
            "error": "vela_minus_sentaurus",
            "mask": {"type": "all_common_nodes"},
            "maximum_rmse": 0.1,
            "maximum_p95_absolute_error": 0.1,
            "maximum_absolute_error": 0.1,
        }
        result = SPATIAL.evaluate_gate([0.0, 0.11], [0.0, 0.0], contract)
        self.assertFalse(result["checks"]["maximum_absolute_error"])
        self.assertFalse(result["pass"])

    def test_transport_vector_metric_is_exact_for_equal_fields(self) -> None:
        result = TRANSPORT.vector_metrics(
            [(1.0, 0.0), (0.0, 2.0)],
            [(1.0, 0.0, 0.0), (0.0, 2.0, 0.0)],
            1e-6,
        )
        self.assertEqual(result["normalized_vector_rmse"], 0.0)
        self.assertAlmostEqual(result["global_vector_cosine_similarity"], 1.0)

    def test_transport_source_metric_preserves_equal_integrals(self) -> None:
        result = TRANSPORT.source_metrics(
            [1.0, -0.5],
            [1.0, -0.5],
            [2.0, 1.0],
            1e-6,
            [(0.0, 0.0), (1.0, 0.0)],
        )
        self.assertEqual(result["normalized_l1_error"], 0.0)
        self.assertEqual(result["absolute_shape_total_variation"], 0.0)
        self.assertEqual(result["peak_location"]["distance_um"], 0.0)
        self.assertAlmostEqual(
            result["signed_integral_A_per_um"]["ratio_vela_over_sentaurus"], 1.0
        )

    def test_overall_acceptance_requires_terminal_and_spatial_gates(self) -> None:
        self.assertTrue(COMPARISON.combine_acceptance(True, True, True))
        self.assertFalse(COMPARISON.combine_acceptance(True, True, True, False))
        self.assertFalse(COMPARISON.combine_acceptance(True, True, False))
        self.assertFalse(COMPARISON.combine_acceptance(True, False, True))
        self.assertFalse(COMPARISON.combine_acceptance(False, True, True))

    def test_spatial_summary_rejects_inconsistent_overall_flag(self) -> None:
        summary = {
            "fields": {
                "potential": {"gate": {"pass": True}},
                "holes": {"gate": {"pass": False}},
            },
            "overall_pass": True,
        }
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            COMPARISON.validated_spatial_state_pass(summary)

    def test_transport_summary_requires_every_asserted_metric(self) -> None:
        summary = {
            "current_density": {
                "electron": {"gate": {"pass": True}},
                "hole": {"gate": {"pass": False}},
            },
            "recombination": {
                "srh": {"gate": {"pass": True}},
                "auger": {"gate": {"pass": True}},
            },
            "overall_pass": False,
        }
        self.assertFalse(COMPARISON.validated_transport_source_pass(summary))

    def test_source_gate_checks_integral_and_shape_while_reporting_peak_location(self) -> None:
        metrics = TRANSPORT.source_metrics(
            [1.0, 0.5],
            [1.0, 0.5],
            [1.0, 1.0],
            1e-6,
            [(0.0, 0.0), (1.0, 0.0)],
        )
        contract = {
            "minimum_absolute_integral_ratio_vela_over_sentaurus": 0.99,
            "maximum_absolute_integral_ratio_vela_over_sentaurus": 1.01,
            "maximum_normalized_l1_error": 0.01,
            "maximum_absolute_shape_total_variation": 0.01,
        }
        result = TRANSPORT.evaluate_source_gate(metrics, contract)
        self.assertTrue(result["pass"])
        self.assertEqual(result["observed"]["peak_location_distance_um"], 0.0)
        self.assertIn("p95_absolute_log10_magnitude_error", result["observed"])
        self.assertNotIn("p95_absolute_log10_magnitude_error", result["checks"])
        self.assertNotIn("peak_location_distance_um", result["checks"])


if __name__ == "__main__":
    unittest.main()
