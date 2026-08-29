from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = (ROOT / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m11_mobility_factorial_contract_v1.json")
SELF_REPORT = (ROOT / "build-release/reference_tcad/simplemos_sentaurus2022"
               / "m11_mobility_factorial/self_consistent_factorial_report.json")
FROZEN_REPORT = (ROOT / "build-release/reference_tcad/simplemos_sentaurus2022"
                 / "m11_mobility_factorial/frozen_vela/frozen_factorial_report.json")


def load_script():
    spec = importlib.util.spec_from_file_location(
        "simplemos_m11", ROOT / "scripts/run_simplemos_m11_mobility_factorial.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M11 = load_script()


class SimpleMosM11MobilityFactorialTest(unittest.TestCase):
    def test_contract_covers_complete_orthogonal_design(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        M11.validate_contract(contract)
        combinations = {tuple(bool(item[factor]) for factor in M11.FACTORS)
                        for item in contract["variants"]}
        self.assertEqual(8, len(combinations))
        self.assertFalse(contract["default_model_policy"]["modify_existing_default"])
        self.assertEqual(32, contract["bias_matrix"]["self_consistent_curve_count"])
        self.assertEqual(1632, contract["bias_matrix"]["direct_bias_points_per_solver"])

    def test_factorial_effect_recovers_main_and_interaction_terms(self) -> None:
        contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
        variants = {item["id"]: item for item in contract["variants"]}
        rows = []
        for variant in contract["variants"]:
            p = 1 if variant["phumob"] else -1
            e = 1 if variant["enormal"] else -1
            h = 1 if variant["hfs"] else -1
            rows.append({"variant": variant["id"],
                         "response": 2.0 * p - 3.0 * e + 4.0 * p * h})
        effects = M11.factorial_effects(rows, "response", variants)
        self.assertAlmostEqual(4.0, effects["phumob"])
        self.assertAlmostEqual(-6.0, effects["enormal"])
        self.assertAlmostEqual(8.0, effects["phumob:hfs"])
        self.assertAlmostEqual(0.0, effects["hfs"])

    def test_generated_reports_are_complete(self) -> None:
        if not SELF_REPORT.is_file() or not FROZEN_REPORT.is_file():
            self.skipTest("M11 execution evidence has not been generated")
        self_consistent = json.loads(SELF_REPORT.read_text(encoding="utf-8"))
        frozen = json.loads(FROZEN_REPORT.read_text(encoding="utf-8"))
        self.assertEqual(32, self_consistent["curve_count"])
        self.assertEqual(1632, self_consistent["direct_bias_points_per_solver"])
        self.assertEqual(16, frozen["state_count"])
        self.assertEqual(128, frozen["variant_evaluations"])
        self.assertFalse(frozen["sentaurus_variant_re_evaluation"])
        for report in (self_consistent, frozen):
            self.assertEqual(set(M11.TERMS), set(next(iter(
                report["aggregate_effects"].values())).keys()))


if __name__ == "__main__":
    unittest.main()
