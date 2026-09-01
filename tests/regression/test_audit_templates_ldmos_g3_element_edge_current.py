import csv
import tempfile
import unittest
from pathlib import Path

from scripts.audit_templates_ldmos_g3_element_edge_current import (
    capability_result,
    compare_field_directories,
    drain_current,
)


def write_field(root: Path, name: str, values: list[tuple[int, float]]) -> None:
    path = root / "fields" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["node_id", "component0"])
        writer.writeheader()
        for node, value in values:
            writer.writerow({"node_id": node, "component0": value})


class ElementEdgeCurrentAuditTest(unittest.TestCase):
    def test_field_directory_identity_is_byte_strict(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            baseline = root / "baseline"
            candidate = root / "candidate"
            write_field(baseline, "eContinuityRhs_region0.csv", [(1, 2.0)])
            write_field(candidate, "eContinuityRhs_region0.csv", [(1, 2.0)])
            result = compare_field_directories(baseline, candidate)
            self.assertTrue(result["all_common_fields_byte_identical"])
            write_field(candidate, "eContinuityRhs_region0.csv", [(1, 3.0)])
            result = compare_field_directories(baseline, candidate)
            self.assertEqual(
                result["byte_different_common_fields"],
                ["eContinuityRhs_region0.csv"],
            )

    def test_capability_probe_requires_native_edge_vector_error(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            root = Path(root_text)
            log = root / "probe.log"
            tcl = root / "probe.tcl"
            log.write_text(
                "writing probe_high_vsv_node3721_minus_edge_api_newton_1_0.tdr\n"
                "Tried to read undefined Edge-Vector eCurrentDensity !\n",
                encoding="utf-8",
            )
            tcl.write_text(
                "set c [$data ReadCoefficient]\n"
                "set j [$data ReadVector $::des_data_edge \"eCurrentDensity\"]\n",
                encoding="utf-8",
            )
            result = capability_result(log, tcl)
            self.assertTrue(result["native_error_observed"])
            self.assertTrue(result["coefficient_positive_control_precedes_target_call"])
            self.assertTrue(result["newton_iteration_completed_before_capability_error"])
            self.assertFalse(result["public_edge_current_dataset_available"])

    def test_drain_current_reads_final_table(self) -> None:
        with tempfile.TemporaryDirectory() as root_text:
            path = Path(root_text) / "run.log"
            path.write_text(
                " drain 1.000E-01 1.626E-06 -2.490E-14 1.626E-06\n",
                encoding="utf-8",
            )
            result = drain_current(path)
            self.assertEqual(result["electron_A"], 1.626e-6)
            self.assertEqual(result["voltage_V"], 0.1)


if __name__ == "__main__":
    unittest.main()
