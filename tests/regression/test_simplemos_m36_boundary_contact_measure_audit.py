import hashlib
import json
from pathlib import Path
import unittest

from tests.regression.simplemos_evidence_chain import (
    assert_source_hashes_current_or_m44,
)


REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/boundary_contact_measure_audit"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m36_boundary_contact_measure_audit_evidence.json")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SimpleMosM36BoundaryContactMeasureAuditTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = json.loads((
            ROOT / "m36_boundary_contact_measure_audit_report.json").read_text())

    def test_complete_mesh_and_classes_are_accepted(self) -> None:
        self.assertTrue(self.report["acceptance"]["all_checks_pass"])
        self.assertEqual(self.report["mesh"]["node_count"], 1482)
        self.assertEqual(self.report["mesh"]["triangle_count"], 2746)
        classes = {row["node_class"] for row in self.report["node_class_summary"]}
        self.assertEqual(classes, {
            "interior", "material_interface", "external_boundary",
            "direct_contact", "contact_first_layer",
        })

    def test_assembled_si_measures_match_sentaurus(self) -> None:
        self.assertLess(self.report["maximum_si_node_relative_error"], 5e-14)
        supports = {(row["contact"], row["support_class"]): row
                    for row in self.report["contact_support_summary"]}
        for contact in ("source", "drain", "substrate"):
            for support_class in ("direct_contact", "first_layer"):
                self.assertLess(
                    supports[(contact, support_class)]["relative_error_maximum"],
                    5e-14)

    def test_local_redistribution_does_not_hide_nodal_closure(self) -> None:
        self.assertGreater(self.report["maximum_local_relative_error"], 1.0)
        self.assertLess(self.report["maximum_si_node_relative_error"], 5e-14)

    def test_default_model_is_unchanged(self) -> None:
        self.assertFalse(self.report["execution"]["default_model_changed"])

    def test_frozen_artifact_and_source_hashes_match(self) -> None:
        evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
        self.assertEqual(evidence["status"], "frozen")
        for artifact in evidence["artifacts"]:
            self.assertEqual(artifact["sha256"], sha256(REPO / artifact["path"]))
        assert_source_hashes_current_or_m44(self, evidence)


if __name__ == "__main__":
    unittest.main()
