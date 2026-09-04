import hashlib
import json
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
DESCRIPTOR = CASE / "genius_bjt_sentaurus2022_reference.json"


class GeniusBjtReferenceFixtureTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = json.loads(DESCRIPTOR.read_text(encoding="utf-8"))

    def test_descriptor_identifies_reusable_checked_in_case(self):
        self.assertEqual(
            self.reference["schema"], "vela.reference_tcad.checked_in.v1"
        )
        self.assertEqual(self.reference["case"], "genius_bjt_sentaurus2022")
        self.assertEqual(self.reference["device"], "two_dimensional_npn_bjt")

    def test_checked_in_case_paths_exist(self):
        relative_paths = []
        relative_paths.extend(self.reference["contracts"])
        relative_paths.extend(self.reference["reports"])
        relative_paths.extend(self.reference["figures"])
        for simulation in self.reference["simulations"]:
            for key in (
                "sdevice_cmd",
                "vela_model_relaxation",
                "vela_config",
                "reference_curve",
            ):
                if key in simulation:
                    relative_paths.append(simulation[key])

        for path in relative_paths:
            with self.subTest(path=path):
                self.assertTrue((CASE / path).is_file())

        for path in self.reference["reproduction_scripts"]:
            with self.subTest(path=path):
                self.assertTrue((REPO / path).is_file())

    def test_recorded_hashes_match_checked_in_inputs(self):
        for path, expected in self.reference["sha256"].items():
            with self.subTest(path=path):
                actual = hashlib.sha256((CASE / path).read_bytes()).hexdigest()
                self.assertEqual(actual, expected)

    def test_acceptance_boundary_is_explicit(self):
        status = self.reference["validation_status"]
        self.assertTrue(status["m1_terminal_31_point"])
        self.assertTrue(status["coarse_spatial_state_vce3"])
        self.assertFalse(status["coarse_default_current_density"])
        self.assertTrue(status["local_refined_vce3_overall"])
        self.assertFalse(status["cell_first_promoted_to_global_default"])


if __name__ == "__main__":
    unittest.main()
