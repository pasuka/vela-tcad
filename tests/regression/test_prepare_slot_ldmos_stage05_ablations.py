from __future__ import annotations

import copy
import unittest

from scripts.prepare_slot_ldmos_stage05_ablations import (
    CASES,
    TARGET_INNER_V,
    build_case,
)


BASE = {
    "output_csv": "stage05/iv.csv",
    "solver": {
        "impact_ionization": {
            "coupling_mode": "self_consistent",
            "source_mapping_mode": "triangle_gss_gradqf_truncated",
        },
        "handoff": {"gummel_max_iter": 50},
    },
    "sweep": {
        "start": 10.0,
        "stop": 60.0,
        "bias_points": [10.0, 60.0],
        "initial_state_file": "outputs/stages/04/final_state.h5",
        "external_circuit": {"mode": "series_resistor"},
        "boundary_control": {"resume": True},
    },
}


class SlotLdmosStage05AblationTest(unittest.TestCase):
    def test_each_case_changes_only_declared_solver_axis(self) -> None:
        for name, settings in CASES.items():
            with self.subTest(name=name):
                document = build_case(copy.deepcopy(BASE), name)
                solver = document["solver"]
                impact = solver["impact_ionization"]
                self.assertEqual(impact["coupling_mode"], settings["coupling_mode"])
                self.assertEqual(
                    impact["source_jacobian"], settings["source_jacobian"]
                )
                self.assertEqual(
                    solver["handoff"]["gummel_max_iter"],
                    settings["gummel_max_iter"],
                )
                self.assertNotIn("external_circuit", document["sweep"])
                self.assertNotIn("boundary_control", document["sweep"])
                self.assertEqual(document["sweep"]["bias_points"], [TARGET_INNER_V])
                self.assertEqual(
                    document["sweep"]["initial_state_file"],
                    BASE["sweep"]["initial_state_file"],
                )


if __name__ == "__main__":
    unittest.main()
