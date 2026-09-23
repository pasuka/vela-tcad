from __future__ import annotations

import copy
import unittest

from scripts.prepare_slot_ldmos_stage06_restart import (
    RestartPreparationError,
    prepare_restart_document,
)


BASE = {
    "sweep": {
        "bias_points": [60.0, 1000.0],
        "initial_state_file": "outputs/stages/05_avalanche_on_60v/final_state.h5",
        "external_circuit": {
            "mode": "series_resistor",
            "initial_inner_voltage_V": 38.5209,
            "max_inner_voltage_step_V": 0.01,
        },
        "boundary_control": {"resume": True},
    }
}


class SlotLdmosStage06RestartTest(unittest.TestCase):
    def test_injects_stage05_inner_voltage_and_reach_step(self) -> None:
        row = {
            "converged": "1",
            "outer_voltage_V": "60.0",
            "inner_voltage_V": "0.73339774542603375",
        }
        result = prepare_restart_document(copy.deepcopy(BASE), row, 0.25)
        circuit = result["sweep"]["external_circuit"]
        self.assertEqual(circuit["initial_inner_voltage_V"], 0.73339774542603375)
        self.assertEqual(circuit["max_inner_voltage_step_V"], 0.25)
        self.assertFalse(result["sweep"]["boundary_control"]["resume"])

    def test_rejects_mismatched_stage_boundary(self) -> None:
        row = {
            "converged": "1",
            "outer_voltage_V": "50.0",
            "inner_voltage_V": "0.6",
        }
        with self.assertRaisesRegex(
            RestartPreparationError, "does not match the first Stage 06 point"
        ):
            prepare_restart_document(copy.deepcopy(BASE), row, 0.25)

    def test_resume_is_explicit(self) -> None:
        row = {
            "converged": "1",
            "outer_voltage_V": "60.0",
            "inner_voltage_V": "0.73339774542603375",
        }
        result = prepare_restart_document(copy.deepcopy(BASE), row, 1.0, True)
        self.assertTrue(result["sweep"]["boundary_control"]["resume"])
        self.assertTrue(result["_stage06_restart"]["resume_boundary_control"])


if __name__ == "__main__":
    unittest.main()
