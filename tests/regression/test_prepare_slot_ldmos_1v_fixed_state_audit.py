import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "prepare_slot_fixed",
    REPO / "scripts" / "prepare_slot_ldmos_1v_fixed_state_audit.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PrepareSlotLdmosFixedStateAuditTest(unittest.TestCase):
    def test_builds_classical_ialmob_on_frozen_replay(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = {
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "materials_file": "materials.json",
                "output_csv": "old.csv",
                "solver": {
                    "method": "gummel_newton",
                    "handoff": {},
                    "impact_ionization": {"model": "van_overstraeten"},
                    "mobility": {"model": "masetti_field"},
                },
                "sweep": {
                    "external_circuit": {"enabled": True},
                    "boundary_control": {},
                    "breakdown": {},
                    "diagnostics": {},
                },
            }
            state = root / "state.csv"
            state.write_text("node_id,psi,phin,phip,electrons_m3,holes_m3\n")
            deck = MODULE.build_config(source, root, state, root / "audit")

            self.assertEqual(deck["solver"]["method"], "frozen_state")
            self.assertNotIn("impact_ionization", deck["solver"])
            self.assertEqual(
                deck["solver"]["mobility"]["model"],
                "masetti_field_lombardi")
            self.assertEqual(
                deck["solver"]["mobility"]["high_field_gradient_discretization"],
                "transport_cell_vector")
            self.assertEqual(
                deck["solver"]["srh_doping_dependence"]["electron"]
                    ["reference_doping_m3"],
                1.0e16)
            self.assertEqual(deck["solver"]["auger_cn_m6_per_s"], 2.90e-31)
            self.assertEqual(deck["solver"]["auger_cp_m6_per_s"], 1.028e-31)
            self.assertTrue(deck["sweep"]["frozen_state_compute_current"])
            self.assertNotIn("external_circuit", deck["sweep"])
            self.assertIn("terminal_balance", deck["sweep"]["diagnostics"])
            self.assertIn("srh_balance", deck["sweep"]["diagnostics"])
            self.assertEqual(
                deck["output_csv"],
                (root / "audit/outputs/iv.csv").resolve().as_posix())
            self.assertEqual(
                deck["sweep"]["diagnostics"]["terminal_balance"]["csv_file"],
                (root / "audit/outputs/terminal_balance.csv").resolve().as_posix())

    def test_supports_a_high_voltage_frozen_replay(self) -> None:
        source = {
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "solver": {"mobility": {}},
            "sweep": {},
        }
        deck = MODULE.build_config(
            source, Path("."), Path("state.csv"), bias_V=60.0)
        self.assertEqual(deck["sweep"]["bias_points"], [60.0])
        self.assertEqual(deck["sweep"]["start"], 60.0)
        self.assertEqual(deck["sweep"]["stop"], 60.0)

    def test_can_disable_enhanced_lombardi_for_an_ablation(self) -> None:
        source = {
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "solver": {"mobility": {}},
            "sweep": {},
        }
        deck = MODULE.build_config(
            source, Path("."), Path("state.csv"),
            mobility_model="masetti_field")
        self.assertEqual(deck["solver"]["mobility"]["model"], "masetti_field")
        self.assertNotIn("surface", deck["solver"]["mobility"])


if __name__ == "__main__":
    unittest.main()
