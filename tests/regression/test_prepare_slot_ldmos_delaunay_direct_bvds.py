from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "prepare_slot_ldmos_delaunay_direct_bvds.py"
SPEC = importlib.util.spec_from_file_location("delaunay_bvds", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class DelaunayDirectBvdsPreparationTest(unittest.TestCase):
    def test_cases_keep_atomic_profile_and_do_not_claim_true_ialmob(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root / "seed.csv"
            seed.write_text("node_id,psi\n", encoding="utf-8")
            template = {
                "mesh_geometry": dict(MODULE.ATOMIC_MESH),
                "solver": {
                    "mobility": {"model": "masetti_field_lombardi"},
                    "impact_ionization": {
                        **MODULE.ATOMIC_IMPACT,
                        "source_jacobian": "local_ad",
                    },
                },
                "contacts": [
                    {"name": "source", "bias": 0.0},
                    {"name": "drain", "bias": 0.0},
                ],
                "sweep": {
                    "boundary_control": {},
                    "external_circuit": {},
                },
            }
            (root / "simulation_05_avalanche_on_60v.json").write_text(
                json.dumps(template), encoding="utf-8")
            manifest = MODULE.prepare(
                root, off_seed=seed, off_inner_v=15.0, off_outer_v=1000.0,
                proxy_seed=seed, proxy_inner_v=1.0, proxy_outer_v=60.0)
            self.assertFalse(manifest["true_sentaurus_ialmob_implemented"])
            for case in ("masetti_field", "enhanced_lombardi_proxy"):
                homotopy_names = manifest["configs"][case]["homotopy"]
                self.assertEqual(len(homotopy_names), len(MODULE.HOMOTOPY_SCALES) + 1)
                homotopy = [json.loads((root / name).read_text()) for name in homotopy_names]
                self.assertEqual(homotopy[0]["solver"]["method"], "poisson_only")
                self.assertEqual(
                    homotopy[0]["solver"]["poisson_only_quasi_fermi_initialization"],
                    "preserve_initial")
                self.assertEqual(homotopy[0]["solver"]["max_update"], 0.1)
                self.assertNotIn("impact_ionization", homotopy[1]["solver"])
                self.assertEqual(
                    homotopy[-1]["solver"]["impact_ionization"]["source_geometry_scale"],
                    1.0)
                self.assertTrue(all(
                    not deck["sweep"]["external_circuit"]["enabled"]
                    for deck in homotopy))
                block_probe_names = manifest["configs"][case]["poisson_block_probes"]
                self.assertEqual(len(block_probe_names), len(MODULE.POISSON_CAP_PROBES) + 1)
                block_probes = [json.loads((root / name).read_text()) for name in block_probe_names]
                self.assertTrue(all(
                    deck["simulation_type"] == "newton_block_step_probe"
                    and deck["block_modes"] == ["poisson_only"]
                    for deck in block_probes))
                self.assertEqual(block_probes[0]["solver"]["max_update"], 0.0)
                self.assertEqual(block_probes[-1]["solver"]["max_update"], 5.0)
                deck = json.loads((root / manifest["configs"][case]["probe"]).read_text())
                self.assertEqual(deck["mesh_geometry"], MODULE.ATOMIC_MESH)
                for key, value in MODULE.ATOMIC_IMPACT.items():
                    self.assertEqual(deck["solver"]["impact_ionization"][key], value)
                self.assertEqual(
                    deck["sweep"]["external_circuit"]["coupled_linear_solver"],
                    "direct_bordered")
                self.assertEqual(deck["solver"]["abstol"], 1.0e-9)
                self.assertEqual(
                    deck["solver"]["stall_residual_floor"], 1.0e-8)
                self.assertEqual(
                    deck["solver"]["quasi_fermi_update_limit_mode"],
                    "componentwise_poisson_recorrect")
                closure = deck["solver"]["global_continuity_closure"]
                self.assertEqual(closure["mode"], "enforce")
                self.assertEqual(closure["tolerance"], 1.0e-2)
                self.assertEqual(closure["source_floor"], 1.0e-11)
                self.assertEqual(closure["absolute_tolerance"], 1.0e-12)
                self.assertTrue(
                    deck["solver"]["carrier_row_qualified_stall_acceptance"])
                self.assertEqual(
                    deck["solver"]["carrier_row_convergence"]["mode"], "enforce")
                self.assertEqual(
                    deck["solver"]["carrier_row_convergence"]["eps_row"], 1.0e-3)
                self.assertEqual(
                    deck["solver"]["carrier_row_convergence"]["min_source_scale"],
                    1.0e-13)
                self.assertEqual(
                    closure["linear_system_stabilization"],
                    "replace_free_row_sum")
                self.assertFalse(
                    deck["_surface_mobility_ablation"]["true_sentaurus_ialmob_implemented"])
            self.assertEqual(
                manifest["global_continuity_contract"]
                ["linear_system_stabilization"],
                "replace_free_row_sum")


if __name__ == "__main__":
    unittest.main()
