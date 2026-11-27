#!/usr/bin/env python3
"""Regression coverage for the frozen TransportModels DD/DG contract."""

from __future__ import annotations

import csv
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import transportmodels_fixed_contract as fixed  # noqa: E402
import run_transportmodels_dd_contact_basin_regression as dd_basin  # noqa: E402
import run_transportmodels_dd_dg_continuous_baseline as continuous  # noqa: E402
from compare_sentaurus_tdr_tdx import iter_dat_datasets  # noqa: E402
from transportmodels_idvd_fixture import stage_transportmodels_inputs  # noqa: E402
from vela_sdevice.mesh_import import _load_dfise  # noqa: E402


class TransportModelsFixedContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.base = {
            "materials_file": "obsolete.json",
            "solver": {
                "method": "newton",
                "mobility": {"model": "obsolete"},
                "electron_quantum_potential": {"enabled": False},
            },
        }

    def test_tracked_material_payload_contains_corrected_sentaurus_values(self) -> None:
        payload = json.loads(fixed.materials_path().read_text(encoding="utf-8"))
        materials = {row["name"]: row for row in payload["materials"]}
        self.assertEqual(14638914958.767616, materials["Si"]["ni"])
        self.assertEqual(4.0727403846153845,
                         materials["Si"]["electron_affinity_eV"])
        self.assertEqual(3.6, materials["Si"]["electron_quantum_gamma"])

    def test_dfise_stage_preserves_active_dopants_and_vela_doping(self) -> None:
        fixture = REPO / "reference_tcad/transportmodels_sentaurus2022/ci_u2022_12"
        with (fixture.parent / "vela/doping.csv").open(encoding="utf-8", newline="") as stream:
            archived = sorted(csv.DictReader(stream), key=lambda row: int(row["node_id"]))
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory) / "staged"
            stage_transportmodels_inputs(fixture, work, "pp6_des")
            datasets = dict(iter_dat_datasets(work / "n1_msh.dat"))
            imported = _load_dfise(work / "n1_msh.grd")
        donors = [float(row["donors_cm3"]) for row in archived]
        acceptors = [float(row["acceptors_cm3"]) for row in archived]
        self.assertTrue(any(d > 0 and a > 0 for d, a in zip(donors, acceptors)))
        for name, values in (("ArsenicActiveConcentration", donors),
                             ("BoronActiveConcentration", acceptors)):
            field = datasets[name]
            self.assertEqual(len(values), field["count"])
            self.assertIn("R.Substrate", field["regions"])
            self.assertEqual(values, field["values"])
        self.assertEqual([d - a for d, a in zip(donors, acceptors)],
                         datasets["DopingConcentration"]["values"])
        self.assertEqual(donors, [row["donors_cm3"] for row in imported["doping"]])
        self.assertEqual(acceptors, [row["acceptors_cm3"] for row in imported["doping"]])

    def test_dd_and_dg_differ_only_by_frozen_quantum_section(self) -> None:
        dd = fixed.apply_contract(self.base, "dd")
        dg = fixed.apply_contract(self.base, "dg")
        dd_solver, dg_classical_solver = fixed.controlled_solver_delta(dd, dg)
        self.assertEqual(dd_solver, dg_classical_solver)
        self.assertNotIn("electron_quantum_potential", dd["solver"])
        self.assertEqual(
            fixed.load_contract()["dg_quantum_contract"],
            dg["solver"]["electron_quantum_potential"],
        )
        self.assertEqual(
            "sentaurus_default", dd["solver"]["srh_density_coupling"]
        )
        self.assertEqual([], fixed.validate_config(dd, "dd"))
        self.assertEqual([], fixed.validate_config(dg, "dg"))

    def test_cold_start_retains_larger_qf_limit(self) -> None:
        cold = json.loads(json.dumps(self.base))
        cold["sweep"] = {"initialization": {"mode": "poisson_block"}}
        dd = fixed.apply_contract(cold, "dd")
        self.assertEqual(0.1, dd["solver"]["quasi_fermi_update_limit_V"])
        self.assertEqual([], fixed.validate_config(dd, "dd"))

    def test_audit_rejects_quantum_or_srh_contract_drift(self) -> None:
        dg = fixed.apply_contract(self.base, "dg")
        dg["solver"]["electron_quantum_potential"]["gamma"] = 1.0
        dg["solver"]["srh_density_coupling"] = "self_consistent"
        violations = fixed.validate_config(dg, "dg")
        self.assertTrue(any("electron_quantum_potential" in row
                            for row in violations))
        self.assertTrue(any("srh_density_coupling" in row
                            for row in violations))

    def test_contract_removes_and_rejects_unlisted_physics_models(self) -> None:
        raw = json.loads(json.dumps(self.base))
        raw["solver"]["impact_ionization"] = {"model": "van_overstraeten"}
        raw["solver"]["band_to_band"] = {"model": "nonlocal"}
        dd = fixed.apply_contract(raw, "dd")
        self.assertNotIn("impact_ionization", dd["solver"])
        self.assertNotIn("band_to_band", dd["solver"])
        self.assertEqual([], fixed.validate_config(dd, "dd"))

        dd["solver"]["impact_ionization"] = {"model": "van_overstraeten"}
        violations = fixed.validate_config(dd, "dd")
        self.assertTrue(any("solver.impact_ionization" in row
                            for row in violations))

    def test_bias_lattices_are_fixed_at_21_points(self) -> None:
        bias = fixed.load_contract()["bias_contract"]
        self.assertEqual(21, len(bias["idvg"]["gate_bias_V"]))
        self.assertEqual([-1.0, -0.84, -0.68],
                         bias["idvg"]["gate_bias_V"][:3])
        self.assertEqual(2.2, bias["idvg"]["gate_bias_V"][-1])
        self.assertEqual(21, len(bias["idvd"]["drain_bias_V"]))
        self.assertEqual(2.0, bias["idvd"]["drain_bias_V"][-1])

    def test_dd_contact_basin_overlay_is_fixed_and_fail_closed(self) -> None:
        dd = fixed.apply_dd_contact_basin_contract(self.base)
        self.assertEqual("contact_basin", dd["solver"]["quasi_fermi_reference"])
        self.assertEqual([], fixed.validate_dd_contact_basin_config(dd))
        self.assertEqual(
            "transportmodels-dd-contact-basin-v1",
            dd["fixed_dd_numerical_contract"]["id"],
        )

        dd["solver"]["quasi_fermi_reference"] = "none"
        violations = fixed.validate_dd_contact_basin_config(dd)
        self.assertTrue(any("quasi_fermi_reference" in row for row in violations))

    def test_contact_basin_internal_biases_are_strictly_increasing(self) -> None:
        exact = fixed.load_contract()["bias_contract"]["idvg"]["gate_bias_V"]
        execution = dd_basin.exact_and_bridge_biases(exact)
        self.assertEqual(sorted(execution), execution)
        self.assertEqual(len(execution), len(set(execution)))
        self.assertEqual(-0.5175, execution[execution.index(-0.52) + 1])
        self.assertEqual(-0.2, execution[execution.index(-0.2025) + 1])

    def test_completed_contact_basin_resume_has_no_sweep_bounds(self) -> None:
        self.assertIsNone(dd_basin.execution_bounds([]))
        self.assertEqual((-1.0, 2.2), dd_basin.execution_bounds([-1.0, 2.2]))

    def test_continuous_overlay_enables_contact_basin_for_both_branches(self) -> None:
        overlay = continuous.load_overlay()
        self.assertEqual(
            "contact_basin", overlay["solver_numerics"]["quasi_fermi_reference"]
        )
        self.assertTrue(overlay["rules"]["previous_accepted_state_only"])
        self.assertTrue(overlay["rules"]["pointwise_reclosure_forbidden"])
        self.assertEqual(
            10.0,
            overlay["deep_off_kcl_acceptance"]["min_current_to_kcl_ratio"],
        )
        self.assertEqual(
            "enforce",
            overlay["deep_off_kcl_acceptance"]["global_continuity_closure"]["mode"],
        )
        self.assertEqual(
            "direct_equilibrium", overlay["idvd_initialization"]["dd"]["mode"]
        )
        self.assertEqual(
            -1.0, overlay["idvd_initialization"]["dg"]["gate_start_V"]
        )
        self.assertEqual(
            1.0, overlay["idvd_initialization"]["dg"]["gate_stop_V"]
        )
        for branch in ("dd", "dg"):
            config = fixed.apply_contract(self.base, branch)
            fixed.deep_merge(config["solver"], overlay["solver_numerics"])
            self.assertEqual(
                "contact_basin", config["solver"]["quasi_fermi_reference"]
            )

    def test_continuous_nominal_lattices_match_fixed_contract(self) -> None:
        contract = fixed.load_contract()["bias_contract"]
        self.assertEqual(
            contract["idvg"]["gate_bias_V"],
            continuous.exact_curve_biases("dd", "idvg"),
        )
        self.assertEqual(
            contract["idvd"]["drain_bias_V"],
            continuous.exact_curve_biases("dg", "idvd"),
        )

    def test_continuous_curve_process_owns_seed_to_first_target_interval(self) -> None:
        idvg = continuous.exact_curve_biases("dg", "idvg")
        idvg_start, idvg_nominal = continuous.curve_execution_lattice(
            "dg", "idvg"
        )
        self.assertEqual(idvg[0], idvg_start)
        self.assertEqual(idvg, idvg_nominal)

        idvd = continuous.exact_curve_biases("dg", "idvd")
        idvd_start, idvd_nominal = continuous.curve_execution_lattice(
            "dg", "idvd"
        )
        self.assertEqual(idvd[1], idvd_start)
        self.assertEqual(idvd[1:], idvd_nominal)

    def test_contact_basin_failure_returns_nonzero(self) -> None:
        self.assertEqual(0, dd_basin.acceptance_exit_code({"overall_pass": True}))
        self.assertEqual(2, dd_basin.acceptance_exit_code({"overall_pass": False}))
        self.assertEqual(2, dd_basin.acceptance_exit_code({}))

    def test_json_serialization_replaces_nonfinite_numbers(self) -> None:
        source = {
            "finite": 1.0,
            "nested": [math.nan, math.inf, -math.inf],
        }
        payload, paths = fixed.strict_json_payload(source)
        self.assertEqual([None, None, None], payload["nested"])
        self.assertEqual(
            ["$.nested[0]", "$.nested[1]", "$.nested[2]"], paths
        )
        encoded = fixed.strict_json_text(source)
        self.assertNotIn("NaN", encoded)
        self.assertNotIn("Infinity", encoded)
        self.assertEqual(payload, json.loads(encoded))

    def test_artifact_root_must_be_explicit_and_existing(self) -> None:
        with self.assertRaises(ValueError):
            fixed.resolve_artifact_root(None, {})
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            self.assertEqual(root, fixed.resolve_artifact_root(root, {}))
            self.assertEqual(
                root,
                fixed.resolve_artifact_root(
                    None, {fixed.ARTIFACT_ROOT_ENV: str(root)}
                ),
            )
        with self.assertRaises(FileNotFoundError):
            fixed.resolve_artifact_root(REPO / "definitely-not-present", {})


if __name__ == "__main__":
    unittest.main()
