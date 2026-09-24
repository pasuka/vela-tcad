"""Qualification and bounded load-line controller tests; no solver executable."""
from __future__ import annotations

import csv
from tests.regression import slot_state_fixture as h5
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[2] / "scripts/run_slot_ldmos_loadline_continuation.py"
M = None
if SCRIPT.exists():
    spec = importlib.util.spec_from_file_location("slot_loadline", SCRIPT)
    M = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(M)


def csv_write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def deck():
    return {"solver": {"method": "newton", "abstol": 1e-9,
        "carrier_row_convergence": {"mode": "enforce", "eps_row": .001, "scale_floor": 1e-30},
        "global_continuity_closure": {"mode": "enforce", "source_floor": 1e-11,
            "tolerance": .01, "absolute_tolerance": 1e-12}},
        "sweep": {"contact": "drain", "current_contact": "drain", "external_circuit": {
            "enabled": False, "mode": "series_resistor", "resistance_ohm_um": 1e12,
            "current_direction": 1, "residual_tolerance_V": 1e-6}},
        "contacts": [{"name": "drain", "bias": 0}]}


def iv_row(bias=1, current=1e-14):
    row = {"converged": "1", "bias_V": str(bias), "bias_contact": "drain",
           "current_contact": "drain", "current_total_A_per_um": str(current),
           "current_total": str(current * 1e6), "solver_method": "newton",
           "newton_convergence_reason": "residual", "newton_iterations": "1",
           "global_continuity_closure_satisfied": "1", "carrier_row_violations": "0",
           "final_psi_residual_norm": "1e-12", "final_electron_continuity_residual_norm": "1e-12",
           "final_hole_continuity_residual_norm": "1e-12"}
    for carrier in ("electron", "hole"):
        row.update({f"global_{carrier}_contact_flux": "2e-13",
                    f"global_{carrier}_integrated_source": "1e-13",
                    f"global_{carrier}_continuity_closure_mismatch": "1e-13",
                    f"global_{carrier}_continuity_closure_ratio": ".01",
                    f"global_{carrier}_continuity_closure_criterion": "absolute",
                    f"global_{carrier}_continuity_closure_qualified": "0"})
    return row


def term_rows():
    rows = []
    for node in range(3):
        row = {"node_id": str(node), "x": "0", "y": "0"}
        for carrier in ("electron", "hole"):
            row[f"{carrier}_continuity_active"] = "1" if node < 2 else "0"
            for name in ("flux", "recombination", "impact", "residual", "gauge", "boundary", "term_sum"):
                row[f"{carrier}_{name}"] = "0"
            row[f"{carrier}_flux_abs_sum"] = "1e-25"
        rows.append(row)
    return rows


class Base(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(M, "load-line controller implementation is missing")


class ScalarTests(Base):
    def test_signed_amp_per_um_units(self):
        self.assertLess(abs(M.residual(37.0327878157634, 1.7212063627360095e-14, 37.05, 1)), 1e-6)
        self.assertAlmostEqual(M.residual(1, -1e-12, 2, -1), 0)
        with self.assertRaises(M.EvidenceError):
            M.residual(1, math.nan, 2, 1)

    def test_backward_secant_and_same_sign_not_bracket(self):
        x, bracket = M.proposal([(1., 2.), (2., 4.)], 2., .5, 0)
        self.assertEqual(x, 1.5)
        self.assertFalse(bracket)

    def test_true_bracket_safeguard(self):
        x, bracket = M.proposal([(0., -1.), (1., 100.)], 0., .5, 0)
        self.assertTrue(bracket)
        self.assertGreater(x, 0)
        self.assertLessEqual(x, .5)

    def test_bracket_flag_never_claimed_for_step_outside_bracket(self):
        x, bracket = M.proposal([(0., -1.), (1., 1.), (100., .1)], 100., .25, 0)
        self.assertTrue(not bracket or 0 < x < 1)

    def test_singular_near_singular_and_nonfinite_slopes(self):
        for pairs in ([(0., 1.), (1., 1.)], [(0., 1.), (1., 1. + 1e-15)],
                      [(0., 1e308), (1., -1e308)]):
            x, _ = M.proposal(pairs, 1., .25, 0)
            self.assertTrue(math.isfinite(x))
            self.assertLessEqual(abs(x - 1), .25)

    def test_recovery_both_directions_and_minimum(self):
        self.assertEqual(M.halfway(2, 4, .01), 3)
        self.assertEqual(M.halfway(2, 0, .01), 1)
        self.assertIsNone(M.halfway(2, 2.0001, .01))


class QualificationTests(Base):
    def setUp(self):
        super().setUp()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.state = self.root / "state.h5"
        h5.write(self.state, [{"node_id": i, "psi": 0, "phin": 0, "phip": 0,
                              "electrons_m3": 1, "holes_m3": 1} for i in range(3)])
        self.iv = self.root / "iv.csv"
        self.terms = self.root / "terms.csv"
        self.row = iv_row()
        self.rows = term_rows()

    def qualify(self, **kwargs):
        csv_write(self.iv, [self.row]); csv_write(self.terms, self.rows)
        return M.qualify_files(self.state, self.iv, self.terms, 1., dict(deck(), mesh_file=str(self.root / "mesh.json")),
                               {0, 1, 2}, expected_count=2, updates=1, **kwargs)

    def test_absolute_criterion_and_all_active_not_source_rows(self):
        point = self.qualify()
        self.assertEqual(point["active_ids"], {"electron": [0, 1], "hole": [0, 1]})
        self.assertEqual(point["current_A_per_um"], 1e-14)
        self.assertEqual(point["state_sha256"], M.sha256(self.state))

    def test_zero_updates_only_at_same_imported_bias(self):
        csv_write(self.iv, [self.row]); csv_write(self.terms, self.rows)
        with self.assertRaisesRegex(M.NumericalFailure, "zero_coupled_update"):
            M.qualify_files(self.state, self.iv, self.terms, 1., dict(deck(), mesh_file=str(self.root / "mesh.json")), {0, 1, 2}, expected_count=2, updates=0)
        M.qualify_files(self.state, self.iv, self.terms, 1., dict(deck(), mesh_file=str(self.root / "mesh.json")), {0, 1, 2}, expected_count=2, updates=0, same_bias=True)

    def test_missing_duplicate_and_changed_active_rows(self):
        for rows in (term_rows()[:-1], term_rows() + [term_rows()[0]]):
            self.rows = rows
            with self.assertRaises(M.EvidenceError): self.qualify()
        self.rows = term_rows()
        self.rows[1]["hole_continuity_active"] = "0"
        self.rows[2]["hole_continuity_active"] = "1"
        with self.assertRaises(M.EvidenceError):
            self.qualify(expected_active={"electron": [0, 1], "hole": [0, 1]})

    def test_weak_source_row_cannot_hide_violation(self):
        self.rows[0]["hole_residual"] = "1e-26"
        with self.assertRaisesRegex(M.NumericalFailure, "local_audit"):
            self.qualify()

    def test_placeholder_wrong_bias_nonfinite_and_missing_current(self):
        for changes in ({"solver_method": "boundary_seed_state"}, {"bias_V": "2"},
                        {"current_total_A_per_um": "nan"}, {"current_total_A_per_um": ""}):
            self.row = {**iv_row(), **changes}
            with self.assertRaises(M.EvidenceError): self.qualify()

    def test_actual_closure_recomputed_even_when_flag_passes(self):
        self.row["global_hole_contact_flux"] = "1e-8"
        self.row["global_hole_continuity_closure_mismatch"] = "9.9999e-9"
        self.row["global_hole_continuity_closure_criterion"] = "relative"
        self.row["global_hole_continuity_closure_qualified"] = "1"
        self.row["global_hole_continuity_closure_ratio"] = ".99999"
        with self.assertRaisesRegex(M.NumericalFailure, "global_closure"):
            self.qualify()

    def test_tampered_state_or_evidence_rejected(self):
        point = self.qualify()
        M.verify_hashes(point["evidence"])
        self.state.write_text("changed")
        with self.assertRaises(M.EvidenceError): M.verify_hashes(point["evidence"])


class IdentityTests(Base):
    def test_relative_nested_geometry_resolved_and_hashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ("mesh.json", "material.json", "doping.csv", "geometry.json", "runner.exe"):
                (root / name).write_text("{}")
            base = deck()
            base.update(mesh_file="mesh.json", materials_file="material.json", node_doping_file="doping.csv")
            base["solver"]["mobility"] = {"ialmob": {"geometry_file": "geometry.json"}}
            resolved, identity = M.freeze_inputs(base, root, root / "runner.exe")
            self.assertEqual(resolved["solver"]["mobility"]["ialmob"]["geometry_file"], str((root / "geometry.json").resolve()))
            M.verify_hashes(identity["files"])
            (root / "geometry.json").write_text("changed")
            with self.assertRaises(M.EvidenceError): M.verify_hashes(identity["files"])

    def test_frozen_profile_rejects_relaxed_gates(self):
        for section, name, value in (("carrier_row_convergence", "eps_row", .1),
                                     ("global_continuity_closure", "mode", "report")):
            base = deck(); base["solver"][section][name] = value
            with self.assertRaises(ValueError): M.validate_profile(base)


class AnalyticBackend:
    """File-backed analytic device, exercising the real qualification adapter."""
    def __init__(self, function=lambda x: .02 * x / 1e12, fail=None):
        self.function = function
        self.fail = fail
        self.calls = []
        self.changed_confirmation = False
        self.interrupt = False

    def solve(self, bias, seed, directory, same_bias=False):
        self.calls.append((bias, seed["inner_V"]))
        if self.interrupt:
            self.interrupt = False
            raise KeyboardInterrupt()
        if self.fail:
            self.fail(bias, seed)
        state = directory / "state.h5"
        h5.write(state, [{"node_id": i, "psi": bias, "phin": 0, "phip": 0,
                          "electrons_m3": 1, "holes_m3": 1} for i in range(3)])
        iv = directory / "iv.csv"; terms = directory / "terms.csv"
        csv_write(iv, [iv_row(bias, self.function(bias))]); csv_write(terms, term_rows())
        return M.qualify_files(state, iv, terms, bias, dict(deck(), mesh_file=str(directory / "mesh.json")), {0, 1, 2},
                               expected_count=2, updates=0 if same_bias else 1, same_bias=same_bias)

    def confirm(self, point, outer, directory):
        if self.changed_confirmation:
            raise M.NumericalFailure("confirmation_state_changed: requires full requalification")
        path = directory / "confirmation.json"
        path.write_text(json.dumps({"inner_V": point["inner_V"], "outer_V": outer,
                                    "current_A_per_um": point["current_A_per_um"]}))
        return {"device_state_sha256": point["state_sha256"], "state_sha256": point["state_sha256"], "evidence": M.hashes([path])}

    def validate(self, point):
        M.verify_hashes(point["evidence"])
        checked = M.qualify_files(Path(point["state_file"]), Path(point["iv_file"]),
                                 Path(point["terms_file"]), point["inner_V"], dict(deck(), mesh_file=str(Path(point["state_file"]).parent / "mesh.json")), {0, 1, 2},
                                 expected_count=2, expected_active=point["active_ids"],
                                 updates=point["coupled_newton_updates"],
                                 same_bias=point["same_bias_requalification"])
        if checked != {k: v for k, v in point.items() if k != "id"}:
            raise M.EvidenceError("point metadata mismatch")


class CampaignTests(Base):
    def setUp(self):
        super().setUp()
        self.assertTrue(hasattr(M, "Campaign"), "persisted bounded controller is missing")
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.seed = self.root / "seed.h5"; self.seed.write_text("imported")
        self.runner = self.root / "runner.exe"; self.runner.write_text("runner")
        self.identity = {"files": M.hashes([self.runner, self.seed]), "config_sha256": "frozen"}
        self.backend = AnalyticBackend()

    def campaign(self, targets=(1.1,), resume=False, **options):
        return M.Campaign(self.root / "out", self.identity, list(targets), self.seed,
                          1., 1, self.backend, resume=resume, **options)

    def test_analytic_root_and_completed_resume_no_duplicate_solve(self):
        run = self.campaign()
        self.assertTrue(run.run())
        last = run.snapshot["outer_ledger"][-1]
        self.assertLessEqual(abs(last["residual_V"]), 1e-6)
        calls = len(self.backend.calls)
        self.assertTrue(self.campaign(resume=True).run())
        self.assertEqual(len(self.backend.calls), calls)

    def test_resume_rejects_changed_solver_environment(self):
        settings = {key: "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "VELA_LINEAR_THREADS", "VELA_BLAS_THREADS")}
        settings["VELA_LINEAR_SOLVER"] = "umfpack"
        with patch.dict("os.environ", settings):
            self.campaign()
            for key in settings:
                with self.subTest(key=key), patch.dict("os.environ", {key: "changed"}):
                    with self.assertRaisesRegex(M.EvidenceError, "identity changed"):
                        self.campaign(resume=True)
        self.assertEqual(self.backend.calls, [])

    def test_environment_drift_stops_before_device_execution(self):
        with patch.dict("os.environ", {"VELA_LINEAR_SOLVER": "umfpack"}):
            run = self.campaign()
        with patch.dict("os.environ", {"VELA_LINEAR_SOLVER": "sparselu"}):
            with self.assertRaisesRegex(M.EvidenceError, "environment"):
                run.run()
        self.assertEqual(self.backend.calls, [])

    def test_orphan_attempt_directories_cannot_exceed_total_budget(self):
        run = self.campaign(max_attempts=2)
        for index in (1, 2):
            (run.output / "attempts" / f"attempt_{index:06d}").mkdir()
        self.assertFalse(run.run())
        self.assertEqual(self.backend.calls, [])
        self.assertFalse((run.output / "attempts/attempt_000003").exists())

    def test_backward_root_and_both_direction_numerical_recovery(self):
        for target in (1.1, .94):
            with self.subTest(target=target):
                self.root = Path(self.tmp.name) / str(target); self.root.mkdir()
                def limited(bias, seed):
                    if abs(bias - seed["inner_V"]) > .02:
                        raise M.NumericalFailure("device_numerical")
                self.backend = AnalyticBackend(fail=limited)
                run = self.campaign((target,), max_attempts=100, max_shots=20)
                self.assertTrue(run.run())
                calls = self.backend.calls
                self.assertTrue(any(0 < abs(x - seed) <= .02 for x, seed in calls))
                self.assertTrue(any(abs(x - seed) > .02 for x, seed in calls))
                self.assertTrue(all(p["coupled_newton_updates"] > 0 for p in run.snapshot["inner_ledger"][1:]))

    def test_outer_budget_halves_and_retains_original_target(self):
        run = self.campaign(max_shots=1, max_attempts=5)
        self.assertFalse(run.run())
        self.assertIn(1.1, run.snapshot["pending_targets_V"])
        self.assertTrue(any(e.get("event") == "outer_step_halved" for e in run.snapshot["events"]))

    def test_failed_candidate_never_promoted_and_budget_stops(self):
        def reject(bias, seed):
            if bias != seed["inner_V"]: raise M.NumericalFailure("local_audit")
        self.backend.fail = reject
        run = self.campaign(max_attempts=7)
        self.assertFalse(run.run())
        self.assertEqual(len(run.snapshot["inner_ledger"]), 1)
        self.assertEqual(len(run.snapshot["outer_ledger"]), 1)
        self.assertLessEqual(run.snapshot["attempt_count"], 7)

    def test_runtime_failure_stops_without_bisection(self):
        self.backend.fail = lambda bias, seed: (_ for _ in ()).throw(RuntimeError("missing runtime DLL"))
        run = self.campaign()
        with self.assertRaisesRegex(RuntimeError, "DLL"): run.run()
        self.assertEqual(len(self.backend.calls), 1)
        self.assertFalse(run.snapshot["inner_ledger"])

    def test_interrupted_attempt_gets_new_directory_on_resume(self):
        self.backend.interrupt = True
        run = self.campaign()
        with self.assertRaises(KeyboardInterrupt): run.run()
        first = sorted((self.root / "out" / "attempts").iterdir())
        resumed = self.campaign(resume=True)
        self.assertTrue(resumed.run())
        self.assertTrue(first[0].exists())
        self.assertGreater(len(list((self.root / "out" / "attempts").iterdir())), len(first))

    def test_changed_runner_and_modified_point_metadata_rejected(self):
        run = self.campaign(); self.assertTrue(run.run())
        self.runner.write_text("different")
        with self.assertRaises(M.EvidenceError): self.campaign(resume=True)
        self.runner.write_text("runner")
        checkpoint = self.root / "out/checkpoint.json"
        data = json.loads(checkpoint.read_text()); data["inner_ledger"][0]["current_A_per_um"] = 1.
        checkpoint.write_text(json.dumps(data))
        with self.assertRaises(M.EvidenceError): self.campaign(resume=True)

    def test_changed_confirmation_remains_unaccepted(self):
        self.backend.changed_confirmation = True
        run = self.campaign(max_attempts=5)
        self.assertFalse(run.run())
        self.assertFalse(run.snapshot["outer_ledger"])
        self.assertIn("confirmation_state_changed", run.snapshot["failure"])

    def test_zero_time_budget_runs_no_solver(self):
        run = self.campaign(max_total_seconds=.001)
        with patch.object(M.time, "monotonic", return_value=run.started + 1):
            self.assertFalse(run.run())
        self.assertFalse(self.backend.calls)

    def test_acceptance_checkpoint_window_resumes_without_reconfirming(self):
        run = self.campaign()
        original = run.save
        def interrupt_after_acceptance():
            original()
            if len(run.snapshot["outer_ledger"]) == 2:
                raise KeyboardInterrupt()
        run.save = interrupt_after_acceptance
        with self.assertRaises(KeyboardInterrupt): run.run()
        confirmation_count = len(list((self.root / "out/attempts").glob("*/confirmation.json")))
        resumed = self.campaign(resume=True)
        self.assertTrue(resumed.run())
        self.assertEqual(len(list((self.root / "out/attempts").glob("*/confirmation.json"))), confirmation_count)

    def test_inner_recovery_survives_interruption(self):
        calls = []
        def failure(bias, seed):
            calls.append(bias)
            if len(calls) == 2: raise M.NumericalFailure("device_numerical")
            if len(calls) == 3: raise KeyboardInterrupt()
        self.backend.fail = failure
        run = self.campaign()
        with self.assertRaises(KeyboardInterrupt): run.run()
        halfway_bias = calls[-1]
        self.backend.fail = None
        resumed = self.campaign(resume=True)
        resumed.run()
        self.assertEqual(self.backend.calls[3][0], halfway_bias)

    def test_changed_resource_during_campaign_stops_before_solver(self):
        run = self.campaign()
        self.runner.write_text("changed")
        with self.assertRaises(M.EvidenceError): run.run()
        self.assertFalse(self.backend.calls)

    def test_inner_baseline_flag_cannot_be_transferred_to_new_bias(self):
        run = self.campaign(); self.assertTrue(run.run())
        checkpoint = self.root / "out/checkpoint.json"
        data = json.loads(checkpoint.read_text())
        data["inner_ledger"][1]["same_bias_requalification"] = True
        checkpoint.write_text(json.dumps(data))
        with self.assertRaises(M.EvidenceError): self.campaign(resume=True)


class SolverBackendTests(Base):
    def setUp(self):
        super().setUp()
        self.assertTrue(hasattr(M, "SolverBackend"), "real helper adapter is missing")
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.base = deck()
        mesh = self.root / "mesh.json"
        mesh.write_text(json.dumps(h5.mesh()))
        self.base["mesh_file"] = str(mesh)
        self.seed = self.root / "seed.h5"
        h5.write(self.seed, [{"node_id": i, "psi": 1, "phin": 0, "phip": 0,
                              "electrons_m3": 1, "holes_m3": 1} for i in range(3)])
        self.calls = []
        self.zero_updates = False
        self.change_confirmation = False
        self.audit_wrong_state = False
        self.backend = M.SolverBackend(self.base, self.root, self.root / "runner.exe",
                                      expected_node_count=3, expected_active_count=2)
        self.backend.raw_execute = self.execute

    def execute(self, runner, repo, config, run_dir, timeout):
        run_dir.mkdir(parents=True, exist_ok=True)
        document = json.loads(config.read_text())
        self.calls.append(document)
        output = (config.parent / document["output_csv"]).resolve()
        if document.get("simulation_type") == "newton_carrier_term_probe":
            csv_write(output, term_rows())
            if self.audit_wrong_state:
                document["state_file"] = str(self.seed)
                config.write_text(json.dumps(document))
            return 0, False, .001
        bias = document["sweep"]["bias_points"][0]
        circuit = document["sweep"]["external_circuit"]
        state = (config.parent / document["sweep"]["write_state_file"]).resolve()
        state.parent.mkdir(parents=True, exist_ok=True)
        if circuit["enabled"]:
            inner = circuit["initial_inner_voltage_V"]
            row = iv_row(inner, inner * .02 / 1e12)
            row.update(solver_method="boundary_seed_state", newton_convergence_reason="boundary_seed_state_reuse",
                       inner_voltage_V=str(inner), outer_voltage_V=str(bias),
                       series_resistance_ohm_um="1e12", load_line_residual_V=str(inner * 1.02 - bias),
                       boundary_control_mode="external_resistor")
            fields, meta = h5.state_archive.inspect(Path(document["sweep"]["initial_state_file"]))
            meta["confirmation_provenance"] = "different container, same fields"
            h5.state_archive.write(state, fields, meta)
            if self.change_confirmation:
                fields, meta = h5.state_archive.inspect(state)
                fields["psi"][0] += 1e-12
                h5.state_archive.write(state, fields, meta)
        else:
            row = iv_row(bias, bias * .02 / 1e12)
            if self.zero_updates: row["newton_iterations"] = "0"
            h5.write(state, [{"node_id": i, "psi": bias, "phin": 0, "phip": 0,
                              "electrons_m3": 1, "holes_m3": 1} for i in range(3)])
        csv_write(output, [row])
        return 0, False, .001

    def solve(self, bias=1.1, same_bias=False):
        directory = self.root / "attempt"; directory.mkdir()
        return self.backend.solve(bias, {"inner_V": 1., "state_file": str(self.seed)}, directory, same_bias=same_bias)

    def test_full_solver_audit_and_byte_identical_confirmation(self):
        point = self.solve()
        point["id"] = 0
        self.backend.validate(point)
        self.assertEqual(len(self.calls), 2)
        confirm = self.backend.confirm(point, 1.122, self.root / "confirm")
        self.assertEqual(confirm["device_state_sha256"], point["state_sha256"])
        self.assertNotEqual(confirm["state_sha256"], point["state_sha256"])
        self.backend.validate_confirmation(point, 1.122, confirm)
        for call in (self.calls[0], self.calls[-1]):
            solver = M.clone(call["solver"])
            solver["carrier_row_convergence"].pop("diagnostic_csv", None)
            self.assertEqual(solver, self.base["solver"])

    def test_changed_confirmation_state_cannot_inherit_proof(self):
        point = self.solve(); self.change_confirmation = True
        with self.assertRaisesRegex(M.NumericalFailure, "confirmation_state_changed"):
            self.backend.confirm(point, 1.122, self.root / "confirm")

    def test_new_voltage_without_update_is_not_qualified(self):
        self.zero_updates = True
        with self.assertRaisesRegex(M.NumericalFailure, "zero_coupled_update"):
            self.solve()

    def test_same_bias_import_is_actually_solved_and_audited(self):
        self.zero_updates = True
        self.solve(1., True)
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(self.calls[1]["simulation_type"], "newton_carrier_term_probe")

    def test_audit_bound_to_wrong_state_is_corrupt_proof(self):
        self.audit_wrong_state = True
        with self.assertRaisesRegex(M.EvidenceError, "audit.*state"):
            self.solve()

    def test_nonzero_exit_without_outputs_is_not_voltage_failure(self):
        self.backend.raw_execute = lambda *args: (2, False, .001)
        with self.assertRaisesRegex(RuntimeError, "runtime|invalid_input"):
            self.solve()

    def test_timeout_is_bounded_retryable_failure(self):
        self.backend.raw_execute = lambda *args: (-1, True, .001)
        with self.assertRaisesRegex(M.NumericalFailure, "timeout"):
            self.solve()

    def test_confirmation_metadata_current_cannot_be_substituted(self):
        point = self.solve()
        confirmation = self.backend.confirm(point, 1.122, self.root / "confirm")
        path = Path(confirmation["iv_file"])
        row = M.read_rows(path)[0]; row["current_total_A_per_um"] = "1"
        csv_write(path, [row])
        with self.assertRaises(M.EvidenceError):
            self.backend.validate_confirmation(point, 1.122, confirmation)

    def test_attempt_paths_do_not_repeat_long_stage_name(self):
        point = self.solve()
        self.assertLess(len(str(Path(point["state_file"]).relative_to(self.root))), 50)

    def test_corrector_counts_only_main_and_full_validation_updates(self):
        original = self.execute
        def main_fails(*args):
            result = original(*args)
            document = json.loads(args[2].read_text())
            if document.get("simulation_type", "dc_sweep") == "dc_sweep" and len(self.calls) == 1:
                output = (args[2].parent / document["output_csv"]).resolve()
                row = M.read_rows(output)[0]; row["converged"] = "0"
                csv_write(output, [row])
            return result
        self.backend.raw_execute = main_fails
        self.backend.corrector_passes = 1
        def correction(document, state, bias, runner, repo, directory, timeout, config_directory, *rest):
            self.zero_updates = True
            rendered = self.backend.helper.render_attempt(document, directory / "validator", bias, state, 2)
            document, config, run_dir, passed = self.backend._write_run(rendered)
            self.assertTrue(passed)
            final = run_dir / "final_state.h5"
            self.backend.helper.audit_local_candidate(document, final, bias, runner, repo, run_dir, timeout, config.parent)
            (run_dir / "config.json").write_bytes(config.read_bytes())
            return {"passed": True, "state_file": str(final), "coupled_newton_iterations": 999}
        self.backend.helper.correct_local_candidate = correction
        point = self.solve()
        self.assertEqual(point["coupled_newton_updates"], 1)
        self.backend.validate(point)


class CliTests(Base):
    def test_cli_requires_explicit_inner_voltage_and_exposes_budgets(self):
        self.assertTrue(hasattr(M, "parser"), "CLI is missing")
        args = M.parser().parse_args(["--repo-root", ".", "--runner", "runner.exe",
            "--base-config", "base.json", "--initial-state", "state.h5", "--initial-inner-bias", "37.03",
            "--outer-targets", "37.075", "37.1", "--output-dir", "out", "--resume",
            "--max-attempts", "80", "--carrier-corrector-passes", "12", "--restore-poly-equilibrium"])
        self.assertEqual(args.outer_targets, [37.075, 37.1])
        self.assertEqual(args.initial_inner_bias, 37.03)
        self.assertEqual(args.max_attempts, 80)

    def test_checkpoint_atomicity_keeps_last_good_on_replace_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "checkpoint.json"
            M.atomic_json(path, {"accepted": 1})
            with patch.object(M.os, "replace", side_effect=OSError("interrupted")):
                with self.assertRaises(OSError): M.atomic_json(path, {"accepted": 2})
            self.assertEqual(json.loads(path.read_text()), {"accepted": 1})


if __name__ == "__main__":
    unittest.main()
