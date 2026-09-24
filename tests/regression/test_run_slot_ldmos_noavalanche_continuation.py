from __future__ import annotations

import importlib.util
import tempfile
import unittest
import csv
import json
import subprocess
from unittest.mock import patch
from pathlib import Path
from tests.regression import slot_state_fixture as h5


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "run_slot_ldmos_noavalanche_continuation.py"
SPEC = importlib.util.spec_from_file_location("slot_adaptive", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SlotAdaptiveContinuationTest(unittest.TestCase):
    def test_hdf5_restart_rejects_wrong_mesh_origin_mode_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "seed.h5"
            rows = [dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)]
            h5.write(path, rows)
            doc = {"mesh_file": "mesh.json"}
            MODULE.states.read(path, doc, root)
            with self.assertRaisesRegex(ValueError, "origin"):
                MODULE.states.read(path, dict(doc, potential_origin_V=1.), root)
            geometry = h5.mesh(1)
            geometry['nodes'][0]['x'] = 1.
            (root / 'mesh.json').write_text(json.dumps(geometry))
            with self.assertRaisesRegex(ValueError, "Mesh identity"):
                MODULE.states.read(path, doc, root)
            h5.write(path, rows)
            with h5.state_archive.h5py.File(path, 'r+') as handle:
                handle['fields/psi'][0] = float('nan')
            with self.assertRaisesRegex(ValueError, "Invalid field"):
                MODULE.states.read(path, doc, root)
            h5.write(path, rows)
            with h5.state_archive.h5py.File(path, 'r+') as handle:
                handle['fields/psi'].attrs['unit'] = 'm'
            with self.assertRaisesRegex(ValueError, "units mismatch"):
                MODULE.states.read(path, doc, root)
            h5.write(path, rows)
            fields, metadata = h5.state_archive.inspect(path)
            fields['temperature_K'] = [300.]
            metadata['mode'] = 'electrothermal'
            h5.state_archive.write(path, fields, metadata)
            with self.assertRaisesRegex(ValueError, "mode/potential origin"):
                MODULE.states.read(path, doc, root)
            legacy = root / 'seed.csv'
            legacy.write_text('node_id,psi\n0,0\n')
            with self.assertRaises(ValueError):
                MODULE.states.read(legacy, doc, root)

    def test_secant_records_unqualified_parents_and_cannot_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous, current, output = (root / name for name in ('a.h5', 'b.h5', 'out.h5'))
            row = dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)
            for path in (previous, current):
                h5.write(path, [row], state_role='accepted', qualified=True)
            doc = dict(mesh_file='mesh.json', sweep={'contact':'drain'}, contacts=[{'name':'drain','bias':0.}])
            MODULE.write_secant_restart(previous, current, output, 0., 1., 2., document=doc, config_directory=root)
            _, metadata = MODULE.states.read(output, doc, root)
            self.assertFalse(metadata['qualified'])
            self.assertEqual(metadata['state_role'], 'secant_candidate')
            self.assertEqual(metadata['bias_V'], 2.)
            self.assertEqual(len(metadata['parent_state_sha256']), 2)
            before = current.read_bytes()
            with self.assertRaises(FileExistsError):
                MODULE.write_secant_restart(previous, current, current, 0., 1., 2., document=doc, config_directory=root)
            self.assertEqual(current.read_bytes(), before)

    def test_avalanche_continuation_preserves_source_and_separates_evidence(self) -> None:
        impact = {"model": "van_overstraeten", "coupling_mode": "self_consistent",
                  "current_approximation": "element_edge_sg_gss_laux", "source_jacobian": "local_ad"}
        base = {"solver": {"impact_ionization": impact}, "sweep": {}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            document, config, stage = MODULE.render_attempt(base, root, 2.0, root / "seed.csv", 1)
            self.assertEqual(document["solver"]["impact_ionization"], impact)
            self.assertEqual(MODULE.physics_branch(document), "avalanche")
            self.assertIn("02_delaunay_avalanche_adaptive", str(stage))
            self.assertEqual(config.name, "adaptive_avalanche_attempt_001.json")
        base["solver"]["impact_ionization"]["coupling_mode"] = "postprocess_only"
        self.assertEqual(MODULE.physics_branch(base), "noavalanche")
        self.assertEqual(MODULE.physics_branch({}), "noavalanche")

    def test_external_circuit_rejects_outer_voltage_as_local_device_bias(self) -> None:
        # Saved SLOT decks retain all resistor settings with enabled=false.
        MODULE.require_direct_device_bias({"sweep": {"external_circuit": {
            "enabled": False, "mode": "series_resistor", "resistance_ohm_um": 1e12}}})
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); config = root / "base.json"
            config.write_text(json.dumps({"sweep": {"external_circuit": {"mode": "series_resistor"}}}))
            import argparse
            with self.assertRaisesRegex(ValueError, "solved inner voltage"):
                MODULE.run(argparse.Namespace(base_config=config))

    def test_attempt_render_keeps_bundle_relative_outputs_and_absolute_restart(self) -> None:
        base = {
            "output_csv": "old.csv",
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old_state.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_attempt(
                {**base, "solver": {
                    "max_iter": 80,
                    "handoff": {"newton_max_iter": 80},
                    "carrier_row_convergence": {"min_newton_max_iter": 80},
                }}, bundle, 0.2, initial, 3, 8, 1.0e-6,
                "uniform_trust_region", "block_filter", "gummel_newton", 11)
            self.assertEqual(document["sweep"]["bias_points"], [0.2])
            self.assertEqual(document["sweep"]["initial_state_file"], str(initial.resolve()))
            self.assertTrue(document["output_csv"].endswith("iv.csv"))
            self.assertEqual(config.parent, bundle)
            self.assertIn("attempt_003_v_0p2", str(run_dir))
            self.assertEqual(document["solver"]["max_iter"], 8)
            self.assertEqual(document["solver"]["stall_residual_floor"], 1.0e-6)
            self.assertEqual(
                document["solver"]["quasi_fermi_update_limit_mode"],
                "uniform_trust_region")
            self.assertEqual(document["solver"]["line_search_mode"], "block_filter")
            self.assertEqual(document["solver"]["method"], "gummel_newton")
            self.assertEqual(document["solver"]["handoff"]["gummel_max_iter"], 11)
            self.assertEqual(document["solver"]["handoff"]["newton_max_iter"], 8)
            self.assertEqual(
                document["solver"]["carrier_row_convergence"]["min_newton_max_iter"], 8)
            self.assertTrue(
                document["solver"]["carrier_row_convergence"]["diagnostic_csv"]
                .endswith("carrier_row_convergence.csv"))

    def test_poisson_prewarm_preserves_qf_and_writes_separate_state(self) -> None:
        base = {
            "output_csv": "old.csv",
            "solver": {
                "method": "newton",
                "max_iter": 80,
                "handoff": {"newton_max_iter": 80},
            },
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old_state.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_poisson_prewarm(
                base, bundle, 0.246875, initial, 8, 17)
            self.assertEqual(document["solver"]["method"], "poisson_only")
            self.assertEqual(
                document["solver"]["poisson_only_quasi_fermi_initialization"],
                "preserve_initial")
            self.assertEqual(document["solver"]["max_iter"], 17)
            self.assertEqual(document["sweep"]["bias_points"], [0.246875])
            self.assertEqual(
                document["sweep"]["initial_state_file"], str(initial.resolve()))
            self.assertIn("prewarm_008_v_0p246875", str(run_dir))
            self.assertIn("prewarm_008", config.name)
            basin, _, _ = MODULE.render_poisson_prewarm(
                base, bundle, 0.246875, initial, 8, 17, "contact_basin")
            self.assertEqual(basin["solver"]["poisson_only_quasi_fermi_initialization"], "contact_basin")
            with self.assertRaises(ValueError):
                MODULE.render_poisson_prewarm(base, bundle, 0.246875, initial, 8, 17, "typo")

    def test_local_update_diagnostics_are_opt_in_and_run_relative(self) -> None:
        document = {"solver": {"max_iter": 7}}
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp) / "attempt"
            MODULE.enable_local_update_diagnostics(
                document, run_dir, [6097, 10163])
            diagnostics = document["solver"]["local_update_diagnostics"]
            self.assertTrue(document["solver"]["diagnostic_history"])
            self.assertTrue(diagnostics["enabled"])
            self.assertEqual(diagnostics["nodes"], [6097, 10163])
            self.assertEqual(diagnostics["first_iterations"], 7)
            self.assertEqual(diagnostics["every_iterations"], 1)
            self.assertEqual(
                Path(diagnostics["csv_file"]),
                (run_dir / "local_updates.csv").resolve())

    def test_qf_limit_overrides_support_explicit_disable(self) -> None:
        document = {
            "solver": {
                "quasi_fermi_update_limit_V": 0.1,
                "quasi_fermi_update_limit_minority_V": 0.05,
            }
        }
        MODULE.apply_qf_limit_overrides(document, 0.0, 0.0)
        self.assertEqual(document["solver"]["quasi_fermi_update_limit_V"], 0.0)
        self.assertEqual(
            document["solver"]["quasi_fermi_update_limit_minority_V"], 0.0)

    def test_carrier_diagonal_floor_override_is_explicit(self) -> None:
        document = {"solver": {}}
        MODULE.apply_carrier_diagonal_floor_override(
            document, 2.0, 0.25)
        floor = document["solver"]["carrier_diagonal_floor"]
        self.assertTrue(floor["enabled"])
        self.assertEqual(floor["scale"], 2.0)
        self.assertEqual(floor["minority_density_ratio"], 0.25)

    def test_rejected_final_state_is_exposed_for_chunked_restart(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            state = root / "rejected_states" / "attempt_final.csv"
            state.parent.mkdir()
            state.write_text("node_id,psi\n0,0\n", encoding="utf-8")
            attempts = root / "attempts.csv"
            with attempts.open("w", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(
                    stream, fieldnames=["status", "rejected_final_state_file"])
                writer.writeheader()
                writer.writerow({
                    "status": "rejected",
                    "rejected_final_state_file": str(state),
                })
            self.assertEqual(
                MODULE.load_rejected_final_state(attempts), state.resolve())

            state.unlink()
            self.assertIsNone(MODULE.load_rejected_final_state(attempts))

    def test_failure_state_policy_overrides_base_deck(self) -> None:
        document = {
            "solver": {"return_best_feasible_on_failure": True},
        }
        MODULE.apply_failure_state_policy(document, False)
        self.assertFalse(
            document["solver"]["return_best_feasible_on_failure"])
        MODULE.apply_failure_state_policy(document, True)
        self.assertTrue(
            document["solver"]["return_best_feasible_on_failure"])

    def test_batch_predictor_keeps_secant_history_in_one_sweep(self) -> None:
        base = {
            "output_csv": "old.csv",
            "solver": {"max_iter": 80},
            "sweep": {
                "bias_points": [60.0],
                "initial_state_file": "old.csv",
                "write_state_every_point_prefix": "old/state",
                "write_state_file": "old/final.csv",
            },
        }
        with tempfile.TemporaryDirectory() as tmp:
            bundle = Path(tmp)
            initial = bundle / "seed.csv"
            document, config, run_dir = MODULE.render_batch_predictor(
                base, bundle, [0.2, 0.35, 0.5], initial, 12, 1.0e-6,
                "componentwise_poisson_recorrect", "merit")
            self.assertEqual(document["sweep"]["bias_points"], [0.2, 0.35, 0.5])
            self.assertEqual(
                document["sweep"]["continuation"]["predictor"]["mode"], "secant")
            self.assertEqual(
                document["sweep"]["continuation"]["predictor"]["fields"],
                ["psi", "phin", "phip"])
            self.assertIn("batch_predictor_v_0p5", str(run_dir))
            self.assertEqual(config.name, "adaptive_noavalanche_batch_predictor.json")

    def test_secant_restart_extrapolates_only_coupled_unknowns(self) -> None:
        fields = [
            "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3",
            "electron_qf_increment_V", "hole_qf_increment_V",
            "electron_qf_reference_V", "hole_qf_reference_V"]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous = root / "previous.h5"
            current = root / "current.h5"
            output = root / "predicted.h5"
            rows = [
                ["0", "1", "2", "3", "10", "20", "2", "3", "0", "0"],
                ["0", "2", "4", "6", "11", "21", "4", "6", "0", "0"],
            ]
            for path, row in ((previous, rows[0]), (current, rows[1])):
                h5.write(path, [dict(zip(fields, row))])
            ratio = MODULE.write_secant_restart(
                previous, current, output, 0.2, 0.35, 0.5,
                document={"mesh_file": "mesh.json"}, config_directory=root)
            self.assertAlmostEqual(ratio, 1.0)
            predicted = h5.read(output)[0]
            self.assertAlmostEqual(float(predicted["psi"]), 3.0)
            self.assertAlmostEqual(float(predicted["phin"]), 6.0)
            self.assertAlmostEqual(float(predicted["phip"]), 9.0)
            self.assertEqual(float(predicted["electrons_m3"]), 11.0)
            self.assertEqual(float(predicted["holes_m3"]), 21.0)
            self.assertAlmostEqual(float(predicted["electron_qf_increment_V"]), 6.0)
            self.assertAlmostEqual(float(predicted["hole_qf_increment_V"]), 9.0)


    def test_secant_preserves_sub_ulp_current_signal_on_a_moving_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            previous, current, output = (root / name for name in ("previous.h5", "current.h5", "out.h5"))
            for path, bias, increment in ((previous, 1.0, 1e-20), (current, 2.0, 2e-20)):
                row = {"node_id": 0, "psi": bias, "phin": bias,
                       "phip": -increment, "electrons_m3": 10, "holes_m3": 20,
                       "electron_qf_reference_V": bias, "electron_qf_increment_V": increment,
                       "hole_qf_reference_V": 0, "hole_qf_increment_V": -increment}
                h5.write(path, [row])
            MODULE.write_secant_restart(previous, current, output, 1.0, 2.0, 3.0,
                document={"mesh_file": "mesh.json"}, config_directory=root)
            row = h5.read(output)[0]
            self.assertEqual(float(row["phin"]), 3.0)
            self.assertEqual(float(row["electron_qf_reference_V"]), 3.0)
            self.assertAlmostEqual(float(row["electron_qf_increment_V"]) / 1e-20, 3.0)
            self.assertAlmostEqual(float(row["hole_qf_increment_V"]) / 1e-20, -3.0)

    def test_loader_failure_stops_bisection_and_records_environment_evidence(self) -> None:
        for code in (0xC0000135, -1073741515, 0xC0000139, 0xC000007B, 0xC0000142):
            with self.subTest(code=code), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                result = subprocess.CompletedProcess([], code, "", "")
                with patch.object(MODULE.subprocess, "run", return_value=result):
                    with self.assertRaisesRegex(RuntimeError, "startup failed"):
                        MODULE.execute_config(root / "runner", root, root / "config.json", root / "run", 1)
                failure = json.loads((root / "run/startup_failure.json").read_text())
                self.assertEqual(failure["failure_class"], "runtime_environment")
                self.assertEqual(failure["exit_code"], code)

    def test_solver_failure_is_still_a_retryable_numerical_attempt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            result = subprocess.CompletedProcess([], 1, "solver failed", "")
            with patch.object(MODULE.subprocess, "run", return_value=result):
                code, timed_out, _ = MODULE.execute_config(
                    root / "runner", root, root / "config.json", root / "run", 1)
            self.assertEqual(code, 1)
            self.assertFalse(timed_out)
            self.assertFalse((root / "run/startup_failure.json").exists())

    def test_timeout_preserves_byte_output_and_returns_a_rejected_attempt(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            error = subprocess.TimeoutExpired("runner", 0.01, output=b"partial output", stderr=b"partial error")
            with patch.object(MODULE.subprocess, "run", side_effect=error):
                code, timed_out, elapsed = MODULE.execute_config(
                    root / "runner", root, root / "config.json", root / "run", 0.01)
            self.assertEqual(code, -1)
            self.assertTrue(timed_out)
            self.assertGreaterEqual(elapsed, 0)
            self.assertEqual((root / "run/stdout.log").read_text(), "partial output")
            self.assertEqual((root / "run/stderr.log").read_text(), "partial error")


    def test_audit_rejects_cancelling_rows_and_preserves_input_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            document = {"mesh_file": "mesh.json", "node_doping_file": "doping.csv",
                        "materials_file": "materials.json", "sweep": {"contact": "drain"},
                        "contacts": [{"name": "drain", "bias": 0.0}],
                        "solver": {"carrier_row_convergence": {
                            "eps_row": 0.001, "min_source_scale": 1e-13, "scale_floor": 1e-30}}}
            def probe(runner, repo, config, run_dir, timeout):
                deck = json.loads(config.read_text())
                self.assertEqual(deck["mesh_file"], str((root / "base/mesh.json").resolve()))
                self.assertEqual(deck["contacts"][0]["bias"], 2.0)
                self.assertNotIn("sweep", deck)
                rows = []
                for node, residual in enumerate((1e-15, -1e-15)):
                    row = {"node_id": node, "x": node, "y": 0}
                    for carrier in ("electron", "hole"):
                        row.update({f"{carrier}_{key}": value for key, value in {
                            "continuity_active": 1, "residual": residual,
                            "flux": residual - 1e-15, "flux_abs_sum": 2e-15,
                            "recombination": 1e-15, "impact": 0, "gauge": 0,
                            "boundary": 0, "term_sum": residual}.items()})
                    rows.append(row)
                with Path(deck["output_csv"]).open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
                    writer.writeheader()
                    writer.writerows(rows)
                return 0, False, 0.1
            with patch.object(MODULE, "execute_config", side_effect=probe):
                result = MODULE.audit_local_candidate(document, root / "state.csv", 2.0,
                    root / "runner", root, root / "attempt", 10, root / "base")
            self.assertFalse(result["passed"])
            self.assertEqual(result["metrics"]["carriers"]["electron"]["residual_signed_sum"], 0)
            self.assertEqual(document["contacts"][0]["bias"], 0)

    def test_cold_avalanche_cannot_advance_without_coupled_and_all_row_qualification(self):
        for coupled_ok in (False, True):
            with self.subTest(coupled_ok=coupled_ok), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                row = dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)
                h5.write(root/'unused_seed.h5', [row])
                impact = dict(model='van_overstraeten', coupling_mode='self_consistent',
                              current_approximation='ialmob_element_density_grad_qf', source_jacobian='local_ad')
                base = root/'base.json'
                base.write_text(json.dumps(dict(mesh_file='mesh.json', state_file='must_not_load.h5',
                    state_fields_dir='must_not_load', solver=dict(method='newton', impact_ionization=impact),
                    contacts=[dict(name='drain', bias=0.)], sweep=dict(contact='drain', initial_state_file='must_not_load.h5'))))
                methods = []
                def execute(runner, repo, config, run_dir, timeout):
                    deck = json.loads(config.read_text())
                    methods.append(deck['solver']['method'])
                    self.assertEqual(deck['solver']['impact_ionization'], impact)
                    self.assertEqual(deck['sweep']['bias_points'], [0.])
                    if methods[-1] == 'poisson_only':
                        self.assertNotIn('state_file', deck)
                        self.assertNotIn('state_fields_dir', deck)
                        self.assertNotIn('initial_state_file', deck['sweep'])
                    run_dir.mkdir(parents=True, exist_ok=True)
                    h5.write(run_dir/'final_state.h5', [row])
                    passed = methods[-1] == 'poisson_only' or coupled_ok
                    (run_dir/'iv.csv').write_text(f'converged,newton_iterations\n{int(passed)},0\n')
                    return 0, False, .1
                argv = ['continuation', '--repo-root', str(root), '--runner', str(root/'runner'),
                        '--base-config', str(base), '--cold-start', '--targets', '.025']
                with patch('sys.argv', argv), patch.object(MODULE, 'execute_config', side_effect=execute), \
                     patch.object(MODULE, 'audit_local_candidate', return_value={'passed':False}) as audit:
                    self.assertEqual(MODULE.main(), 1)
                self.assertEqual(methods, ['poisson_only', 'newton'])
                self.assertEqual(audit.call_count, int(coupled_ok))
                result = json.loads((root/'reports/adaptive_avalanche_summary.json').read_text())
                self.assertIsNone(result['last_successful_bias_V'])
                self.assertIsNone(result['last_successful_state'])
                self.assertEqual(result['remaining_targets_V'], [0., .025])

    def test_timed_out_solve_keeps_target_pending_without_voltage_bisection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root / 'seed.h5'
            h5.write(seed, [dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)])
            base = root / 'base.json'
            base.write_text(json.dumps(dict(mesh_file='mesh.json', solver={},
                contacts=[dict(name='drain', bias=0.)], sweep=dict(contact='drain'))))
            argv = ['continuation', '--repo-root', str(root), '--runner', str(root/'runner'),
                    '--base-config', str(base), '--initial-state', str(seed),
                    '--initial-bias', '0', '--targets', '1', '--timeout-seconds', '1']
            with patch('sys.argv', argv), patch.object(MODULE, 'execute_config',
                    return_value=(-1, True, 1.)) as execute:
                self.assertEqual(MODULE.main(), 1)
            self.assertEqual(execute.call_count, 1)
            result = json.loads((root/'reports/adaptive_noavalanche_summary.json').read_text())
            self.assertEqual(result['termination_reason'], 'timeout')
            self.assertEqual(result['remaining_targets_V'], [1.])
            self.assertIsNone(result['last_successful_bias_V'])
            self.assertIsNone(result['last_successful_state'])
            self.assertEqual(result['initial_bias_V'], 0.)
            self.assertEqual(Path(result['initial_state']), seed)

    def test_corrector_timeout_stops_at_each_stage(self):
        for stage in ('block_probe', 'candidate_audit', 'coupled_validation', 'final_audit'):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                document = dict(simulation_type='dc_sweep', solver={},
                    sweep=dict(contact='drain'), contacts=[dict(name='drain', bias=0.)])
                def execute(runner, repo, config, run_dir, timeout):
                    deck = json.loads(config.read_text())
                    probe = deck['simulation_type'] == 'newton_block_step_probe'
                    if stage == ('block_probe' if probe else 'coupled_validation'):
                        return -1, True, 1.
                    if probe:
                        (run_dir/'trial_carrier_only.h5').write_text('candidate')
                    else:
                        (run_dir/'final_state.h5').write_text('coupled candidate')
                        (run_dir/'iv.csv').write_text('converged,newton_iterations\n1,1\n')
                    return 0, False, .1
                audits = [dict(passed=stage != 'candidate_audit', timed_out=stage == 'candidate_audit'),
                          dict(passed=stage != 'final_audit', timed_out=stage == 'final_audit')]
                with patch.object(MODULE, 'execute_config', side_effect=execute), \
                     patch.object(MODULE, 'audit_local_candidate', side_effect=audits):
                    result = MODULE.correct_local_candidate(document, root/'seed.h5', 1.,
                        root/'runner', root, root/'attempt', 1., root, 4)
                self.assertTrue(result['timed_out'])
                self.assertEqual(result['timeout_stage'], stage)
                self.assertFalse(result['passed'])
                self.assertEqual(result['state_file'], '')
                self.assertEqual(len(result['steps']), 1)

    def test_audit_timeout_keeps_target_without_correction_or_bisection(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root/'seed.h5'
            h5.write(seed, [dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)])
            base = root/'base.json'
            base.write_text(json.dumps(dict(mesh_file='mesh.json', solver={},
                contacts=[dict(name='drain', bias=0.)], sweep=dict(contact='drain'))))
            argv = ['continuation', '--repo-root', str(root), '--runner', str(root/'runner'),
                    '--base-config', str(base), '--initial-state', str(seed), '--initial-bias', '0',
                    '--targets', '1', '--audit-local-continuity', '--carrier-corrector-passes', '4']
            def execute(runner, repo, config, run_dir, timeout):
                run_dir.mkdir(parents=True, exist_ok=True)
                (run_dir/'final_state.h5').write_text('candidate')
                (run_dir/'iv.csv').write_text('converged,newton_iterations\n1,1\n')
                return 0, False, .1
            with patch('sys.argv', argv), patch.object(MODULE, 'execute_config', side_effect=execute) as run, \
                 patch.object(MODULE, 'audit_local_candidate', return_value=dict(passed=False, timed_out=True)), \
                 patch.object(MODULE, 'correct_local_candidate') as correct:
                self.assertEqual(MODULE.main(), 1)
            self.assertEqual(run.call_count, 1)
            correct.assert_not_called()
            result = json.loads((root/'reports/adaptive_noavalanche_summary.json').read_text())
            self.assertEqual(result['termination_reason'], 'timeout')
            self.assertEqual(result['remaining_targets_V'], [1.])
            self.assertIsNone(result['last_successful_state'])
            self.assertEqual(result['attempts'][0]['timeout_stage'], 'local_audit')

    def test_interrupted_next_attempt_preserves_qualified_recovery_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seed = root/'seed.h5'
            h5.write(seed, [dict(node_id=0, psi=0., phin=0., phip=0., electrons_m3=1., holes_m3=1.)])
            base = root/'base.json'
            base.write_text(json.dumps(dict(mesh_file='mesh.json', solver={},
                contacts=[dict(name='drain', bias=0.)], sweep=dict(contact='drain'))))
            argv = ['continuation', '--repo-root', str(root), '--runner', str(root/'runner'),
                    '--base-config', str(base), '--initial-state', str(seed), '--initial-bias', '0',
                    '--targets', '1', '2', '--audit-local-continuity']
            def execute(runner, repo, config, run_dir, timeout):
                target = json.loads(config.read_text())['sweep']['bias_points'][0]
                if target == 2:
                    raise KeyboardInterrupt('simulated interruption of next solve')
                run_dir.mkdir(parents=True, exist_ok=True)
                (run_dir/'final_state.h5').write_bytes(seed.read_bytes())
                (run_dir/'iv.csv').write_text('converged,newton_iterations\n1,1\n')
                return 0, False, .1
            with patch('sys.argv', argv), patch.object(MODULE, 'execute_config', side_effect=execute), \
                 patch.object(MODULE, 'audit_local_candidate', return_value=dict(passed=True, timed_out=False)):
                with self.assertRaises(KeyboardInterrupt):
                    MODULE.main()
            checkpoint = json.loads((root/'reports/adaptive_checkpoint.json').read_text())
            self.assertEqual(checkpoint['last_successful_bias_V'], 1.)
            self.assertEqual(checkpoint['remaining_targets_V'], [2.])
            self.assertTrue(checkpoint['last_finished_attempt']['local_audit_passed'])
            self.assertFalse(checkpoint['completed'])
            state = Path(checkpoint['last_successful_state'])
            self.assertEqual(checkpoint['state_sha256'], MODULE.hashlib.sha256(state.read_bytes()).hexdigest())
            self.assertEqual(Path(checkpoint['previous_state']), seed)

    def test_missing_hdf5_candidate_stops_instead_of_bisecting(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            document = {"simulation_type": "dc_sweep", "sweep": {"contact": "drain"},
                        "contacts": [{"name": "drain", "bias": 0}], "solver": {}}
            with patch.object(MODULE, "execute_config", return_value=(0, False, 0.1)):
                with self.assertRaisesRegex(RuntimeError, "HDF5 candidate"):
                    MODULE.correct_local_candidate(document, root / "seed.h5", 1.,
                        root / "runner", root, root / "attempt", 10, root, 1)
            summary = json.loads((root / "attempt/local_corrector/summary.json").read_text())
            self.assertEqual(summary["failure_reason"], "missing_hdf5_candidate")

    def test_carrier_corrector_requires_full_coupled_and_post_validation_audit(self) -> None:
        for coupled_pass, final_audit, expected in ((False, True, False),
                                                   (True, False, False), (True, True, True)):
            with self.subTest(coupled=coupled_pass, audit=final_audit), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                document = {"simulation_type": "dc_sweep", "mesh_file": "mesh.json",
                            "output_csv": "old.csv", "sweep": {"contact": "drain"},
                            "contacts": [{"name": "drain", "bias": 0}],
                            "solver": {"abstol": 1e-9, "carrier_diagonal_floor": {"enabled": True}}}
                def execute(runner, repo, config, run_dir, timeout):
                    deck = json.loads(config.read_text())
                    self.assertEqual(deck["solver"]["abstol"], 1e-9)
                    if deck["simulation_type"] == "newton_block_step_probe":
                        self.assertFalse(deck["solver"]["carrier_diagonal_floor"]["enabled"])
                        (run_dir / "trial_carrier_only.h5").write_text("trial")
                    else:
                        self.assertTrue(deck["solver"]["carrier_diagonal_floor"]["enabled"])
                        self.assertEqual(deck["mesh_file"], str((root / "mesh.json").resolve()))
                        (run_dir / "final_state.h5").write_text("candidate")
                        (run_dir / "iv.csv").write_text("converged,newton_iterations\n" +
                                                        f"{int(coupled_pass)},2\n")
                    return 0, False, 0.1
                audits = [{"passed": True}, {"passed": final_audit}]
                with patch.object(MODULE, "execute_config", side_effect=execute), \
                     patch.object(MODULE, "audit_local_candidate", side_effect=audits):
                    result = MODULE.correct_local_candidate(document, root / "seed.csv", 1.0,
                        root / "runner", root, root / "attempt", 10, root, 1)
                self.assertEqual(result["passed"], expected)
                self.assertEqual(bool(result["state_file"]), expected)
                self.assertTrue(document["solver"]["carrier_diagonal_floor"]["enabled"])
                self.assertEqual(document["output_csv"], "old.csv")

    def test_poly_restart_is_only_a_candidate_and_retains_both_validation_gates(self) -> None:
        for coupled_pass, final_audit, expected in ((False, True, False),
                                                   (True, False, False), (True, True, True)):
            with self.subTest(coupled=coupled_pass, audit=final_audit), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                document = {"simulation_type": "dc_sweep", "sweep": {"contact": "drain"},
                            "contacts": [{"name": "drain", "bias": 0.0}],
                            "solver": {"abstol": 1e-9}}
                original = json.dumps(document, sort_keys=True)
                def execute(runner, repo, config, run_dir, timeout):
                    deck = json.loads(config.read_text())
                    self.assertEqual(deck["solver"]["abstol"], 1e-9)
                    if deck["simulation_type"] == "newton_block_step_probe":
                        (run_dir / "trial_carrier_only.h5").write_text("original trial")
                    else:
                        self.assertEqual(Path(deck["sweep"]["initial_state_file"]).name,
                                         "poly_equilibrium_candidate.h5")
                        (run_dir / "final_state.h5").write_text("validated state")
                        (run_dir / "iv.csv").write_text(
                            f"converged,newton_iterations\n{int(coupled_pass)},2\n")
                    return 0, False, 0.1
                def restore(doc, source, output, config_directory, bias):
                    self.assertEqual(bias, 15.0)
                    self.assertEqual(source.read_text(), "original trial")
                    self.assertEqual(json.dumps(document, sort_keys=True), original)
                    output.write_text("equilibrium candidate")
                    return {"max_qf_change_V": 1e-20, "components": [{"nodes": 3, "qf_V": 0.0}]}
                with patch.object(MODULE, "execute_config", side_effect=execute), \
                     patch.object(MODULE, "audit_local_candidate", side_effect=[
                         {"passed": True}, {"passed": final_audit}]), \
                     patch.object(MODULE, "restore_poly_candidate", side_effect=restore, create=True) as projection:
                    result = MODULE.correct_local_candidate(document, root / "seed.csv", 15.0,
                        root / "runner", root, root / "attempt", 10, root, 1,
                        restore_poly_equilibrium=True)
                self.assertEqual(result["passed"], expected)
                self.assertEqual(bool(result["state_file"]), expected)
                self.assertEqual(projection.call_count, 1)
                self.assertEqual(result["steps"][0]["poly_equilibrium"]["max_qf_change_V"], 1e-20)
                self.assertEqual(json.dumps(document, sort_keys=True), original)

    def test_corrector_retries_coupled_state_without_accepting_intermediate_candidate(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            document = {"simulation_type": "dc_sweep", "sweep": {"contact": "drain"},
                        "contacts": [{"name": "drain", "bias": 0}], "solver": {}}
            validation_count = 0
            def execute(runner, repo, config, run_dir, timeout):
                nonlocal validation_count
                deck = json.loads(config.read_text())
                if deck["simulation_type"] == "newton_block_step_probe":
                    if validation_count:
                        self.assertEqual(Path(deck["state_file"]).name, "saved_rejected.h5")
                    (run_dir / "trial_carrier_only.h5").write_text("trial")
                else:
                    validation_count += 1
                    passed = validation_count == 2
                    (run_dir / "iv.csv").write_text(
                        f"converged,newton_iterations\n{int(passed)},3\n")
                    if passed:
                        (run_dir / "final_state.h5").write_text("accepted")
                    else:
                        rejected = run_dir / "saved_rejected.h5"
                        rejected.write_text("rejected")
                        with (run_dir / "attempts.csv").open("w", newline="") as stream:
                            writer = csv.DictWriter(stream, fieldnames=["rejected_final_state_file"])
                            writer.writeheader()
                            writer.writerow({"rejected_final_state_file": str(rejected)})
                return 0, False, 0.1
            with patch.object(MODULE, "execute_config", side_effect=execute), \
                 patch.object(MODULE, "audit_local_candidate", return_value={"passed": True}):
                result = MODULE.correct_local_candidate(document, root / "seed.csv", 1.0,
                    root / "runner", root, root / "attempt", 10, root, 2)
            self.assertTrue(result["passed"])
            self.assertEqual(result["coupled_newton_iterations"], 6)
            self.assertEqual(len(result["steps"]), 2)
            self.assertNotEqual(result["steps"][0]["validation_directory"],
                                result["steps"][1]["validation_directory"])

    def test_attempt_diagnostics_do_not_overwrite_base_or_other_attempts(self) -> None:
        base = {"solver": {
                    "global_continuity_closure": {"audit_jsonl_file": "old.jsonl"},
                    "carrier_row_convergence": {"trace_csv": "old_trace.csv"}},
                "sweep": {"diagnostics": {"newton_history": {
                    "attempts_csv_file": "old_attempts.csv", "iterations_csv_file": "old_iterations.csv",
                    "rejected_state_directory": "old_states"}}}}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            one, _, _ = MODULE.render_attempt(base, root, 1, root / "state.csv", 1)
            two, _, _ = MODULE.render_attempt(base, root, 2, root / "state.csv", 2)
            first = one["solver"]["global_continuity_closure"]["audit_jsonl_file"]
            second = two["solver"]["global_continuity_closure"]["audit_jsonl_file"]
            self.assertNotEqual(first, second)
            self.assertEqual(base["solver"]["global_continuity_closure"]["audit_jsonl_file"], "old.jsonl")
            self.assertNotEqual(one["sweep"]["diagnostics"]["newton_history"]["iterations_csv_file"],
                                two["sweep"]["diagnostics"]["newton_history"]["iterations_csv_file"])


if __name__ == "__main__":
    unittest.main()
