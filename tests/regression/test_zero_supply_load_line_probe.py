"""Bounded zero-supply diagnostics must never turn into a voltage sweep."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

RUNNER = None


class ZeroSupplyProbe(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        mesh = {
            "nodes": [{"id": i, "x": x, "y": y} for i, (x, y) in enumerate(
                [(0, 0), (1, 0), (0, 1), (1, 1), (.5, .5)])],
            "triangles": [{"id": i, "node_ids": ids, "region_id": 0}
                          for i, ids in enumerate([[0, 1, 4], [1, 3, 4], [3, 2, 4], [2, 0, 4]])],
            "regions": [{"id": 0, "name": "Si", "material": "Si", "cell_ids": [0, 1, 2, 3]}],
            "contacts": [{"id": 0, "name": "source", "region_id": 0, "node_ids": [0, 2]},
                         {"id": 1, "name": "drain", "region_id": 0, "node_ids": [1, 3]}]}
        (self.root / "mesh.json").write_text(json.dumps(mesh))
        f = self.root / "fields"
        f.mkdir()
        for key in ["ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential"]:
            (f / (key + "_region0.csv")).write_text(
                "node_id,component0\n" + "".join(f"{i},0\n" for i in range(5)))
        self.cfg = {"simulation_type": "newton_zero_supply_load_line_iteration_probe",
                    "mesh_file": "mesh.json", "state_fields_dir": "fields",
                    "output_csv": "rows.csv", "linear_system_output_dir": "data",
                    "contact": "drain", "bordered_linear_backend": "sparselu",
                    "supply_voltage_V": 0, "resistance_ohm_um": 1e12,
                    "scaling": {"mode": "unit_scaling"},
                    "contacts": [{"name": "source", "bias": 0}, {"name": "drain", "bias": 0}],
                    "doping": [{"region": "Si", "donors": 0, "acceptors": 0}],
                    "solver": {"method": "newton", "warm_start": True}}

    def tearDown(self):
        self.tmp.cleanup()

    def run_probe(self):
        p = self.root / "config.json"
        p.write_text(json.dumps(self.cfg))
        return subprocess.run([str(RUNNER), "--config", str(p)], capture_output=True,
                              text=True, timeout=60)

    def test_continuation_retains_internal_voltage_and_is_bounded(self):
        self.cfg["contacts"][1]["bias"] = 1e-6
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        result = json.loads((self.root / "data/summary.json").read_text())
        self.assertEqual(result["initial_voltage_V"], 1e-6)
        self.assertLessEqual(result["accepted_updates"], 1)
        self.assertFalse(result["qualified"])
        self.assertEqual(result["new_sweep_points"], 0)
        self.assertLessEqual(abs(result["final_voltage_V"]), .01)

    def test_split_voltage_column_matches_normalized_boundary(self):
        self.cfg["solver"]["split_dd_state"] = True
        self.cfg["solver"]["contact_boundary_minority_electron_relaxation"] = False
        self.cfg["solver"]["mobility"] = {
            "model": "phumob", "edge_averaging": "element_box_phumob",
            "doping_concentration_basis": "total_impurity"}
        self.cfg["mesh_geometry"] = {"poisson_permittivity_policy": "cell_material"}
        self.cfg["solver"]["region_resolved_interface_assembly"] = {
            "transport_edge_geometry": "element_box"}
        self.cfg["solver"]["recombination"] = ["srh"]
        self.cfg["contacts"][1]["bias"] = 1e-24
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        report = json.loads((self.root / "data/summary.json").read_text())
        path = self.root / "data/physical_bordered.mtx"
        rows = [l for l in path.read_text().splitlines() if l and not l.startswith("%")]
        size = int(rows[0].split()[0])
        column = {int(i)-1: float(v) for i, j, v in (l.split() for l in rows[1:])
                  if int(j) == size and int(i) < size}
        self.assertEqual(set(column), {1, 3, 6, 8, 11, 13})
        for value in column.values():
            self.assertAlmostEqual(value * report["potential_scale_V"], -1., places=13)
        self.assertFalse(report["qualified"])

        # Changing only the finite boundary update must leave the frozen
        # physical/stabilized matrices, RHS and returned raw direction intact.
        self.cfg.update(simulation_type="newton_bounded_supply_load_line_iteration_probe",
                        trial_boundary_policy="rounded_voltage_projection", trace_weak_rows=True,
                        linear_system_output_dir="projected")
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        projected = self.root / "projected"
        result = json.loads((projected / "summary.json").read_text())
        self.assertEqual(result["trial_boundary_policy"], "rounded_voltage_projection")
        self.assertFalse(result["qualified"])
        for name in ["physical_bordered.mtx", "bordered.mtx", "vectors.csv"]:
            self.assertEqual((self.root / "data" / name).read_bytes(), (projected / name).read_bytes())
        self.assertTrue(result["trials"])
        self.assertTrue(all("boundary_rhs_inf" in trial for trial in result["trials"]))

        # Both policies consume one actual solve, not two independently
        # regenerated directions. All alpha candidates survive rejection.
        self.cfg.update(paired_trial_audit=True, linear_system_output_dir="paired")
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        paired = json.loads((self.root / "paired/summary.json").read_text())
        self.assertEqual(paired["linear_solve_count"], 1)
        self.assertEqual(paired["accepted_updates"], 0)
        self.assertEqual(len(paired["audit_candidates"]), 26)
        for candidate in paired["audit_candidates"]:
            self.assertTrue((self.root / "paired" / candidate["state_file"]).is_file())
            self.assertFalse(candidate["qualified"])
            self.assertIn("all_rows_satisfied", candidate)
            self.assertEqual(set(candidate["raw_reaction_ports_A_per_um"]), {"source", "drain"})
        for a, b in zip(paired["audit_candidates"][::2], paired["audit_candidates"][1::2]):
            self.assertEqual(a["alpha"], b["alpha"])
            self.assertEqual(a["requested_voltage_delta_V"], b["requested_voltage_delta_V"])
            self.assertNotEqual(a["policy"], b["policy"])

        self.cfg.update(trial_boundary_policy="split_contact_voltage", linear_system_output_dir="split_voltage")
        self.cfg["solver"]["split_contact_bias_low_V"] = {"drain": 1e-42}
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        split = json.loads((self.root / "split_voltage/summary.json").read_text())
        self.assertEqual(split["initial_voltage_low_V"], 1e-42)
        self.assertEqual(split["final_voltage_low_V"], 1e-42)
        self.assertEqual(len(split["audit_candidates"]), 13)
        self.assertEqual(split["linear_solve_count"], 1)
        self.assertFalse(split["qualified"])
        self.assertTrue(all("voltage_low_V" in t for t in split["audit_candidates"]))

        self.cfg.update(linear_system_output_dir="invalid_low")
        self.cfg["solver"]["split_contact_bias_low_V"] = {"drain": .1}
        q = self.run_probe()
        self.assertNotEqual(q.returncode, 0)
        self.assertIn("normalized split", q.stderr)

    def test_nonzero_supply_or_other_terminal_is_rejected(self):
        self.cfg["supply_voltage_V"] = .1
        q = self.run_probe()
        self.assertNotEqual(q.returncode, 0)
        self.cfg["supply_voltage_V"] = 0
        self.cfg["contacts"][0]["bias"] = .1
        q = self.run_probe()
        self.assertNotEqual(q.returncode, 0)
        self.assertIn("other terminals", q.stderr)

    def test_bounded_supply_enters_only_the_circuit_rhs(self):
        self.cfg["simulation_type"] = "newton_bounded_supply_load_line_iteration_probe"
        snapshots = []
        for index, supply in enumerate([0., .001]):
            self.cfg["supply_voltage_V"] = supply
            folder = self.root / f"bounded_{index}"
            self.cfg["linear_system_output_dir"] = str(folder)
            q = self.run_probe()
            self.assertEqual(q.returncode, 0, q.stderr)
            report = json.loads((folder / "summary.json").read_text())
            self.assertFalse(report["qualified"])
            self.assertEqual(report["new_sweep_points"], 0)
            self.assertLessEqual(report["accepted_updates"], 1)
            import csv
            with (folder / "vectors.csv").open() as stream:
                vectors = list(csv.DictReader(stream))
            matrix = (folder / "physical_bordered.mtx").read_text()
            lines = [line for line in matrix.splitlines()
                     if line and not line.startswith("%")]
            size = int(lines[0].split()[0])
            diagonal = next(float(value) for i, j, value in
                            (line.split() for line in lines[1:])
                            if int(i) == size and int(j) == size)
            snapshots.append((matrix, vectors, diagonal, report))
        self.assertEqual(snapshots[0][0], snapshots[1][0])
        for left, right in zip(snapshots[0][1][:-1], snapshots[1][1][:-1]):
            self.assertEqual(left["physical_rhs"], right["physical_rhs"])
        rhs_change = (float(snapshots[1][1][-1]["physical_rhs"]) -
                      float(snapshots[0][1][-1]["physical_rhs"]))
        self.assertAlmostEqual(rhs_change / snapshots[1][2], .001, places=14)
        self.assertAlmostEqual(snapshots[1][3]["initial_load_line_V"] -
                               snapshots[0][3]["initial_load_line_V"], -.001, places=14)

    def test_bounded_supply_rejects_unbounded_or_invalid_requests(self):
        self.cfg["simulation_type"] = "newton_bounded_supply_load_line_iteration_probe"
        for index, updates in enumerate([{"supply_voltage_V": .201},
                                         {"max_abs_device_voltage_V": .201},
                                         {"resistance_ohm_um": -1}]):
            original = self.cfg.copy()
            self.cfg.update(updates)
            self.cfg["linear_system_output_dir"] = str(self.root / f"invalid_{index}")
            q = self.run_probe()
            self.assertNotEqual(q.returncode, 0)
            self.assertIn("positive resistance and a voltage box", q.stderr)
            self.cfg = original

    def test_projection_policy_requires_split_and_rejects_unknown(self):
        self.cfg["simulation_type"] = "newton_bounded_supply_load_line_iteration_probe"
        for index, policy in enumerate(["rounded_voltage_projection", "unknown"]):
            self.cfg["trial_boundary_policy"] = policy
            self.cfg["linear_system_output_dir"] = str(self.root / f"policy_{index}")
            q = self.run_probe()
            self.assertNotEqual(q.returncode, 0)
            self.assertIn("requires split state" if index == 0 else "invalid trial boundary policy", q.stderr)

    def test_qualification_never_performs_a_newton_update(self):
        self.cfg["simulation_type"] = "newton_zero_update_qualification_probe"
        self.cfg["doping"][0]["donors"] = 1e17
        self.cfg["solver"]["carrier_row_convergence"] = {
            "mode": "enforce", "eps_row": .001, "scale_floor": 1e-30,
            "min_carrier_density_m3": 1e-200}
        q = self.run_probe()
        # The probe can report a failed qualification, but must execute no step.
        records = [json.loads(l) for l in q.stdout.splitlines() if l.startswith("{")]
        self.assertTrue(records, q.stderr)
        self.assertEqual(records[-1]["iterations"], 0)
        self.assertFalse(records[-1]["converged"])
        rows = records[-1]["carrier_row_convergence"]
        groups = rows["scale_groups"]
        self.assertTrue(groups["diagnostic_only"])
        self.assertEqual(sum(groups[g]["row_count"] for g in
                             ("floor_limited", "non_floor_limited")), rows["qualified_row_count"])
        self.assertEqual(sum(groups[g]["violation_count"] for g in
                             ("floor_limited", "non_floor_limited")), rows["violation_count"])
        self.assertAlmostEqual(rows["floor_absolute_limit"] / 1e-33, 1)
        for row in rows["violations"]:
            self.assertEqual(row["floor_limited"], row["physical_scale"] <= rows["scale_floor"])
        self.cfg["solver"]["carrier_row_convergence"]["floor_limited_action"] = "report"
        q = self.run_probe()
        reported = [json.loads(l) for l in q.stdout.splitlines() if l.startswith("{")][-1]["carrier_row_convergence"]
        self.assertEqual(reported["floor_limited_action"], "report")
        self.assertEqual(reported["blocking_violation_count"], reported["scale_groups"]["non_floor_limited"]["violation_count"])
        self.assertEqual(reported["satisfied"], reported["blocking_violation_count"] == 0)
        for row in reported["violations"]:
            self.assertEqual(row["blocks_acceptance"], not row["floor_limited"])
        self.cfg["solver"]["carrier_row_convergence"]["floor_limited_action"] = "invalid"
        self.assertNotEqual(self.run_probe().returncode, 0)

    def test_split_contact_audit_is_read_only(self):
        self.cfg["simulation_type"] = "newton_candidate_conservation_probe"
        self.cfg["candidate_port_diagnostics"] = True
        self.cfg["solver"]["split_dd_state"] = True
        self.cfg["solver"]["mobility"] = {
            "model": "phumob", "edge_averaging": "element_box_phumob",
            "doping_concentration_basis": "total_impurity"}
        self.cfg["mesh_geometry"] = {"poisson_permittivity_policy": "cell_material"}
        self.cfg["solver"]["region_resolved_interface_assembly"] = {
            "transport_edge_geometry": "element_box"}
        self.cfg["solver"]["recombination"] = ["srh"]
        self.cfg["solver"]["global_continuity_closure"] = {"mode": "enforce"}
        q = self.run_probe()
        self.assertEqual(q.returncode, 0, q.stderr)
        result = json.loads((self.root / "rows.csv").read_text())
        self.assertFalse(result["device_qualified"])
        self.assertEqual(set(result["contact_residual_functional_A_per_um"]), {"source", "drain"})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", required=True, type=Path)
    args, rest = parser.parse_known_args()
    RUNNER = args.runner.resolve()
    unittest.main(argv=[sys.argv[0]] + rest)
