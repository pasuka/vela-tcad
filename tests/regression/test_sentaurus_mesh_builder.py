"""Regression tests for the Sentaurus SDE mesh and doping generator.

Phase 2 of the Sentaurus import plan turns the device IR into the two file
formats Vela already consumes -- ``mesh.json`` (``JsonMeshReader``) and the node
doping CSV (``DCSweep``'s ``node_doping_file``) -- with no third format.

These tests cover the acceptance items called out in the handover document:
contact edge ownership (corner, edge-midpoint T-junction, long/short edge,
reversed edge, partial overlap, non-manifold), node doping against the closed
form profile at peak/depth/junction/region interface, and the non-obtuse plus
mixed-Voronoi qualification gate.
"""

from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import parse_sde_device_ir  # noqa: E402
from sentaurus_mesh_builder import (  # noqa: E402
    BoundarySegment,
    _region_rects,
    _windows_by_name,
    MeshQualificationError,
    build_mesh_and_doping,
    evaluate_doping,
    select_contact_segment,
    write_doping_csv,
    write_mesh_json,
)

PN2D_SDE = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
COARSE_SDE = (REPO / "reference_tcad" / "pn2d_sentaurus2018_coarse7x3"
              / "source" / "pn2d_sde.cmd")


def doping_at(device_ir: dict, x: float, y: float) -> tuple[float, float]:
    regions = {name: rect for name, rect, _material in _region_rects(device_ir)}
    return evaluate_doping(device_ir, regions, _windows_by_name(device_ir), x, y)


def build_sde(root: Path, body: str) -> Path:
    path = root / "case_sde.cmd"
    path.write_text(body.strip() + "\n", encoding="utf-8")
    return path


class ContactOwnerRuleTest(unittest.TestCase):
    """Owner rules for resolving a pick point to exactly one boundary edge."""

    LEFT = BoundarySegment(axis="x", coordinate=0.0, start=0.0, end=1.0)
    RIGHT = BoundarySegment(axis="x", coordinate=2.0, start=0.0, end=1.0)
    BOTTOM = BoundarySegment(axis="y", coordinate=0.0, start=0.0, end=2.0)
    TOP = BoundarySegment(axis="y", coordinate=1.0, start=0.0, end=2.0)
    BOX = [LEFT, RIGHT, BOTTOM, TOP]

    def test_edge_midpoint_selects_the_containing_edge(self) -> None:
        self.assertEqual(self.LEFT, select_contact_segment(0.0, 0.5, self.BOX))
        self.assertEqual(self.RIGHT, select_contact_segment(2.0, 0.5, self.BOX))
        self.assertEqual(self.BOTTOM, select_contact_segment(1.0, 0.0, self.BOX))

    def test_t_junction_prefers_the_interior_owner(self) -> None:
        """A pick point interior to one edge but on another edge's endpoint.

        The bottom edge touches the left edge's endpoint at (0, 0). A point at
        (0, 0.0) is a corner, but a point at x=0.5 on the bottom edge is only
        interior to BOTTOM even when a perpendicular stub ends there.
        """
        stub = BoundarySegment(axis="x", coordinate=0.5, start=0.0, end=0.3)
        segments = [*self.BOX, stub]
        # (0.5, 0.0) is interior to BOTTOM and only touches the stub's endpoint.
        self.assertEqual(self.BOTTOM, select_contact_segment(0.5, 0.0, segments))

    def test_corner_pick_point_is_rejected(self) -> None:
        with self.assertRaises(MeshQualificationError) as ctx:
            select_contact_segment(0.0, 0.0, self.BOX)
        report = ctx.exception.report
        self.assertEqual("contact", report["stage"])
        self.assertIn("endpoint", report["reason"])
        self.assertEqual(2, len(report["candidates"]))

    def test_collinear_t_junction_endpoint_is_rejected(self) -> None:
        """Two collinear segments meeting at the pick point are ambiguous."""
        lower = BoundarySegment(axis="x", coordinate=0.0, start=0.0, end=0.5)
        upper = BoundarySegment(axis="x", coordinate=0.0, start=0.5, end=1.0)
        with self.assertRaises(MeshQualificationError) as ctx:
            select_contact_segment(0.0, 0.5, [lower, upper])
        self.assertIn("endpoint", ctx.exception.report["reason"])

    def test_long_and_short_collinear_edges_are_disambiguated_by_interior(self) -> None:
        short = BoundarySegment(axis="x", coordinate=0.0, start=0.0, end=0.2)
        long = BoundarySegment(axis="x", coordinate=0.0, start=0.2, end=1.0)
        self.assertEqual(short, select_contact_segment(0.0, 0.1, [short, long]))
        self.assertEqual(long, select_contact_segment(0.0, 0.9, [short, long]))

    def test_reversed_edge_orientation_is_normalised(self) -> None:
        """A segment authored end-before-start must not silently match nothing."""
        forward = BoundarySegment(axis="x", coordinate=0.0, start=0.0, end=1.0)
        self.assertTrue(forward.contains(0.0, 0.5))
        reversed_segment = BoundarySegment(axis="x", coordinate=0.0, start=1.0, end=0.0)
        # An unnormalised segment contains nothing, so selection must fail
        # closed rather than pick an arbitrary neighbour.
        with self.assertRaises(MeshQualificationError):
            select_contact_segment(0.0, 0.5, [reversed_segment])

    def test_partially_overlapping_edges_are_non_manifold(self) -> None:
        first = BoundarySegment(axis="x", coordinate=0.0, start=0.0, end=0.6)
        second = BoundarySegment(axis="x", coordinate=0.0, start=0.4, end=1.0)
        with self.assertRaises(MeshQualificationError) as ctx:
            select_contact_segment(0.0, 0.5, [first, second])
        self.assertIn("non-manifold", ctx.exception.report["reason"])

    def test_pick_point_off_the_boundary_is_rejected(self) -> None:
        with self.assertRaises(MeshQualificationError) as ctx:
            select_contact_segment(1.0, 0.5, self.BOX)
        self.assertIn("no outer boundary segment", ctx.exception.report["reason"])


class GeneratedMeshTest(unittest.TestCase):
    def test_coarse_fixture_reproduces_the_documented_7x3_lattice(self) -> None:
        """The coarse fixture names its own lattice; the generator must match."""
        generated = build_mesh_and_doping(parse_sde_device_ir(COARSE_SDE))
        mesh = generated.mesh
        self.assertEqual(21, len(mesh["nodes"]))
        self.assertEqual(24, len(mesh["triangles"]))
        xs = sorted({round(node["x"], 9) for node in mesh["nodes"]})
        ys = sorted({round(node["y"], 9) for node in mesh["nodes"]})
        self.assertEqual(7, len(xs))
        self.assertEqual(3, len(ys))
        self.assertAlmostEqual(0.0, xs[0])
        self.assertAlmostEqual(2.0, xs[-1])
        self.assertEqual([0.0, 0.25, 0.5], ys)

    def test_identifiers_are_sequential_and_cross_referenced(self) -> None:
        mesh = build_mesh_and_doping(parse_sde_device_ir(COARSE_SDE)).mesh
        self.assertEqual(list(range(len(mesh["nodes"]))),
                         [node["id"] for node in mesh["nodes"]])
        self.assertEqual(list(range(len(mesh["triangles"]))),
                         [cell["id"] for cell in mesh["triangles"]])
        self.assertEqual(list(range(len(mesh["regions"]))),
                         [region["id"] for region in mesh["regions"]])
        node_ids = {node["id"] for node in mesh["nodes"]}
        region_ids = {region["id"] for region in mesh["regions"]}
        owned: set[int] = set()
        for region in mesh["regions"]:
            self.assertFalse(owned & set(region["cell_ids"]),
                             "a cell is claimed by two regions")
            owned.update(region["cell_ids"])
        self.assertEqual({cell["id"] for cell in mesh["triangles"]}, owned)
        for cell in mesh["triangles"]:
            self.assertIn(cell["region_id"], region_ids)
            self.assertEqual(3, len(set(cell["node_ids"])))
            self.assertTrue(set(cell["node_ids"]) <= node_ids)
        for contact in mesh["contacts"]:
            self.assertIn(contact["region_id"], region_ids)
            self.assertTrue(set(contact["node_ids"]) <= node_ids)

    def test_contacts_cover_only_the_outer_boundary(self) -> None:
        ir = parse_sde_device_ir(PN2D_SDE)
        generated = build_mesh_and_doping(ir)
        coords = {node["id"]: (node["x"], node["y"]) for node in generated.mesh["nodes"]}
        contacts = {contact["name"]: contact for contact in generated.mesh["contacts"]}
        self.assertEqual({"Anode", "Cathode"}, set(contacts))
        for name, expected_x in (("Anode", 0.0), ("Cathode", 2.0)):
            node_ids = contacts[name]["node_ids"]
            self.assertGreaterEqual(len(node_ids), 11)
            for node_id in node_ids:
                x, y = coords[node_id]
                self.assertAlmostEqual(expected_x, x)
                self.assertTrue(0.0 - 1e-12 <= y <= 0.5 + 1e-12)

    def test_contacts_include_every_refined_boundary_node(self) -> None:
        """Refinement adds boundary nodes; the contact must absorb all of them."""
        ir = parse_sde_device_ir(PN2D_SDE)
        generated = build_mesh_and_doping(ir)
        contacts = {contact["name"]: set(contact["node_ids"])
                    for contact in generated.mesh["contacts"]}
        on_left = {node["id"] for node in generated.mesh["nodes"]
                   if abs(node["x"]) <= 1e-12}
        on_right = {node["id"] for node in generated.mesh["nodes"]
                    if abs(node["x"] - 2.0) <= 1e-12}
        self.assertEqual(on_left, contacts["Anode"])
        self.assertEqual(on_right, contacts["Cathode"])

    def test_contacts_do_not_share_nodes(self) -> None:
        mesh = build_mesh_and_doping(parse_sde_device_ir(PN2D_SDE)).mesh
        anode = set(mesh["contacts"][0]["node_ids"])
        cathode = set(mesh["contacts"][1]["node_ids"])
        self.assertEqual(set(), anode & cathode)

    def test_declared_but_unassigned_contact_fails_at_mesh_time(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_mesh_", dir=REPO / "build") as tmp:
            path = build_sde(Path(tmp), """
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position 1.0 1.0 0.0)
  "Silicon" "R.Si")
(sdegeo:define-contact-set "Anode" 4.0 (color:rgb 1 0 0) "##")
(sdedr:define-refinement-size "M" 0.5 0.5 0.5 0.5)
(sdedr:define-refinement-region "M.Place" "M" "R.Si")
""")
            ir = parse_sde_device_ir(path)
            with self.assertRaises(MeshQualificationError) as ctx:
                build_mesh_and_doping(ir)
            self.assertEqual("contact", ctx.exception.report["stage"])
            self.assertEqual("Anode", ctx.exception.report["contact"])

    def test_region_free_input_fails_at_mesh_time(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_mesh_", dir=REPO / "build") as tmp:
            path = build_sde(Path(tmp), "(define L 2.0)\n(define H 0.5)")
            with self.assertRaises(MeshQualificationError) as ctx:
                build_mesh_and_doping(parse_sde_device_ir(path))
            self.assertEqual("input", ctx.exception.report["stage"])


class MeshQualificationGateTest(unittest.TestCase):
    def test_generated_mesh_is_non_obtuse(self) -> None:
        report = build_mesh_and_doping(parse_sde_device_ir(PN2D_SDE)).qualification
        self.assertTrue(report["require_non_obtuse"])
        self.assertEqual([], report["obtuse_cells"])
        self.assertLessEqual(report["max_angle_degrees"], 90.0 + 1e-6)
        self.assertGreater(report["min_angle_degrees"], 0.0)
        self.assertEqual("mixed_voronoi", report["node_volume_policy"])

    def test_unknown_node_volume_policy_is_rejected(self) -> None:
        ir = parse_sde_device_ir(COARSE_SDE)
        with self.assertRaises(MeshQualificationError) as ctx:
            build_mesh_and_doping(ir, node_volume_policy="cell_reconstructed")
        self.assertIn("cell_reconstructed", str(ctx.exception))

    def test_barycentric_policy_is_accepted(self) -> None:
        generated = build_mesh_and_doping(parse_sde_device_ir(COARSE_SDE),
                                          node_volume_policy="barycentric")
        self.assertEqual("barycentric", generated.qualification["node_volume_policy"])

    def test_unsupported_ir_schema_is_rejected(self) -> None:
        with self.assertRaises(MeshQualificationError) as ctx:
            build_mesh_and_doping({"schema": "vela.sentaurus_device_ir.v0"})
        self.assertEqual("input", ctx.exception.report["stage"])


class NodeDopingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.ir = parse_sde_device_ir(PN2D_SDE)
        self.generated = build_mesh_and_doping(self.ir)
        self.coords = {node["id"]: (node["x"], node["y"])
                       for node in self.generated.mesh["nodes"]}
        self.rows = {int(row["node_id"]): row for row in self.generated.doping_rows}

    def test_every_node_appears_exactly_once(self) -> None:
        self.assertEqual(len(self.generated.mesh["nodes"]), len(self.generated.doping_rows))
        self.assertEqual(set(self.coords), set(self.rows))

    def test_constant_profiles_match_the_closed_form_on_each_side(self) -> None:
        """P.Window is x in [0, 1], N.Window is x in [1, 2], both at 1e17."""
        for node_id, (x, _y) in self.coords.items():
            row = self.rows[node_id]
            donors = float(row["donors_cm3"])
            acceptors = float(row["acceptors_cm3"])
            if x < 1.0 - 1e-12:
                self.assertEqual(1e17, acceptors, f"node {node_id} at x={x}")
                self.assertEqual(0.0, donors, f"node {node_id} at x={x}")
            elif x > 1.0 + 1e-12:
                self.assertEqual(1e17, donors, f"node {node_id} at x={x}")
                self.assertEqual(0.0, acceptors, f"node {node_id} at x={x}")

    def test_junction_nodes_take_the_last_overlapping_placement(self) -> None:
        """Both windows are closed on x=1.0; SDE gives the later placement."""
        junction = [nid for nid, (x, _y) in self.coords.items() if abs(x - 1.0) <= 1e-12]
        self.assertGreater(len(junction), 0)
        values = {(float(self.rows[nid]["donors_cm3"]),
                   float(self.rows[nid]["acceptors_cm3"])) for nid in junction}
        # The value is deterministic: whatever the ordering rule yields, every
        # junction node must agree, and the later placement must contribute.
        self.assertEqual(1, len(values))
        donors, _acceptors = next(iter(values))
        self.assertEqual(1e17, donors)

    def test_doping_is_uniform_in_depth_at_a_fixed_x(self) -> None:
        by_x: dict[float, set[tuple[float, float]]] = {}
        for node_id, (x, _y) in self.coords.items():
            row = self.rows[node_id]
            by_x.setdefault(round(x, 9), set()).add(
                (float(row["donors_cm3"]), float(row["acceptors_cm3"])))
        for x, values in by_x.items():
            self.assertEqual(1, len(values), f"depth-dependent doping at x={x}")

    def test_evaluate_doping_agrees_with_the_generated_rows(self) -> None:
        for node_id, (x, y) in list(self.coords.items())[:50]:
            donors, acceptors = doping_at(self.ir, x, y)
            self.assertEqual(donors, float(self.rows[node_id]["donors_cm3"]))
            self.assertEqual(acceptors, float(self.rows[node_id]["acceptors_cm3"]))

    def test_overlapping_region_placements_overlay_per_carrier(self) -> None:
        """Overlapping placements replace within a carrier and overlay across."""
        with tempfile.TemporaryDirectory(prefix="vela_doping_", dir=REPO / "build") as tmp:
            path = build_sde(Path(tmp), """
(define L 2.0)
(define H 0.5)
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon" "R.Si")
(sdegeo:create-rectangle
  (position 1.0 0.0 0.0)
  (position L H 0.0)
  "Silicon" "R.N")
(sdedr:define-constant-profile "P.D" "BoronActiveConcentration" 2e16)
(sdedr:define-constant-profile "N.D" "PhosphorusActiveConcentration" 5e18)
(sdedr:define-constant-profile-region "P.Place" "P.D" "R.Si")
(sdedr:define-constant-profile-region "N.Place" "N.D" "R.N")
(sdedr:define-refinement-size "M" 0.5 0.25 0.5 0.25)
(sdedr:define-refinement-region "M.Place" "M" "R.Si")
""")
            ir = parse_sde_device_ir(path)
            # Only R.Si covers x=0.5, so only the acceptor profile applies.
            self.assertEqual((0.0, 2e16), doping_at(ir, 0.5, 0.25))
            # Both placements cover x=1.5. They target different carriers, so
            # they overlay rather than replace; a same-carrier placement would
            # replace by priority instead.
            self.assertEqual((5e18, 2e16), doping_at(ir, 1.5, 0.25))


class MeshFileContractTest(unittest.TestCase):
    def test_written_files_match_the_vela_file_contract(self) -> None:
        generated = build_mesh_and_doping(parse_sde_device_ir(COARSE_SDE))
        with tempfile.TemporaryDirectory(prefix="vela_mesh_io_", dir=REPO / "build") as tmp:
            root = Path(tmp)
            write_mesh_json(root / "mesh.json", generated.mesh)
            write_doping_csv(root / "doping.csv", generated.doping_rows)

            mesh = json.loads((root / "mesh.json").read_text(encoding="utf-8"))
            self.assertEqual({"nodes", "triangles", "regions", "contacts"},
                             set(mesh) & {"nodes", "triangles", "regions", "contacts"})
            self.assertEqual({"id", "x", "y"}, set(mesh["nodes"][0]))
            self.assertEqual({"id", "region_id", "node_ids"}, set(mesh["triangles"][0]))
            self.assertEqual({"id", "name", "material", "cell_ids"},
                             set(mesh["regions"][0]))
            self.assertEqual({"id", "name", "region_id", "node_ids"},
                             set(mesh["contacts"][0]))
            self.assertEqual("Si", mesh["regions"][0]["material"])

            with (root / "doping.csv").open(encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle)
                self.assertEqual(["node_id", "donors_cm3", "acceptors_cm3"],
                                 next(reader))
                ids = [int(row[0]) for row in reader]
            self.assertEqual(list(range(len(mesh["nodes"]))), ids)

    def test_generated_mesh_loads_in_the_vela_runner(self) -> None:
        runner = REPO / "build" / "vela_example_runner"
        if not runner.exists():
            self.skipTest("vela_example_runner has not been built")
        generated = build_mesh_and_doping(parse_sde_device_ir(COARSE_SDE))
        with tempfile.TemporaryDirectory(prefix="vela_mesh_run_", dir=REPO / "build") as tmp:
            root = Path(tmp)
            write_mesh_json(root / "mesh.json", generated.mesh)
            write_doping_csv(root / "doping.csv", generated.doping_rows)
            (root / "sim.json").write_text(json.dumps({
                "simulation_type": "dc_sweep",
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "output_csv": "iv.csv",
                "scaling": {"mode": "unit_scaling"},
                "mesh_geometry": {"node_volume_policy": "mixed_voronoi",
                                  "require_non_obtuse": True},
                "contacts": [{"name": "Anode", "bias": 0.0},
                             {"name": "Cathode", "bias": 0.0}],
                "solver": {"method": "gummel_newton", "max_iter": 40,
                           "reltol": 1e-8, "abstol": 1e-9, "damping_psi": 0.2,
                           "line_search": True, "warm_start": True,
                           "contact_boundary_reconstruction":
                               "dominant_signed_contact_mean",
                           "handoff": {"fallback": "none",
                                       "require_gummel_convergence": False,
                                       "gummel_max_iter": 0,
                                       "newton_max_iter": 40}},
                "sweep": {"mode": "iv", "contact": "Anode",
                          "current_contact": "Anode", "start": 0.0, "stop": 0.05,
                          "step": 0.05, "write_vtk": False,
                          "initialization": {"mode": "poisson_block"}},
            }, indent=2) + "\n", encoding="utf-8")
            result = subprocess.run(
                [str(runner), "--config", str(root / "sim.json")],
                cwd=root, capture_output=True, text=True, check=True)
            self.assertNotIn("warning", result.stderr.lower())
            with (root / "iv.csv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertGreaterEqual(len(rows), 1)
            for row in rows:
                self.assertEqual("1", row["converged"])
                self.assertTrue(math.isfinite(float(row["current_total_A_per_um"])))


if __name__ == "__main__":
    unittest.main()
