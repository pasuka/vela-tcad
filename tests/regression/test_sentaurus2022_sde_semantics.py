#!/usr/bin/env python3
"""Golden tests for the Sentaurus Structure Editor T-2022.03 semantics."""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_device_ir import SdeParseError, parse_sde_device_ir  # noqa: E402
from sentaurus_mesh_builder import (  # noqa: E402
    _region_rects,
    _owner_at,
    _windows_by_name,
    build_mesh_and_doping,
    evaluate_doping,
)


def parse_text(text: str) -> dict:
    root = tempfile.TemporaryDirectory(prefix="vela_sde2022_")
    try:
        path = Path(root.name) / "case_sde.cmd"
        path.write_text(text.strip() + "\n", encoding="utf-8")
        return parse_sde_device_ir(path)
    finally:
        root.cleanup()


RECT = """
(sdegeo:create-rectangle (position 0 0 0) (position 2 1 0) "Silicon" "R.Base")
"""


def doping_at(ir: dict, x: float, y: float) -> tuple[float, float]:
    regions = {name: rect for name, rect, _material in _region_rects(ir)}
    return evaluate_doping(ir, regions, _windows_by_name(ir), x, y)


class BooleanSemanticsTest(unittest.TestCase):
    def test_aba_new_replaces_old_and_bab_old_replaces_new(self) -> None:
        common = """
(sdegeo:create-rectangle (position 0 0 0) (position 2 1 0) "Silicon" "R.Old")
(sdegeo:set-default-boolean "%s")
(sdegeo:create-rectangle (position 1 0 0) (position 3 1 0) "Silicon" "R.New")
(sdegeo:define-contact-set "left")
(sdegeo:set-contact (find-edge-id (position 0 0.5 0)) "left")
(sdegeo:define-contact-set "right")
(sdegeo:set-contact (find-edge-id (position 3 0.5 0)) "right")
"""
        aba = build_mesh_and_doping(parse_text(common % "ABA")).mesh
        bab = build_mesh_and_doping(parse_text(common % "BAB")).mesh
        def owner(mesh: dict, x: float) -> str:
            coords = {node["id"]: node for node in mesh["nodes"]}
            regions = {region["id"]: region["name"] for region in mesh["regions"]}
            for cell in mesh["triangles"]:
                cx = sum(coords[n]["x"] for n in cell["node_ids"]) / 3.0
                if abs(cx - x) < 0.34:
                    return regions[cell["region_id"]]
            raise AssertionError(x)
        self.assertEqual("R.New", owner(aba, 1.5))
        self.assertEqual("R.Old", owner(bab, 1.5))

    def test_unimplemented_boolean_modes_fail_closed(self) -> None:
        with self.assertRaises(SdeParseError):
            parse_text('(sdegeo:set-default-boolean "ABiA")')

    def test_offline_control_flow_reports_evaluated_frontend_hint(self) -> None:
        with self.assertRaises(SdeParseError) as ctx:
            parse_text('(set! L 2)')
        self.assertIn("sentaurus2022-evaluated", str(ctx.exception))

    def test_complete_containment_and_multiple_old_regions(self) -> None:
        aba = parse_text("""
(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R.Left")
(sdegeo:create-rectangle (position 1 0 0) (position 2 1 0) "Silicon" "R.Right")
(sdegeo:set-default-boolean "ABA")
(sdegeo:create-rectangle (position 0.5 0 0) (position 1.5 1 0) "Silicon" "R.New")
""")
        rects = {name: rect for name, rect, _material in _region_rects(aba)}
        self.assertEqual("R.Left", _owner_at(aba, rects, 0.25, 0.5))
        self.assertEqual("R.New", _owner_at(aba, rects, 0.75, 0.5))
        self.assertEqual("R.New", _owner_at(aba, rects, 1.25, 0.5))
        self.assertEqual("R.Right", _owner_at(aba, rects, 1.75, 0.5))

        bab = parse_text("""
(sdegeo:create-rectangle (position 0 0 0) (position 2 1 0) "Silicon" "R.Old")
(sdegeo:set-default-boolean "BAB")
(sdegeo:create-rectangle (position 0.5 0.2 0) (position 1.5 0.8 0) "Silicon" "R.Contained")
""")
        rects = {name: rect for name, rect, _material in _region_rects(bab)}
        self.assertEqual("R.Old", _owner_at(bab, rects, 1.0, 0.5))


class ContactSemanticsTest(unittest.TestCase):
    def test_insert_vertex_limits_contact_to_complete_split_edge(self) -> None:
        ir = parse_text(RECT + """
(sdegeo:insert-vertex (position 0 0.5 0))
(sdegeo:set-contact (find-edge-id (position 0 0.25 0)) "lower-left")
(sdegeo:set-contact (find-edge-id (position 2 0.5 0)) "right")
""")
        mesh = build_mesh_and_doping(ir).mesh
        coords = {node["id"]: (node["x"], node["y"]) for node in mesh["nodes"]}
        contact = next(item for item in mesh["contacts"] if item["name"] == "lower-left")
        ys = [coords[node][1] for node in contact["node_ids"]]
        self.assertAlmostEqual(0.0, min(ys))
        self.assertAlmostEqual(0.5, max(ys))


class ReplacementAndGaussianTest(unittest.TestCase):
    def test_constant_noreplace_is_omitted_and_replacement_requires_decay_slot(self) -> None:
        with self.assertRaises(SdeParseError):
            parse_text(RECT + """
(sdedr:define-constant-profile "D" "PhosphorusActiveConcentration" 1e16)
(sdedr:define-constant-profile-region "P" "D" "R.Base" "NoReplace")
""")
        with self.assertRaises(SdeParseError):
            parse_text(RECT + """
(sdedr:define-constant-profile "D" "PhosphorusActiveConcentration" 1e16)
(sdedr:define-constant-profile-region "P" "D" "R.Base" "Replace")
""")

    def test_three_replacement_modes_are_species_aware(self) -> None:
        base = RECT + """
(sdedr:define-constant-profile "D1" "PhosphorusActiveConcentration" 1e16)
(sdedr:define-constant-profile "D2" "PhosphorusActiveConcentration" 2e16)
(sdedr:define-constant-profile "A" "BoronActiveConcentration" 3e16)
(sdedr:define-constant-profile-region "P1" "D1" "R.Base")
(sdedr:define-constant-profile-region "P2" "D2" "R.Base")
"""
        self.assertEqual((3e16, 0.0), doping_at(parse_text(base), 1, 0.5))
        local = base.replace('"P2" "D2" "R.Base"',
                             '"P2" "D2" "R.Base" 0 "LocalReplace"')
        self.assertEqual((2e16, 0.0), doping_at(parse_text(local), 1, 0.5))
        replace = base + '(sdedr:define-constant-profile-region "P3" "A" "R.Base" 0 "Replace")\n'
        self.assertEqual((0.0, 3e16), doping_at(parse_text(replace), 1, 0.5))

    def test_gaussian_positive_is_left_normal_and_noeval_excludes_baseline(self) -> None:
        ir = parse_text("""
(sdegeo:create-rectangle (position 0 -1 0) (position 1 1 0) "Silicon" "R.Base")
(sdedr:define-refeval-window "B" "Line" (position 0 0 0) (position 1 0 0))
(sdedr:define-gaussian-profile "G" "PhosphorusActiveConcentration"
 "PeakPos" 0.2 "PeakVal" 1e18 "ValueAtDepth" 1e16 "Depth" 0.4 "Gauss" "Factor" 1)
(sdedr:define-analytical-profile-placement "P" "G" "B" "Positive" "NoReplace" "NoEval")
""")
        self.assertTrue(math.isclose(doping_at(ir, 0.5, 0.2)[0], 1e18, rel_tol=1e-12))
        self.assertEqual(0.0, doping_at(ir, 0.5, -0.2)[0])
        self.assertEqual(0.0, doping_at(ir, 0.5, 0.0)[0])


if __name__ == "__main__":
    unittest.main()
