"""WP4a: pin the Sentaurus gate-contact (Eq. 104) -> Vela ``metal_gate`` mapping.

Formula tests run everywhere; the flat-band fixture needs ``vela_example_runner``
and is skipped when the binary is absent.  Evidence and derivation:
``docs/validation/sdevice_gate_contact_reference_2026-09-29.md``.
"""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import check_sdevice_gate_contact_mapping as gate  # noqa: E402


def _tmpdir() -> tempfile.TemporaryDirectory:
    base = REPO / "build"
    return tempfile.TemporaryDirectory(dir=base if base.is_dir() else None, prefix="sdevice_gate_contact_")


class MappingFormulaTest(unittest.TestCase):
    def test_intrinsic_work_function_uses_vela_silicon(self) -> None:
        expected = 4.05 + 0.56 + 0.5 * gate.VT * math.log(2.8e19 / 1.04e19)
        self.assertAlmostEqual(expected, gate.intrinsic_work_function_eV(), delta=1e-15)
        self.assertAlmostEqual(4.6228018935, gate.intrinsic_work_function_eV(), delta=1e-9)

    def test_barrier_maps_one_to_one_and_workfunction_through_phi_s_i(self) -> None:
        self.assertAlmostEqual(0.55, gate.phi_ms_from_workfunction(0.55 + gate.intrinsic_work_function_eV()),
                               delta=1e-15)
        self.assertAlmostEqual(5.1 - 4.6228018935, gate.phi_ms_from_workfunction(5.1), delta=1e-9)
        # A metal whose work function equals the electron affinity sits Eg/2 (+ DOS term) above E_i.
        self.assertAlmostEqual(-(0.56 + 0.5 * gate.VT * math.log(2.8e19 / 1.04e19)),
                               gate.phi_ms_from_workfunction(4.05), delta=1e-15)

    def test_flatband_prediction_matches_textbook_limit(self) -> None:
        # p-type: V_FB = Phi_MS - phi_Fp with phi_Fp = Vt ln(NA/ni) (asinh -> ln for NA >> ni).
        phi_fp = gate.VT * math.log(1.0e17 / 1.0e10)
        self.assertAlmostEqual(0.55 - phi_fp, gate.predicted_flatband_V(0.55, 0.0, 1.0e17), delta=1e-9)
        self.assertAlmostEqual(0.55 + phi_fp, gate.predicted_flatband_V(0.55, 1.0e17, 0.0), delta=1e-9)
        self.assertEqual(0.55, gate.predicted_flatband_V(0.55, 0.0, 0.0))

    def test_fixture_mesh_is_a_two_region_strip_with_probes(self) -> None:
        fixture = gate.build_mesh()
        mesh = fixture["mesh"]
        self.assertEqual(["R.Si", "R.Ox"], [r["name"] for r in mesh["regions"]])
        self.assertEqual(["Si", "SiO2"], [r["material"] for r in mesh["regions"]])
        self.assertEqual({"substrate", "gate"}, {c["name"] for c in mesh["contacts"]})
        nodes = {n["id"]: n for n in mesh["nodes"]}
        for node_id in fixture["probes"]["surface"]:
            self.assertAlmostEqual(1.0, nodes[node_id]["y"], delta=1e-12)
        for node_id in fixture["probes"]["gate"]:
            self.assertAlmostEqual(1.01, nodes[node_id]["y"], delta=1e-12)
        self.assertEqual(len(mesh["triangles"]), 2 * (len(mesh["nodes"]) // 2 - 1))


class FlatbandFixtureTest(unittest.TestCase):
    def test_flatband_matches_eq104_within_1mV_for_three_barriers(self) -> None:
        runner = gate.find_runner()
        if runner is None:
            self.skipTest("vela_example_runner not built; gate-contact flat-band fixture skipped")
        fixture = gate.build_mesh()
        with _tmpdir() as tmp:
            workdir = Path(tmp)
            (workdir / "mesh.json").write_text(json.dumps(fixture["mesh"]), encoding="utf-8")
            deviations = {}
            for phi_ms in (0.0, 0.55, -0.3):
                result = gate.check_barrier(runner, workdir, phi_ms, fixture["probes"], 1.0e17, 0.0, refine=1)
                deviations[phi_ms] = result["deviation_V"]
                self.assertLessEqual(abs(result["deviation_V"]), 1.0e-3,
                                     f"Phi_MS={phi_ms}: V_FB deviates {result['deviation_V']:.3e} V from Eq. 104")
                # Successive refinements agree, i.e. the extraction is converged well below 1 mV.
                self.assertLess(abs(result["estimates"][-1] - result["estimates"][0]), 1.0e-4)
            # The mapping is linear in Phi_MS: flat band shifts one-to-one with the barrier.
            self.assertAlmostEqual(0.55, deviations[0.55] - deviations[0.0] + 0.55, delta=1.0e-3)

    def test_vtu_export_carries_bias_and_potential(self) -> None:
        runner = gate.find_runner()
        if runner is None:
            self.skipTest("vela_example_runner not built")
        fixture = gate.build_mesh()
        with _tmpdir() as tmp:
            workdir = Path(tmp)
            (workdir / "mesh.json").write_text(json.dumps(fixture["mesh"]), encoding="utf-8")
            samples = gate.run_points(runner, workdir, 0.0, [-0.5, -0.4], fixture["probes"], 1.0e17, 0.0, "probe")
            self.assertEqual({-0.5, -0.4}, set(samples))
            # Below flat band the p-type surface is accumulated (surface potential below bulk).
            self.assertLess(samples[-0.5], 0.0)
            self.assertGreater(samples[-0.4], 0.0)


if __name__ == "__main__":
    unittest.main()
