"""M4 public Gmsh device orchestration and internal TDX gate tests."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import sentaurus_import as device  # noqa: E402
import sentaurus_gmsh_tdx_gate as tdx_gate  # noqa: E402


SDE = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_sde.cmd"
SDEVICE_IV = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" /
    "pn2d_iv_sdevice.cmd"
)


class M4ContractUnitTest(unittest.TestCase):
    def test_frontend_backend_matrix_rejects_public_official_tdr(self) -> None:
        for frontend in ("offline-static", "sentaurus2022-evaluated"):
            for backend in ("gmsh", "tensor-legacy"):
                device._validate_device_backend(frontend, backend)
        device._validate_device_backend("sentaurus2022-evaluated", "official-tdr")
        with self.assertRaisesRegex(device.ExecutionIrError, "requires the internal"):
            device._validate_device_backend("offline-static", "official-tdr")

    def test_gmsh_mesh_contract_override_accepts_barycentric_obtuse_mesh(self) -> None:
        case = device.load_e2e_contract()["cases"]["pn2d_sentaurus2018_iv"]
        qualification = {
            "node_count": 1000, "cell_count": 1900,
            "max_angle_degrees": 135.0, "node_volume_policy": "barycentric",
            "contacts": [{"node_count": 11}, {"node_count": 11}],
        }
        with tempfile.TemporaryDirectory(prefix="vela_m4_contract_") as tmp:
            root = Path(tmp)
            curve = root / "curve.csv"
            terminal = root / "terminal.csv"
            with curve.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "mode", "bias_V", "converged", "current_total_A_per_um"])
                writer.writeheader()
                for index in range(10):
                    writer.writerow({
                        "mode": "iv", "bias_V": 0.3 * index / 9,
                        "converged": 1, "current_total_A_per_um": index * 1.0e-12,
                    })
            with terminal.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=[
                    "bias_V", "contact", "I_sgflux_A_per_um"])
                writer.writeheader()
                for index in range(10):
                    bias = 0.3 * index / 9
                    current = index * 1.0e-12
                    writer.writerow({"bias_V": bias, "contact": "Anode",
                                     "I_sgflux_A_per_um": current})
                    writer.writerow({"bias_V": bias, "contact": "Cathode",
                                     "I_sgflux_A_per_um": -current})
            report = device.evaluate_e2e_acceptance(
                "pn2d_sentaurus2018_iv", case, "iv", qualification,
                curve, terminal, mesh_backend="gmsh")
        self.assertEqual("accepted", report["status"])
        self.assertEqual("gmsh", report["mesh_backend"])

    def test_parser_defaults_to_public_gmsh_report_only(self) -> None:
        parser = device.build_parser()
        args = parser.parse_args([
            "device", "--sde", str(SDE), "--sdevice", str(SDEVICE_IV),
            "--output-dir", "out",
        ])
        self.assertEqual("offline-static", args.sde_frontend)
        self.assertEqual("gmsh", args.mesh_backend)
        self.assertEqual("report-only", args.mesh_quality_profile)
        self.assertIsNone(args.node_volume_policy)


class M4TdxGateUnitTest(unittest.TestCase):
    def test_dry_run_hash_seals_commands_without_invoking_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m4_tdx_") as tmp:
            root = Path(tmp)
            preprocess = root / "preprocess"
            output = root / "output"
            preprocess.mkdir()
            (preprocess / "device.grd").write_text("grd", encoding="utf-8")
            (preprocess / "device.dat").write_text("dat", encoding="utf-8")
            (preprocess / "mesh.json").write_text("{}", encoding="utf-8")
            (preprocess / "doping.csv").write_text("node_id\n", encoding="utf-8")
            manifest = {
                "mesh_backend": "gmsh", "gmsh_version": "4.15.2",
                "cross_format": {"status": "passed"},
                "sha256": {
                    name: tdx_gate.sha256_file(preprocess / name)
                    for name in ("device.grd", "device.dat")
                },
            }
            (preprocess / "gmsh_preprocess_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8")
            report = tdx_gate.run_gate(argparse.Namespace(
                preprocess_dir=preprocess, output_dir=output,
                run_id="m4_dry", ssh_target="sentaurus",
                remote_root="~/sentaurus_runs/vela_gmsh_tdx",
                dry_run=True, ssh_bin=None, scp_bin=None,
                tdr_importer=REPO / "build-release" / "sentaurus_import.exe",
            ))
        self.assertEqual("prepared", report["status"])
        self.assertEqual("tdx -h", report["commands"][0])
        self.assertIn("tdx -d device.grd device.dat device.tdr", report["commands"][2])


class M4GmshIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        prefix = Path(sys.prefix)
        cls.wrapper = prefix / "lib" / "gmsh.py"
        cls.library = prefix / "bin" / "libgmsh.dll"
        cls.executable = prefix / "bin" / "gmsh.exe"
        if not all(path.is_file() for path in (cls.wrapper, cls.library, cls.executable)):
            raise unittest.SkipTest("Gmsh 4.15.2 UCRT64 runtime is not installed")

    def test_default_device_command_generates_complete_gmsh_pn2d_bundle(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m4_device_") as tmp:
            root = Path(tmp)
            output = root / "output"
            completed = subprocess.run([
                sys.executable, str(REPO / "scripts" / "sentaurus_import.py"),
                "device", "--sde", str(SDE), "--sdevice", str(SDEVICE_IV),
                "--output-dir", str(output),
                "--gmsh-runtime-dir", str(root / "runtime"),
                "--gmsh-python", str(self.wrapper),
                "--gmsh-library", str(self.library),
                "--gmsh-executable", str(self.executable),
            ], cwd=REPO, text=True, capture_output=True, check=False)
            self.assertEqual(0, completed.returncode, completed.stderr)
            manifest = json.loads((output / "run_manifest.json").read_text())
            deck = json.loads((output / "simulation.json").read_text())
            expected = {
                "device.msh", "device.vtu", "device.vtk", "device.grd", "device.dat",
                "mesh.json", "doping.csv", "id_mapping.json",
                "gmsh_preprocess_manifest.json", "cross_format_report.json",
            }
            self.assertTrue(expected.issubset({path.name for path in output.iterdir()}))
            self.assertEqual("generated", manifest["status"])
            self.assertEqual("gmsh", manifest["mesh_backend"])
            self.assertEqual("passed", manifest["mesh_source"]["cross_format_status"])
            self.assertEqual("barycentric", deck["mesh_geometry"]["node_volume_policy"])
            self.assertFalse(deck["mesh_geometry"]["require_non_obtuse"])


if __name__ == "__main__":
    unittest.main()
