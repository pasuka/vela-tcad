"""Tests for ``scripts/reference_artifacts.py`` (WP9: register / fetch / verify).

A temporary directory mirrors the VM archive layout
(``<store>/<case>/v1/{source,mesh,plt,tdr,log}``) with small synthetic files and
is served through :class:`LocalDirBackend`; no Sentaurus VM access is needed.
The seven registered manifests under ``reference_tcad/`` are exercised in the
manifest-only ``--offline`` mode.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

import reference_artifacts as ra  # noqa: E402
from vela_sdevice import schema_validator as sv  # noqa: E402
from vela_sdevice.plt import build_plt_document  # noqa: E402

SCRIPT = REPO / "scripts" / "reference_artifacts.py"
C7 = REPO / "schemas" / "vela.reference_artifact_manifest.v1.schema.json"
CASE = "wp9_synthetic"
VERSION = "v1"

CMD_TEXT = """* synthetic sdevice deck for the WP9 tests
File {
  Grid      = "wp9_msh.tdr"
  Plot      = "wp9_des.tdr"
  Current   = "wp9_des.plt"
  Output    = "wp9_des.log"
  Parameter = "wp9.par"
}

Electrode {
  { Name="anode" Voltage=0.0 }   # swept electrode
  { Name="cathode" Voltage=0.0 }
}

Physics { Mobility(DopingDependence) Recombination(SRH) }

Math { Digits=5 Iterations=50 }

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled { Poisson Electron Hole }
  Quasistationary(
    InitialStep=0.1 MaxStep=0.25 MinStep=1e-5
    Goal { Name="anode" Voltage=1.0 }
  ) { Coupled { Poisson Electron Hole } }
}
"""

# Same deck as the repository would keep it: different comments and indentation, same tokens.
REPO_CMD_TEXT = """# repository copy (re-indented by Workbench, comments differ)
File {
    Grid = "wp9_msh.tdr"   Plot = "wp9_des.tdr"
    Current = "wp9_des.plt"   Output = "wp9_des.log"
    Parameter = "wp9.par"
}
Electrode {
    { Name="anode" Voltage=0.0 }
    { Name="cathode" Voltage=0.0 }
}
Physics { Mobility(DopingDependence) Recombination(SRH) }
Math { Digits=5 Iterations=50 }
Solve {
    Coupled(Iterations=100) { Poisson }
    Coupled { Poisson Electron Hole }
    Quasistationary( InitialStep=0.1 MaxStep=0.25 MinStep=1e-5
        Goal { Name="anode" Voltage=1.0 } ) { Coupled { Poisson Electron Hole } }
}
"""

PAR_TEXT = """Material = "Silicon" {
  * synthetic parameter file
  Epsilon { epsilon = 11.7 }
}
"""

MSH_CMD_TEXT = """Title ""
Controls { }
Definitions { }
Placements { }
"""

LOG_TEXT = """****************************************************************************
***                           Sentaurus Device                           ***
***                        Version T-2022.03-SP2                         ***
***                      (0.7745337, x86_64, Linux)                      ***
****************************************************************************

Tue Aug 18 10:27:03 2026: checked out 1 sdevice license(s)
Loading parameter file 'wp9.par' ...
Writing plot 'wp9_des.tdr' (TDR format) ... done.

Sentaurus Device simulation finished (Date: Tue Aug 18 10:27:03 2026  (CST)).
********************************* Good Bye ! *********************************
"""

FAKE_TDR = b"\x89HDF\r\n\x1a\n" + bytes(range(256)) * 4


def plt_text(final_bias: float, points: int = 5) -> str:
    rows = []
    for index in range(points):
        time = index / (points - 1)
        bias = final_bias * time
        rows.append({
            "time": time,
            "anode": {"outer_voltage_V": bias, "inner_voltage_V": bias, "e_current_A_per_um": -1e-12 * (index + 1),
                      "h_current_A_per_um": 1e-13 * (index + 1), "total_current_A_per_um": -9e-13 * (index + 1)},
            "cathode": {"outer_voltage_V": 0.0, "inner_voltage_V": 0.0, "e_current_A_per_um": 1e-12 * (index + 1),
                        "h_current_A_per_um": -1e-13 * (index + 1), "total_current_A_per_um": 9e-13 * (index + 1)},
        })
    return build_plt_document(["anode", "cathode"], rows)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lf_sha256(text: str) -> str:
    return hashlib.sha256(text.replace("\r\n", "\n").encode("utf-8")).hexdigest()


def temporary_directory(prefix: str) -> tempfile.TemporaryDirectory:
    build = REPO / "build"
    if build.is_dir():
        return tempfile.TemporaryDirectory(prefix=prefix, dir=str(build))
    return tempfile.TemporaryDirectory(prefix=prefix)


def build_mirror(store: Path, *, final_bias: float = 1.0) -> dict[str, Path]:
    base = store / CASE / VERSION
    files = {
        "source/wp9_sdevice.cmd": CMD_TEXT.encode("utf-8"),
        "source/wp9.par": PAR_TEXT.encode("utf-8"),
        "source/wp9_msh.cmd": MSH_CMD_TEXT.encode("utf-8"),
        "mesh/wp9_msh.tdr": FAKE_TDR,
        "plt/wp9_des.plt": plt_text(final_bias).encode("utf-8"),
        "tdr/wp9_des.tdr": FAKE_TDR[::-1],
        "log/wp9_des.log": LOG_TEXT.encode("utf-8"),
    }
    written = {}
    for rel, data in files.items():
        path = base.joinpath(*rel.split("/"))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        written[rel] = path
    return written


def build_repo_source(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "wp9_sdevice.cmd").write_text(REPO_CMD_TEXT, encoding="utf-8", newline="\n")
    (directory / "wp9.par").write_text(PAR_TEXT, encoding="utf-8", newline="\n")


def run_cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], capture_output=True, text=True, cwd=str(REPO))


class SyntheticStoreTest(unittest.TestCase):
    """register -> fetch -> verify against a local mirror of the VM archive."""

    def setUp(self) -> None:
        self._tmp = temporary_directory("wp9_ref_")
        self.tmp = Path(self._tmp.name)
        self.store = self.tmp / "store"
        self.repo_source = self.tmp / "repo_source"
        self.staging = self.tmp / "staging"
        self.reports = self.tmp / "reports"
        self.manifest_path = self.tmp / "artifact_manifest.json"
        self.missing_store = self.tmp / "missing_store"
        self.files = build_mirror(self.store)
        build_repo_source(self.repo_source)

    def tearDown(self) -> None:
        self._tmp.cleanup()

    # -- helpers -----------------------------------------------------------------

    def register(self, **extra: str) -> dict:
        args = ["register", "--case", CASE, "--version", VERSION, "--backend", f"local:{self.store}",
                "--manifest-out", str(self.manifest_path), "--repo-source-dir", str(self.repo_source),
                "--tdr-importer", "none", "--origin", "/vm/runs/wp9", "--note", "synthetic"]
        for key, value in extra.items():
            args.extend([key, value])
        proc = run_cli(*args)
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        return json.loads(self.manifest_path.read_text(encoding="utf-8"))

    def fetch(self) -> subprocess.CompletedProcess:
        return run_cli("fetch", "--case", CASE, "--manifest", str(self.manifest_path), "--dest", str(self.staging),
                       "--tdr-importer", "none")

    def verify(self, *extra: str, backend: str | None = None, dest: Path | None = None,
               report_dir: Path | None = None) -> tuple[subprocess.CompletedProcess, dict]:
        report_dir = report_dir or self.reports
        args = ["verify", "--case", CASE, "--manifest", str(self.manifest_path), "--dest", str(dest or self.staging),
                "--report-dir", str(report_dir), "--repo-source-dir", str(self.repo_source), "--tdr-importer", "none"]
        if backend:
            args.extend(["--backend", backend])
        args.extend(extra)
        proc = run_cli(*args)
        report_file = report_dir / "verify_report.json"
        report = json.loads(report_file.read_text(encoding="utf-8")) if report_file.is_file() else {}
        return proc, report

    def staged(self, rel: str) -> Path:
        return self.staging.joinpath(CASE, VERSION, *rel.split("/"))

    @staticmethod
    def check(report: dict, number: int) -> dict:
        return next(check for check in report["checks"] if check["number"] == number)

    # -- tests -------------------------------------------------------------------

    def test_register_fetch_verify_round_trip(self) -> None:
        manifest = self.register()
        self.assertEqual([], sv.collect_errors(sv.load_schema(C7), manifest))
        self.assertEqual(CASE, manifest["case"])
        self.assertEqual(VERSION, manifest["version"])
        self.assertEqual({"kind": "local", "root": str(self.store)}, manifest["store"])
        self.assertEqual("T-2022.03-SP2 (0.7745337)", manifest["sentaurus_version"])
        self.assertTrue(manifest["deck_matches_repository"])
        self.assertEqual({"vm_run_directory": "/vm/runs/wp9", "note": "synthetic"}, manifest["origin"])

        artifacts = {artifact["path"]: artifact for artifact in manifest["artifacts"]}
        self.assertEqual(set(self.files), set(artifacts))
        for rel, path in self.files.items():
            artifact = artifacts[rel]
            self.assertEqual(sha256_file(path), artifact["sha256"], rel)
            self.assertEqual(path.stat().st_size, artifact["bytes"], rel)
            self.assertEqual(ra.CATEGORY_BY_DIR[rel.split("/")[0]], artifact["category"], rel)
            self.assertEqual("registered", artifact["verify"]["status"], rel)
            self.assertEqual("/vm/runs/wp9", artifact["producer"]["run_directory"], rel)

        cmd_sha = lf_sha256(CMD_TEXT)
        par_sha = lf_sha256(PAR_TEXT)
        self.assertEqual("other", artifacts["source/wp9_sdevice.cmd"]["producer"]["tool"])   # comments differ
        self.assertEqual("repository", artifacts["source/wp9.par"]["producer"]["tool"])       # byte-identical
        self.assertEqual("other", artifacts["source/wp9_msh.cmd"]["producer"]["tool"])

        for rel in ("plt/wp9_des.plt", "tdr/wp9_des.tdr", "log/wp9_des.log"):
            producer = artifacts[rel]["producer"]
            self.assertEqual("sdevice", producer["tool"], rel)
            self.assertEqual("wp9_sdevice.cmd", producer["cmd"], rel)
            self.assertEqual(cmd_sha, producer["cmd_sha256"], rel)
            self.assertEqual("wp9.par", producer["parameter_file"], rel)
            self.assertEqual(par_sha, producer["parameter_sha256"], rel)
            self.assertRegex(producer["finished_utc"], r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")

        mesh = artifacts["mesh/wp9_msh.tdr"]["producer"]
        self.assertEqual("snmesh", mesh["tool"])
        self.assertEqual("wp9_msh.cmd", mesh["cmd"])
        self.assertEqual(lf_sha256(MSH_CMD_TEXT), mesh["cmd_sha256"])
        self.assertNotIn("summary", artifacts["mesh/wp9_msh.tdr"])
        self.assertNotIn("summary", artifacts["tdr/wp9_des.tdr"])

        plt_summary = artifacts["plt/wp9_des.plt"]["summary"]
        self.assertIn("anode OuterVoltage", plt_summary["datasets"])
        self.assertIn("cathode TotalCurrent", plt_summary["datasets"])
        self.assertEqual("time", plt_summary["datasets"][0])
        self.assertEqual(5, plt_summary["rows"])
        self.assertEqual(0.0, plt_summary["first_bias_V"])
        self.assertEqual(1.0, plt_summary["last_bias_V"])
        self.assertTrue(plt_summary["goal_reached"])
        self.assertEqual({"exit_status": "Good Bye"}, artifacts["log/wp9_des.log"]["summary"])

        proc = self.fetch()
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        for rel, path in self.files.items():
            self.assertTrue(self.staged(rel).is_file(), rel)
            self.assertEqual(sha256_file(path), sha256_file(self.staged(rel)), rel)
        proc = self.fetch()   # second run: everything already present
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        self.assertIn("'present': 7", proc.stdout)

        proc, report = self.verify()
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("pass", report["status"])
        self.assertEqual("online", report["mode"])
        self.assertTrue(report["compared_with_source"])
        self.assertTrue(report["store_reachable"])
        self.assertEqual(["pass"] * 6, [check["status"] for check in report["checks"]])
        self.assertEqual([], report["failures"])
        self.assertEqual([], report["partial_references"])
        structure = self.check(report, 3)
        self.assertTrue(any("reaches Goal 1 V" in item["detail"] for item in structure["items"]), structure)
        completeness = self.check(report, 4)
        self.assertTrue(any("not inspected (HDF5)" in item["detail"] for item in completeness["items"]), completeness)
        self.assertTrue(any(item["subject"] == "source/wp9_sdevice.cmd" and "equals repository" in item["detail"]
                            for item in self.check(report, 2)["items"]), self.check(report, 2))
        markdown = list(self.reports.glob(f"verify_{CASE}_*.md"))
        self.assertEqual(1, len(markdown), markdown)
        self.assertIn("**PASS**", markdown[0].read_text(encoding="utf-8"))

        proc, report = self.verify("--offline", report_dir=self.tmp / "offline_reports")
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("offline", report["mode"])
        self.assertEqual(["skip", "pass", "skip", "skip", "pass", "pass"],
                         [check["status"] for check in report["checks"]])

    def test_flipped_byte_in_staged_copy_fails_check_1_naming_the_file(self) -> None:
        self.register()
        self.assertEqual(0, self.fetch().returncode)
        staged_plt = self.staged("plt/wp9_des.plt")
        data = bytearray(staged_plt.read_bytes())
        data[len(data) // 2] ^= 0x01
        staged_plt.write_bytes(bytes(data))

        proc, report = self.verify()
        self.assertEqual(1, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("fail", report["status"])
        self.assertEqual("fail", self.check(report, 1)["status"])
        failed = [item for item in self.check(report, 1)["items"] if item["status"] == "fail"]
        self.assertEqual(["plt/wp9_des.plt"], [item["subject"] for item in failed], failed)
        self.assertIn("staged copy", failed[0]["detail"])
        self.assertIn("plt/wp9_des.plt", proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)

    def test_edited_staged_cmd_fails_source_consistency(self) -> None:
        self.register()
        self.assertEqual(0, self.fetch().returncode)
        staged_cmd = self.staged("source/wp9_sdevice.cmd")
        staged_cmd.write_text(CMD_TEXT.replace("Voltage=1.0", "Voltage=2.0"), encoding="utf-8", newline="\n")

        # Store unreachable -> the local copies are what gets verified.
        proc, report = self.verify(backend=f"local:{self.missing_store}")
        self.assertEqual(1, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("local_copy", report["mode"])
        self.assertFalse(report["compared_with_source"])
        self.assertEqual("fail", self.check(report, 2)["status"])
        failures = [item for item in self.check(report, 2)["items"] if item["status"] == "fail"]
        self.assertTrue(failures, self.check(report, 2))
        self.assertTrue(all(item["subject"] == "source/wp9_sdevice.cmd" for item in failures), failures)
        self.assertTrue(any("cmd_sha256" in item["detail"] for item in failures), failures)
        self.assertTrue(any("differs from repository" in item["detail"] for item in failures), failures)
        self.assertTrue(any("wp9_sdevice.cmd" in line for line in report["failures"]), report["failures"])

    def test_tampered_cmd_sha256_fails_check_2_only(self) -> None:
        manifest = self.register()
        self.assertEqual(0, self.fetch().returncode)
        for artifact in manifest["artifacts"]:
            if "cmd_sha256" in artifact["producer"] and artifact["producer"]["cmd"] == "wp9_sdevice.cmd":
                artifact["producer"]["cmd_sha256"] = "0" * 64
        self.manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        proc, report = self.verify()
        self.assertEqual(1, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual(["pass", "fail", "pass", "pass", "pass", "pass"],
                         [check["status"] for check in report["checks"]])
        self.assertTrue(any("producer.cmd_sha256" in item["detail"] and item["status"] == "fail"
                            for item in self.check(report, 2)["items"]), self.check(report, 2))

    def test_unreachable_store_without_local_copy_is_unable_to_verify(self) -> None:
        self.register()
        proc, report = self.verify(backend=f"local:{self.missing_store}", dest=self.tmp / "empty_staging")
        self.assertEqual(2, proc.returncode, proc.stdout + proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertNotIn("Traceback", proc.stdout)
        self.assertEqual("unable_to_verify", report["status"])
        self.assertEqual(2, report["exit_code"])
        self.assertFalse(report["store_reachable"])
        self.assertEqual(["skip"] * 5 + ["pass"], [check["status"] for check in report["checks"]])
        console = (proc.stdout + proc.stderr).lower()
        self.assertTrue("unable to verify" in console or "unable_to_verify" in console, console)

        # The Python API reports the same without touching the store.
        backend = ra.LocalDirBackend(self.missing_store)
        self.assertFalse(backend.is_reachable())
        api_report = ra.verify_case(json.loads(self.manifest_path.read_text(encoding="utf-8")),
                                    manifest_path=self.manifest_path, backend=backend,
                                    staging_root=self.tmp / "empty_staging", report_dir=self.tmp / "api_reports",
                                    repo_source_dir=self.repo_source)
        self.assertEqual("unable_to_verify", api_report["status"])
        self.assertEqual(ra.EXIT_UNAVAILABLE, api_report["exit_code"])

    def test_fetch_reports_store_copy_that_no_longer_matches(self) -> None:
        self.register()
        store_plt = self.files["plt/wp9_des.plt"]
        store_plt.write_bytes(store_plt.read_bytes() + b"\n")
        proc = self.fetch()
        self.assertEqual(1, proc.returncode, proc.stdout + proc.stderr)
        self.assertIn("mismatch: plt/wp9_des.plt", proc.stderr)
        self.assertTrue(self.staged("log/wp9_des.log").is_file())

    def test_register_unreachable_store_exits_2(self) -> None:
        proc = run_cli("register", "--case", CASE, "--backend", f"local:{self.missing_store}",
                       "--manifest-out", str(self.manifest_path), "--tdr-importer", "none")
        self.assertEqual(2, proc.returncode, proc.stdout + proc.stderr)
        self.assertNotIn("Traceback", proc.stderr)
        self.assertFalse(self.manifest_path.exists())

    def test_partial_reference_when_goal_not_reached(self) -> None:
        self.files["plt/wp9_des.plt"].write_text(plt_text(0.5), encoding="utf-8", newline="\n")
        manifest = self.register()
        summary = next(a for a in manifest["artifacts"] if a["path"] == "plt/wp9_des.plt")["summary"]
        self.assertEqual(0.5, summary["last_bias_V"])
        self.assertFalse(summary["goal_reached"])
        self.assertEqual(0, self.fetch().returncode)
        proc, report = self.verify()
        self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("pass", report["status"])
        self.assertEqual(1, len(report["partial_references"]), report["partial_references"])
        self.assertIn("plt/wp9_des.plt", report["partial_references"][0])
        self.assertIn("Goal 1 V not reached", report["partial_references"][0])

    def test_missing_electrode_dataset_fails_structure_check(self) -> None:
        self.register()
        self.assertEqual(0, self.fetch().returncode)
        staged_plt = self.staged("plt/wp9_des.plt")
        text = staged_plt.read_text(encoding="utf-8").replace('"cathode OuterVoltage"', '"gate OuterVoltage"')
        staged_plt.write_text(text, encoding="utf-8", newline="\n")
        self.files["plt/wp9_des.plt"].write_text(text, encoding="utf-8", newline="\n")
        proc, report = self.verify(backend=f"local:{self.missing_store}")
        self.assertEqual(1, proc.returncode, proc.stdout + proc.stderr)
        self.assertEqual("fail", self.check(report, 3)["status"])
        self.assertTrue(any("cathode" in item["detail"] for item in self.check(report, 3)["items"]
                            if item["status"] == "fail"), self.check(report, 3))


class RegisteredManifestsOfflineTest(unittest.TestCase):
    """The seven VM-registered manifests pass the manifest-only checks (2 and 5)."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.manifests = sorted((REPO / "reference_tcad").glob("*/artifact_manifest.json"))

    def test_seven_manifests_are_registered(self) -> None:
        self.assertGreaterEqual(len(self.manifests), 7, self.manifests)

    def test_offline_verify_passes_for_every_registered_manifest(self) -> None:
        with temporary_directory("wp9_offline_") as tmp:
            for path in self.manifests:
                case = path.parent.name
                with self.subTest(case=case):
                    manifest = ra.load_manifest(path)
                    report = ra.verify_case(manifest, manifest_path=path, backend=None,
                                            staging_root=Path(tmp) / "staging", report_dir=Path(tmp) / case,
                                            offline=True)
                    self.assertEqual("pass", report["status"], report["failures"])
                    self.assertEqual(ra.EXIT_OK, report["exit_code"])
                    self.assertEqual("offline", report["mode"])
                    statuses = {check["number"]: check["status"] for check in report["checks"]}
                    self.assertEqual({1: "skip", 2: "pass", 3: "skip", 4: "skip", 5: "pass", 6: "pass"}, statuses)
                    self.assertTrue((Path(tmp) / case / "verify_report.json").is_file())
                    self.assertTrue(list((Path(tmp) / case).glob(f"verify_{case}_*.md")))
                    check_2 = next(check for check in report["checks"] if check["number"] == 2)
                    self.assertTrue(check_2["items"], case)
                    self.assertFalse([item for item in check_2["items"] if item["status"] == "fail"], check_2)

    def test_offline_cli_exit_code(self) -> None:
        with temporary_directory("wp9_offline_cli_") as tmp:
            proc = run_cli("verify", "--case", "pn2d_sentaurus2018", "--offline", "--report-dir", tmp,
                           "--dest", str(Path(tmp) / "staging"))
            self.assertEqual(0, proc.returncode, proc.stdout + proc.stderr)
            self.assertIn("[offline]: pass", proc.stdout)
            report = json.loads((Path(tmp) / "verify_report.json").read_text(encoding="utf-8"))
            self.assertEqual("pn2d_sentaurus2018", report["case"])
            check_2 = next(check for check in report["checks"] if check["number"] == 2)
            subjects = {item["subject"] for item in check_2["items"]}
            self.assertIn("source/pn2d_0v_sdevice.cmd", subjects)


class HelperTest(unittest.TestCase):
    def test_sdevice_output_name_rules(self) -> None:
        self.assertEqual("pn2d_0v.plt", ra.sdevice_output_name("pn2d_0v.plt", "plt"))
        self.assertEqual("schottky_n_iv_des.plt", ra.sdevice_output_name("schottky_n_iv", "plt"))
        self.assertEqual("n6_des.tdr", ra.sdevice_output_name("n6", "tdr"))
        self.assertEqual("pn2d_0v.log_des.log", ra.sdevice_log_name("pn2d_0v.log"))
        self.assertEqual("schottky_n_des.log", ra.sdevice_log_name("schottky_n_des.log"))

    def test_deck_analysis(self) -> None:
        deck = ra.Deck("schottky_n_des.cmd", (REPO / "reference_tcad" / "schottky_charon_sentaurus2018" / "source"
                                             / "schottky_n_des.cmd").read_bytes())
        self.assertEqual("sdevice", deck.kind)
        self.assertEqual(["anode", "cathode"], deck.electrodes)
        self.assertEqual({"Grid": "schottky_n_msh.tdr", "Plot": "schottky_n_des.tdr", "Current": "schottky_n_iv",
                          "Output": "schottky_n_des.log"}, deck.file_names)
        self.assertIsNone(deck.parameter_file)
        self.assertEqual(["forward_"], deck.current_prefixes)
        self.assertEqual([ra.Goal("forward_", "anode", 1.0)], deck.goals)
        outputs = deck.expected_outputs()
        self.assertEqual("plt", outputs["schottky_n_iv_des.plt"])
        self.assertEqual("plt", outputs["forward_schottky_n_iv_des.plt"])
        self.assertEqual("tdr", outputs["schottky_n_des.tdr"])
        self.assertEqual("log", outputs["schottky_n_des.log"])
        self.assertEqual(ra.Goal("forward_", "anode", 1.0), deck.goal_for_output("forward_schottky_n_iv_des.plt"))
        self.assertIsNone(deck.goal_for_output("schottky_n_iv_des.plt"))
        self.assertEqual(4, deck.match_score("forward_schottky_n_iv_des.plt"))
        self.assertEqual(0, deck.match_score("unrelated.plt"))

    def test_deck_kinds(self) -> None:
        self.assertEqual("snmesh", ra.Deck("x_msh.cmd", MSH_CMD_TEXT.encode()).kind)
        self.assertEqual("sde", ra.Deck("x_sde.cmd", b'(sdegeo:create-rectangle (position 0 0 0) (position 1 1 0) "Silicon" "R")\n(sde:build-mesh "" "x")\n').kind)
        self.assertEqual("sprocess", ra.Deck("sprocess.cmd", b"line x location=0.0 spacing=1<nm> tag=top\ninit silicon\nstruct tdr=n1\n").kind)
        self.assertEqual("other", ra.Deck("notes.cmd", b"just text\n").kind)

    def test_deck_equivalence_ignores_comments_and_whitespace(self) -> None:
        self.assertTrue(ra.decks_equivalent(CMD_TEXT.encode(), REPO_CMD_TEXT.encode()))
        self.assertTrue(ra.decks_equivalent(CMD_TEXT.encode(), CMD_TEXT.replace("\n", "\r\n").encode()))
        self.assertFalse(ra.decks_equivalent(CMD_TEXT.encode(), CMD_TEXT.replace("Voltage=1.0", "Voltage=1.5").encode()))
        self.assertEqual(lf_sha256(CMD_TEXT), ra.lf_sha256(CMD_TEXT.replace("\n", "\r\n").encode("utf-8")))
        self.assertEqual(lf_sha256(CMD_TEXT), ra.lf_sha256(b"\xef\xbb\xbf" + CMD_TEXT.encode("utf-8")))

    def test_version_banner_and_exit_status(self) -> None:
        self.assertEqual("T-2022.03-SP2 (0.7745337)", ra.parse_version_banner(LOG_TEXT))
        self.assertEqual("T-2022.03-SP2", ra.parse_version_banner("*** Version T-2022.03-SP2 ***"))
        self.assertIsNone(ra.parse_version_banner("no banner"))
        self.assertEqual("Good Bye", ra.log_exit_status(LOG_TEXT))
        self.assertEqual("no Good Bye marker", ra.log_exit_status("Sentaurus Device ... error: stopped"))

    def test_plt_summary(self) -> None:
        summary = ra.summarise_plt(plt_text(1.0), ra.Goal("", "anode", 1.0))
        self.assertEqual(5, summary["rows"])
        self.assertEqual((0.0, 1.0, True), (summary["first_bias_V"], summary["last_bias_V"], summary["goal_reached"]))
        summary = ra.summarise_plt(plt_text(0.75))
        self.assertEqual(0.75, summary["last_bias_V"])
        self.assertNotIn("goal_reached", summary)

    def test_tdr_inventory_summary(self) -> None:
        inventory = {
            "vertex_count": 6,
            "regions": [
                {"name": "Si", "material": "Silicon", "type": 0, "triangles": 4, "edges": 0, "points": 0},
                {"name": "anode", "material": "Contact", "type": 1, "triangles": 0, "edges": 2, "points": 0},
            ],
            "geometry": {"vertices": [], "regions": [
                {"name": "Si", "edges": [], "points": []},
                {"name": "anode", "edges": [[0, 1], [1, 2]], "points": []},
            ]},
            "fields": [{"name": "ElectrostaticPotential", "region": 0}, {"name": "eDensity", "region": 0},
                       {"name": "ElectrostaticPotential", "region": 1}],
        }
        self.assertEqual({"nodes": 6, "elements": 4, "contacts": [{"name": "anode", "node_count": 3}],
                          "materials": ["Silicon"], "datasets": ["ElectrostaticPotential", "eDensity"]},
                         ra.summarise_tdr_inventory(inventory))

    def test_backend_factory_and_ssh_defaults(self) -> None:
        local = ra.make_backend(f"local:{REPO}")
        self.assertIsInstance(local, ra.LocalDirBackend)
        self.assertTrue(local.is_reachable())
        ssh = ra.make_backend("ssh:sentaurus:/home/tcad/vela_reference_artifacts")
        self.assertIsInstance(ssh, ra.SshBackend)
        self.assertEqual({"kind": "ssh", "host": "sentaurus", "root": "/home/tcad/vela_reference_artifacts"},
                         ssh.store_descriptor())
        self.assertTrue(ssh.ssh.endswith("ssh") or ssh.ssh.lower().endswith("ssh.exe"))
        if os.name == "nt" and (Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "OpenSSH" / "ssh.exe").is_file():
            self.assertIn("OpenSSH", ssh.ssh)
        from_store = ra.make_backend(None, {"kind": "ssh", "host": "sentaurus", "root": "/r"})
        self.assertIsInstance(from_store, ra.SshBackend)
        self.assertIsNone(ra.make_backend(None, {"kind": "github_release", "root": "owner/repo"}))
        with self.assertRaises(ValueError):
            ra.make_backend("ftp:host")
        with self.assertRaises(ValueError):
            ra.make_backend("local:")

    def test_local_backend_operations(self) -> None:
        with temporary_directory("wp9_backend_") as tmp:
            store = Path(tmp) / "store"
            files = build_mirror(store)
            backend = ra.LocalDirBackend(store)
            self.assertTrue(backend.is_reachable())
            listed = backend.list_files(f"{CASE}/{VERSION}/source")
            self.assertEqual([f"{CASE}/{VERSION}/source/wp9.par", f"{CASE}/{VERSION}/source/wp9_msh.cmd",
                              f"{CASE}/{VERSION}/source/wp9_sdevice.cmd"], listed)
            self.assertEqual([], backend.list_files(f"{CASE}/{VERSION}/nothing"))
            rel = f"{CASE}/{VERSION}/plt/wp9_des.plt"
            self.assertEqual(sha256_file(files["plt/wp9_des.plt"]), backend.sha256(rel))
            self.assertEqual(files["plt/wp9_des.plt"].stat().st_size, backend.size(rel))
            self.assertEqual((backend.size(rel), backend.sha256(rel)), backend.stat_and_hash(rel))
            target = Path(tmp) / "copy" / "wp9_des.plt"
            backend.fetch(rel, target)
            self.assertEqual(files["plt/wp9_des.plt"].read_bytes(), target.read_bytes())
            self.assertEqual(files["plt/wp9_des.plt"].read_bytes(), backend.read_bytes(rel))
            with self.assertRaises(ra.TransportError):
                backend.sha256(f"{CASE}/{VERSION}/plt/absent.plt")
            self.assertFalse(ra.LocalDirBackend(Path(tmp) / "absent").is_reachable())


if __name__ == "__main__":
    unittest.main()
