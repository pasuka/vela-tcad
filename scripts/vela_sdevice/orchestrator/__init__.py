"""SDevice frontend gates and stage execution (WP5)."""
from __future__ import annotations
import json
import os
from pathlib import Path
import queue
import re
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time

from ..cmd_parser import parse_cmd_file
from ..deck_emit import DeckEmitError, EmitOptions, bind_current_start, emit_stage_plan, rebase_stage, write_stage_plan
from ..materials import MaterialsError, active_models_from_ir, emit_materials, write_materials
from ..mesh_import import MeshImportError, import_grid, check_electrodes
from ..resources import resource_root
from .. import schema_validator as sv
from .outputs import OutputPublisher, accepted_records
from .support import RunError, RunLock, write_json


def executable(name, bin_dir):
    variable = "VELA_EXAMPLE_RUNNER" if name == "vela_example_runner" else "VELA_SENTAURUS_IMPORT"
    override = os.environ.get(variable)
    if override:
        path = Path(override).resolve()
        if not path.is_file():
            raise RunError("E202", f"Option '{variable}' names a missing executable")
        return [sys.executable, str(path)] if path.suffix == ".py" else [str(path)]
    suffix = ".exe" if os.name == "nt" else ""
    for root in (bin_dir, bin_dir.parent):
        path = root / (name + suffix)
        if path.is_file():
            return [str(path.resolve())]
    found = shutil.which(name + suffix)
    return [found] if found else None


def version_info():
    root = resource_root()
    if (root / "version.json").is_file():
        return json.loads((root / "version.json").read_text(encoding="utf-8"))
    version = re.search(r"project\([^\s)]+\s+VERSION\s+(\S+)",
                        (root / "CMakeLists.txt").read_text()).group(1)
    done = subprocess.run(["git", "-c", "core.fsmonitor=false", "rev-parse", "HEAD"],
                          cwd=root, capture_output=True, text=True)
    return {"version": version, "commit": done.stdout.strip() if done.returncode == 0 else "unknown"}


class DeviceRun:
    def __init__(self, cmd, options):
        self.cmd, self.options = Path(cmd).resolve(), options
        self.ir = parse_cmd_file(self.cmd)
        stem = self.cmd.stem
        if not re.fullmatch(r"[A-Za-z0-9._-]+", stem):
            stem = "device_" + self.ir["source"]["sha256"][:8]
        self.work = self.cmd.parent / (stem + "_vela")
        self.log_path = self.cmd.parent / self.ir["files"]["output"]["resolved"]
        self.lock = RunLock(self.work / ".lock")
        self.stream = None
        self.stage = None
        self.proc = None
        self.interrupted = False
        self.signal_number = None
        self.handlers = {}
        self.result = None
        self.publisher = None
        self.manifest = {"schema": "vela.sdevice_run_manifest.v1", "cmd": str(self.cmd),
            "sha256": self.ir["source"]["sha256"], "version": version_info(), "check_only": options.check,
            "stages": [], "outputs": [], "decisions": [], "sidecar_overrides": [], "exit_code": None}

    def log(self, level, text):
        if self.stream:
            self.stream.write(f"[{level}] {text}\n")
            self.stream.flush()

    def _signals(self):
        def handler(number, frame):
            self.interrupted, self.signal_number = True, number
        for name in ("SIGINT", "SIGTERM", "SIGBREAK"):
            number = getattr(signal, name, None)
            if number is not None:
                self.handlers[number] = signal.signal(number, handler)

    def _restore_signals(self):
        for number, handler in self.handlers.items():
            signal.signal(number, handler)

    def _stop_child(self):
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.proc.kill()
                self.proc.wait(timeout=5)

    def run(self):
        exit_code = 0
        try:
            self.work.mkdir(parents=True, exist_ok=True)
            self.lock.acquire()
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            self.stream = self.log_path.open("a", encoding="utf-8")
            self._signals()
            self.log("HEADER", json.dumps({**self.manifest["version"], "cmd": str(self.cmd),
                "sha256": self.manifest["sha256"], "time": time.time(), "host": socket.gethostname(), "pid": os.getpid()}))
            self.log("CMD", "Command file echo begins")
            self.stream.write(self.cmd.read_text(encoding="utf-8-sig", errors="replace") + "\n")
            self.log("CMD", "Command file echo ends")
            if self.lock.stale:
                owner = self.lock.previous
                self.log("WARNING", f"W701: Stale lock from PID {owner['pid']} on {owner['host']} replaced")
            for warning in self.ir["warnings"]:
                self.log("WARNING", json.dumps(warning))
            if not self.ir["verdict"]["accepted"]:
                for error in self.ir["errors"]:
                    self.log("FATAL", json.dumps(error))
                identifier = min(self.ir["verdict"]["fatal_ids"], key=lambda x: int(x[1]))
                raise RunError(identifier, "Command file failed input gates", findings=self.ir["errors"])
            write_json(self.work / "execution_ir.json", self.ir)
            (self.work / "fail_report.json").unlink(missing_ok=True)
            for key, identifier in (("grid", "E302"), ("parameter", "E303")):
                path = self.cmd.parent / self.ir["files"][key]["resolved"]
                if not path.is_file():
                    raise RunError(identifier, f"{key.capitalize()} file '{path}' not found")
            location = Path(os.environ.get("VELA_BIN", Path(sys.argv[0]).resolve().parent))
            runner = executable("vela_example_runner", location)
            importer = executable("sentaurus_import", location)
            mesh = import_grid(self.cmd.parent / self.ir["files"]["grid"]["resolved"], self.work,
                               importer=importer[0] if importer else None)
            findings = [f for f in check_electrodes(mesh, [e["name"] for e in self.ir["electrodes"]]) if not f["matched"]]
            if findings:
                raise RunError("E401", "Electrodes have no matching mesh contacts", findings=findings)
            sv.validate(sv.load_schema(resource_root() / "schemas/vela.mesh_qualification.v1.schema.json"), mesh)
            write_json(self.work / "mesh_qualification.json", mesh)
            electrode_materials = [e["attributes"]["material"] for e in self.ir["electrodes"] if "material" in e["attributes"]]
            materials = emit_materials(self.cmd.parent / self.ir["files"]["parameter"]["resolved"],
                active_models=active_models_from_ir(self.ir), temperature_K=self.ir["physics"]["temperature_K"] or 300,
                mesh_materials=mesh["materials"], electrode_materials=electrode_materials,
                effective_models=self.ir["physics"]["effective"], allow_lossy=self.options.allow_lossy,
                band_edge_sensitive=any(e["attributes"].get("schottky") for e in self.ir["electrodes"]))
            write_materials(self.work, materials)
            for kind in ("lossy", "inert"):
                for parameter in materials.report[kind]:
                    self.log("MATERIAL", json.dumps({"status": kind, **parameter}))
            sidecar_path = self.cmd.with_suffix(".vela.json")
            try:
                sidecar = json.loads(sidecar_path.read_text(encoding="utf-8")) if sidecar_path.exists() else None
            except (OSError, ValueError) as error:
                raise RunError("E409", f"Sidecar '{sidecar_path}' is invalid: {error}") from error
            self.result = emit_stage_plan(self.ir, mesh_qualification=mesh, materials_report=materials.report,
                work_dir=self.work.name, cmd_basename=self.cmd.name, cmd_sha256=self.ir["source"]["sha256"],
                options=EmitOptions(exit_on_failure=self.options.exit_on_failure, threads=self.options.threads,
                    plot_final_only=self.options.plot_final_only, vtu_binary=self.options.vtu_binary,
                    keep_intermediate=self.options.keep_intermediate, allow_lossy=self.options.allow_lossy,
                    defer_current_times=True), sidecar=sidecar)
            self._write_plan()
            self.log("INFO", "Multi-stage PLT time resets per stage; sign factor +1 (WP10 pending)")
            if self.interrupted:
                raise RunError("E801", f"Interrupted by signal {self.signal_number} during input gates")
            if self.options.check:
                self.log("INFO", "Check-only gates passed; runner not invoked")
            else:
                if runner is None:
                    raise RunError("E202", "Option 'runner' cannot be resolved in VELA_BIN or PATH")
                self.publisher = OutputPublisher(self.cmd.parent, self.work, self.result)
                count = 1 if self.options.initial_only else len(self.result.decks)
                for index in range(count):
                    self.stage = index
                    if self.interrupted:
                        raise RunError("E801", f"Interrupted by signal {self.signal_number} before stage {index}")
                    stage = self.result.stage_plan["stages"][index]
                    if stage.get("sweep", {}).get("mode") == "current":
                        contact = stage["sweep"]["contact"]
                        if not index or not records:
                            raise RunError("E407", "Current stage has no preceding accepted state")
                        current = records[-1].get(f"contact_{contact}_total_A_per_um")
                        if current is None and records[-1].get("current_contact") == contact:
                            current = records[-1].get("current_total_A_per_um")
                        if current is None:
                            raise RunError("E407", f"Previous CSV has no measured current for '{contact}'")
                        self.result = bind_current_start(self.result, index, float(current), voltages)
                        self.publisher.result = self.result
                        self._write_plan()
                        stage = self.result.stage_plan["stages"][index]
                        self.log("INFO", f"Current Time mapping bound to {current} A/um at stage {index}")
                    self.log("STAGE", f"Begin stage {index}: {stage['kind']}")
                    status, records, stderr = self._execute(index, runner)
                    reached, voltages = self._reached(index, records)
                    runtime = {"index": index, "runner_exit_code": status, "goal_reached": reached,
                               "accepted_points": len(records), "reached_voltages_V": voltages}
                    self.manifest["stages"].append(runtime)
                    self.log("STAGE", json.dumps(runtime))
                    if self.interrupted:
                        raise RunError("E801", f"Interrupted by signal {self.signal_number} at stage {index}; partial outputs kept")
                    if status == 2:
                        raise RunError("E406", f"Generated deck '{stage['deck']}' was rejected by runner", stderr=stderr)
                    if status not in (0, 1):
                        raise RunError("E503", f"Solver terminated abnormally (exit code {status})", stderr=stderr)
                    if index == 0 and (status != 0 or not records):
                        raise RunError("E501", "Equilibrium stage failed", stderr=stderr)
                    if (status != 0 or not reached) and self.result.stage_plan["exit_on_failure"]:
                        raise RunError("E502", f"Stage {index} failed with ExitOnFailure", stderr=stderr)
                    if not records or not (self.work / stage["output_state"]).is_file():
                        raise RunError("E503", f"Stage {index} produced no accepted state for the next stage", stderr=stderr)
                    if index + 1 < count and (not reached or stage.get("sweep", {}).get("mode") == "current"):
                        self.result = rebase_stage(self.result, index + 1, voltages, update_saved_state=True)
                        self.publisher.result = self.result
                        self.log("WARNING", f"W206: stage {index} reached {voltages}; rebasing later decks")
                        self._write_plan()
        except (DeckEmitError, MaterialsError, MeshImportError) as error:
            exit_code = self._failure(RunError(error.message_id, getattr(error, "message", str(error))))
        except RunError as error:
            exit_code = self._failure(error)
        except (sv.SchemaError, sv.ValidationError, ValueError, KeyError, TypeError) as error:
            exit_code = self._failure(RunError("E406", f"Input/deck validation failed: {error}"))
        except OSError as error:
            exit_code = self._failure(RunError("E601", f"Cannot write output file: {error}"))
        finally:
            self._stop_child()
            if self.interrupted and exit_code == 0:
                exit_code = self._failure(RunError("E801", f"Interrupted by signal {self.signal_number}; partial outputs kept"))
            if self.lock.owned:
                self.manifest["exit_code"] = exit_code
                self.manifest["outputs"] = sorted(self.publisher.outputs) if self.publisher else []
                try:
                    write_json(self.work / "run_manifest.json", self.manifest)
                    self.log("SUMMARY", json.dumps({"exit_code": exit_code, "stages": self.manifest["stages"],
                        "outputs": self.manifest["outputs"]}))
                    if exit_code == 0 and not self.options.keep_intermediate and not self.options.check:
                        # Keep the manifest and diagnostics; remove only known generated runner products.
                        for path in self.work.glob("stage*.csv"):
                            path.unlink(missing_ok=True)
                        for path in self.work.glob("stage*.vtu"):
                            path.unlink(missing_ok=True)
                except OSError as error:
                    exit_code = self._failure(RunError("E601", f"Cannot finalize output files: {error}"))
                    self.manifest["exit_code"] = exit_code
                    try:
                        write_json(self.work / "run_manifest.json", self.manifest)
                    except OSError:
                        pass
                finally:
                    self.lock.release()
            self._restore_signals()
            if self.stream:
                self.stream.close()
        return exit_code

    def _failure(self, error):
        report = {"stage": self.stage, "message_id": error.message_id, "message": error.message,
                  "details": error.details, "exit_code": error.exit_code}
        try:
            self.log("FATAL", str(error))
        except OSError:
            pass
        print(str(error), file=sys.stderr)
        if self.lock.owned:
            try:
                write_json(self.work / "fail_report.json", report)
            except OSError:
                pass
        return error.exit_code

    def _write_plan(self):
        for index, deck in enumerate(self.result.decks):
            stage = self.result.stage_plan["stages"][index]
            metadata = {"CmdSHA256": self.ir["source"]["sha256"], "Stage": index}
            swept = deck["sweep"]["contact"]
            for name, bias in stage["electrode_voltages_start_V"].items():
                if name != swept:
                    metadata[name + "_OuterVoltage_V"] = bias
                    metadata[name + "_InnerVoltage_V"] = bias
            deck["sweep"]["vtk_field_data"] = metadata
            deck.setdefault("runtime_log", {})["append"] = True
        write_stage_plan(self.work, self.result)
        self.manifest["decisions"] = self.result.decisions
        self.manifest["sidecar_overrides"] = [{"stage": s["index"], "keys": s["sidecar_overrides"]}
                                             for s in self.result.stage_plan["stages"]]
        for decision in self.result.decisions:
            self.log("DECISION", json.dumps(decision))

    def _execute(self, index, runner):
        stage = self.result.stage_plan["stages"][index]
        for path in [self.work / stage["csv"], *self.work.glob(f"stage{index:02d}_*.vtu")]:
            path.unlink(missing_ok=True)
        if stage["output_state"] != stage["input_state"]:
            (self.work / stage["output_state"]).unlink(missing_ok=True)
        env = dict(os.environ)
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "EIGEN_NUM_THREADS"):
            env[key] = str(self.result.stage_plan["threads"])
        env["HDF5_USE_FILE_LOCKING"] = "FALSE"
        profile = "debug" if self.options.verbose else self.result.decks[index].get("runtime_log", {}).get("profile", "default")
        # Keep runner arguments ASCII. Its narrow filesystem API cannot accept
        # Unicode Windows paths; Python publishes the user's log and filenames.
        runner_log = self.work / f"stage{index:02d}_runner.log"
        runner_log.unlink(missing_ok=True)
        args = runner + ["--config", stage["deck"], "--log", runner_log.name, "--log-profile", profile]
        try:
            self.proc = subprocess.Popen(args, cwd=self.work, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace", creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if os.name == "nt" else 0)
        except OSError as error:
            raise RunError("E503", f"Solver could not start: {error}") from error
        lines, errors = queue.Queue(), queue.Queue()
        stop = threading.Event()
        progress_count = [0]
        def read_output():
            for line in self.proc.stdout:
                lines.put(line.rstrip())
        def publish():
            seen = 0
            while not stop.wait(.05):
                try:
                    records = self.publisher.update(index)
                    if not self.options.quiet:
                        for row in records[seen:]:
                            print(f"stage={index} bias={row['bias_V']} V accepted", flush=True)
                    seen = len(records)
                    progress_count[0] = seen
                except (FileNotFoundError, ValueError, PermissionError):
                    continue  # wait for complete CSV/file publication; final read is strict
                except Exception as error:
                    errors.put(error)
                    return
        reader = threading.Thread(target=read_output, daemon=True)
        writer = threading.Thread(target=publish, daemon=True)
        reader.start(); writer.start()
        captured = []
        try:
            while self.proc.poll() is None:
                if self.interrupted or not errors.empty():
                    self._stop_child()
                    break
                try:
                    line = lines.get(timeout=.05)
                    captured.append(line)
                    self.log("RUNNER", line)
                except queue.Empty:
                    pass
        finally:
            stop.set()
            # A slow filesystem can keep update() busy beyond five seconds.
            # Final publication must never overlap this publisher's mutation
            # of PLT/PVD temporary files and its accepted-record bookkeeping.
            writer.join()
            reader.join()
            self.proc.stdout.close()
        while not lines.empty():
            line = lines.get_nowait(); captured.append(line); self.log("RUNNER", line)
        for line in captured:
            if "vtk_fields_skipped" in line:
                self.log("WARNING", "W202: " + line)
        if runner_log.exists():
            for line in runner_log.read_text(encoding="utf-8", errors="replace").splitlines():
                self.log("RUNTIME", line)
                if "vtk_fields_skipped" in line:
                    self.log("WARNING", "W202: " + line)
        if not errors.empty():
            raise RunError("E602", f"Output writer failed: {errors.get()}")
        try:
            records = self.publisher.update(index, final=True)
        except (OSError, ValueError) as error:
            raise RunError("E602", f"Output writer failed: {error}") from error
        if not self.options.quiet:
            for row in records[progress_count[0]:]:
                print(f"stage={index} bias={row['bias_V']} V accepted", flush=True)
        return self.proc.returncode, records, "\n".join(captured)

    def _reached(self, index, records):
        stage = self.result.stage_plan["stages"][index]
        voltages = dict(stage["electrode_voltages_start_V"])
        if not records:
            return False, voltages
        sweep = stage.get("sweep", {})
        contact = self.result.decks[index]["sweep"]["contact"]
        last = records[-1]
        voltages[contact] = float(last.get("inner_voltage_V") or last["bias_V"])
        if sweep.get("mode") == "current":
            return abs(float(last["current_total_A_per_um"]) - sweep["stop_current_A_per_um"]) <= 1e-10, voltages
        outer = float(last.get("outer_voltage_V") or last["bias_V"])
        return abs(outer - stage["electrode_voltages_goal_V"][contact]) <= 1e-12, voltages
