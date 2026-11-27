"""SDevice frontend gates and stage execution (WP5)."""
from __future__ import annotations
import csv
import json
from decimal import Decimal
from hashlib import sha256
import math
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
from ..cmd_parser.syntax import SourceNumber
from ..deck_emit import DeckEmitError, EmitOptions, bind_current_start, emit_stage_plan, rebase_stage, write_stage_plan
from ..materials import MaterialsError, active_models_from_ir, emit_materials, write_materials
from ..log_levels import WORKBENCH_COMMAND_ENV, log_enabled
from ..mesh_import import MeshImportError, import_grid, check_electrodes
from ..mesh_import.doping import DopingPolicy
from ..resources import resource_root
from .. import schema_validator as sv
from .outputs import OutputPublisher, accepted_records
from .extended_stage import control_decimal, read_accepted_extended_stage, verify_extended_stage_parent
from .arithmetic_input import collect_effective_doping, collect_input, collect_slot_input
from .slot_ialmob import bind_slot_stage_inputs, check_slot_effective_doping, stage_slot_decimal_material
from .support import RunError, RunLock, write_json


def _is_slot_idvd_profile(ir, mesh):
    keys = {entry["key"] for entry in ir["physics"]["effective"]}
    goals = [stage["goal"] for stage in ir.get("stages", []) if stage.get("goal")]
    idsat = ([goal["contact"] for goal in goals] == ["drain", "gate"]
             and all("voltage_V" in goal for goal in goals))
    return (("Recombination.Avalanche" in keys or idsat)
            and {"Fermi", "Mobility.PhuMob", "Mobility.Enormal", "Mobility.HighFieldSaturation",
                 "Recombination.SRH", "Recombination.Auger"} <= keys
            and {"Si", "PolySilicon"} <= set(mesh["materials"])
            and {e["name"] for e in ir["electrodes"]} == {"source", "drain", "gate", "substrate"})


def _slot_idvd_ohmic_polarity(work, contacts):
    ohmic = {contact["name"] for contact in contacts if contact.get("type") == "ohmic"}
    mesh = json.loads((work / "mesh.json").read_text(encoding="utf-8"))
    nodes = {contact["name"]: contact["node_ids"] for contact in mesh["contacts"]
             if contact["name"] in ohmic}
    with (work / "doping.csv").open(encoding="utf-8", newline="") as source:
        signs = {}
        for row in csv.DictReader(source):
            donors, acceptors = Decimal(row["donors_cm3"]), Decimal(row["acceptors_cm3"])
            signs[int(row["node_id"])] = int(donors > acceptors) - int(donors < acceptors)
    result = {}
    for name, ids in nodes.items():
        polarity = [signs[node] for node in ids]
        result[name] = {"positive": polarity.count(1), "negative": polarity.count(-1),
                        "zero": polarity.count(0)}
    return result


def _apply_slot_idvd_bgn_default(ir):
    """Resolve the native SLOT IdVd default before PAR projection and deck emission."""
    effective = ir["physics"]["effective"]
    selectors = {entry["key"] for entry in effective
                 if entry["key"].startswith("EffectiveIntrinsicDensity.")}
    if selectors - {"EffectiveIntrinsicDensity.NoFermi"}:
        return False
    effective.append({
        "key": "EffectiveIntrinsicDensity.BennettWilson",
        "parameters": {},
        "status": "supported",
        "vela": "solver.bandgap_narrowing=bennett_wilson",
    })
    effective.sort(key=lambda entry: entry["key"])
    return True

def executable(name, bin_dir):
    variable = {
        "vela_example_runner": "VELA_EXAMPLE_RUNNER",
        "vela_example_runner_extended80": "VELA_EXAMPLE_RUNNER_EXTENDED80",
        "sentaurus_import": "VELA_SENTAURUS_IMPORT",
    }[name]
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


def arithmetic_capabilities(runner, cwd):
    try:
        done = subprocess.run(runner + ["--arithmetic-capabilities"], cwd=cwd,
                              capture_output=True, encoding="utf-8", timeout=10, check=True)
        capabilities = json.loads(done.stdout)
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        raise RunError("E202", "Selected runner did not provide arithmetic capabilities") from error
    if not isinstance(capabilities, dict):
        raise RunError("E202", "Selected runner returned invalid arithmetic capabilities")
    return capabilities


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
        source_arithmetic = {}
        self.ir = parse_cmd_file(self.cmd, arithmetic_input=source_arithmetic)
        self.source_temperature = source_arithmetic.get("temperature_K")
        self.arithmetic_input = source_arithmetic if getattr(options, "check_arithmetic_input", False) else None
        stem = self.cmd.stem
        if not re.fullmatch(r"[A-Za-z0-9._-]+", stem):
            stem = "device_" + self.ir["source"]["sha256"][:8]
        self.work = self.cmd.parent / (stem + "_vela")
        self.log_path = self.cmd.parent / self.ir["files"]["output"]["resolved"]
        self.invocation = os.environ.get(WORKBENCH_COMMAND_ENV)
        self.lock = RunLock(self.work / ".lock")
        self.stream = None
        self.stage = None
        self.proc = None
        self.interrupted = False
        self.signal_number = None
        self.handlers = {}
        self.result = None
        self.publisher = None
        self.extended_states = {}
        self.manifest = {"schema": "vela.sdevice_run_manifest.v1", "cmd": str(self.cmd),
            "sha256": self.ir["source"]["sha256"], "version": version_info(), "check_only": options.check,
            "stages": [], "outputs": [], "decisions": [], "sidecar_overrides": [], "exit_code": None}

    def log(self, level, text):
        if self.stream and log_enabled(self.options.log_level, level, text):
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
            if self.invocation is not None:
                self.log("INVOCATION", self.invocation)
            self.log("CMD", "Command file echo begins")
            if log_enabled(self.options.log_level, "CMD"):
                self.stream.write(self.cmd.read_text(encoding="utf-8-sig", errors="replace") + "\n")
            self.log("CMD", "Command file echo ends")
            if self.lock.stale:
                owner = self.lock.previous
                self.log("WARNING", f"W701: Stale lock from PID {owner['pid']} on {owner['host']} replaced")
            for warning in self.ir["warnings"]:
                self.log("WARNING", json.dumps(warning))
            sidecar_path = self.cmd.with_suffix(".vela.json")
            slot_context = slot_report = None
            try:
                sidecar = (json.loads(sidecar_path.read_text(encoding="utf-8"), parse_float=SourceNumber)
                           if sidecar_path.exists() else None)
                if isinstance(sidecar, dict) and "slot_ialmob" in sidecar:
                    sv.validate(sv.load_schema(resource_root() / "schemas/vela.sdevice_sidecar.v1.schema.json"), sidecar)
                    if set(sidecar) - {"schema", "description", "slot_ialmob"} or getattr(self.options, "pn2d_population", False):
                        raise ValueError("SLOT binding cannot be combined with other overrides or PN2D")
                    self.ir, slot_report, slot_context = bind_slot_stage_inputs(self.ir, self.cmd, sidecar["slot_ialmob"])
                    write_json(self.work / "slot_ialmob_binding.json", slot_report)
                    self.manifest["slot_ialmob_binding"] = "slot_ialmob_binding.json"
            except (OSError, ValueError, KeyError, TypeError, sv.SchemaError, sv.ValidationError) as error:
                raise RunError("E409", f"Sidecar '{sidecar_path}' is invalid: {error}") from error
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
            precision = self.ir["math"].get("arithmetic_precision", "double")
            if self.arithmetic_input is not None and precision != "extended80":
                raise RunError("E202", "--check-arithmetic-input requires ExtendedPrecision in Math")
            runner_name = "vela_example_runner_extended80" if precision == "extended80" else "vela_example_runner"
            self.manifest["arithmetic"] = {"requested": precision, "runner": runner_name,
                "capability_check": "not_evaluated",
                "linear_solver_override": os.environ.get("VELA_LINEAR_SOLVER") or None}
            try:
                runner = (None if self.options.check and self.arithmetic_input is None
                          else executable(runner_name, location))
            except RunError:
                self.manifest["arithmetic"]["capability_check"] = "failed"
                raise
            importer = executable("sentaurus_import", location)
            grid_path = self.cmd.parent / self.ir["files"]["grid"]["resolved"]
            if doping_file := self.ir["files"].get("doping"):
                declared = self.cmd.parent / doping_file["resolved"]
                sources = {grid_path.resolve(), grid_path.with_suffix(".dat").resolve()} if grid_path.suffix.lower() == ".grd" else {grid_path.resolve()}
                if declared.resolve() not in sources:
                    raise RunError("E402", "File{Doping} names a separate mesh input; only Grid or its DF-ISE DAT is supported")
            registry = list(getattr(self.options, "datexcodes", ()))
            local_registry = self.cmd.parent / "datexcodes.txt"
            if local_registry.is_file() and local_registry.resolve() not in {path.resolve() for path in registry}:
                registry.append(local_registry)
            mode = getattr(self.options, "doping_mode", "auto")
            compute_net = self.ir["math"].get("compute_doping_concentration", False)
            if compute_net and mode != "auto":
                raise RunError("E402", "ComputeDopingConcentration requires --doping-mode auto")
            try:
                policy = DopingPolicy.from_inputs(mode, compute_net, registry if mode == "auto" else ())
            except (OSError, ValueError) as error:
                raise RunError("E402", f"Cannot load datexcodes: {error}") from error
            decimal_mesh = precision == "extended80" and slot_context is None
            mesh = import_grid(grid_path, self.work,
                               importer=importer[0] if importer else None,
                               decimal_doping=decimal_mesh, decimal_coordinates=decimal_mesh, doping_policy=policy)
            findings = [f for f in check_electrodes(mesh, [e["name"] for e in self.ir["electrodes"]]) if not f["matched"]]
            if findings:
                raise RunError("E401", "Electrodes have no matching mesh contacts", findings=findings)
            sv.validate(sv.load_schema(resource_root() / "schemas/vela.mesh_qualification.v1.schema.json"), mesh)
            constants_profile = getattr(self.options, "physical_constants_profile", None) or "sentaurus_u2022_12"
            electrode_materials = [e["attributes"]["material"] for e in self.ir["electrodes"] if "material" in e["attributes"]]
            default_par = getattr(self.options, "default_par", None)
            if slot_context:
                if default_par is not None:
                    raise RunError("E407", "--default-par applies to the generic material importer")
                constants_profile = slot_context["reference"].get("physical_constants_profile", "codata2018")
                slot_report["effective_doping"] = check_slot_effective_doping(
                    self.work / "doping.csv", Path(slot_report["inputs"]["doping"]["path"]))
                for entry in slot_report["inputs"].values():
                    if "copy_as" in entry:
                        target = self.work / entry["copy_as"]
                        shutil.copyfile(entry["path"], target)
                        if sha256(target.read_bytes()).hexdigest() != entry["sha256"]:
                            raise RunError("E409", f"SLOT input changed during copy: {target.name}")
                if precision == "extended80":
                    stage_slot_decimal_material(self.work / "materials.json", slot_report,
                                                slot_context["reference"], self.source_temperature)
                write_json(self.work / "slot_ialmob_binding.json", slot_report)
                material_report = {"schema": "vela.sdevice_materials_report.v1", "lossy": [], "inert": []}
            else:
                if default_par is not None and not default_par.is_file():
                    raise RunError("E303", f"Material defaults file '{default_par}' is unavailable")
                slot_idvd = _is_slot_idvd_profile(self.ir, mesh)
                if slot_idvd and _apply_slot_idvd_bgn_default(self.ir):
                    write_json(self.work / "execution_ir.json", self.ir)
                    self.manifest["bandgap_narrowing_default"] = "bennett_wilson"
                    self.log("PHYSICS", "Bennett–Wilson selected by the SLOT Fermi/PhuMob profile")
                materials = emit_materials(self.cmd.parent / self.ir["files"]["parameter"]["resolved"],
                    default_par=default_par,
                    constants_profile=constants_profile,
                    active_models=active_models_from_ir(self.ir), temperature_K=self.ir["physics"]["temperature_K"] or 300,
                    mesh_materials=mesh["materials"], electrode_materials=electrode_materials,
                    effective_models=self.ir["physics"]["effective"], allow_lossy=self.options.allow_lossy,
                    band_edge_sensitive=any(e["attributes"].get("schottky") for e in self.ir["electrodes"]),
                    capture_arithmetic_inputs=self.arithmetic_input is not None,
                    extended80_temperature_literal=self.source_temperature if precision == "extended80" else None)
                if slot_idvd and precision == "extended80":
                    materials.report["slot_idvd_population_coordinates"] = True
                write_materials(self.work, materials)
                material_report = materials.report
            self.manifest["doping"] = mesh["doping"]["resolution"]
            write_json(self.work / "mesh_qualification.json", mesh)
            for kind in ("lossy", "inert"):
                for parameter in material_report[kind]:
                    self.log("MATERIAL", json.dumps({"status": kind, **parameter}))
            self.result = emit_stage_plan(self.ir, mesh_qualification=mesh, materials_report=material_report,
                work_dir=self.work.name, cmd_basename=self.cmd.name, cmd_sha256=self.ir["source"]["sha256"],
                options=EmitOptions(exit_on_failure=self.options.exit_on_failure, threads=self.options.threads,
                    plot_final_only=self.options.plot_final_only, vtu_binary=self.options.vtu_binary,
                    keep_intermediate=self.options.keep_intermediate, allow_lossy=self.options.allow_lossy,
                    defer_current_times=True, pn2d_population=getattr(self.options, "pn2d_population", False)),
                sidecar=sidecar, slot_context=slot_context)
            if material_report.get("slot_idvd_population_coordinates"):
                contacts = self.result.decks[0]["contacts"]
                audit = _slot_idvd_ohmic_polarity(self.work, contacts)
                self.manifest["ohmic_contact_doping"] = audit
                self.log("PHYSICS", f"Node-local Ohmic contact doping counts: {json.dumps(audit)}")
            self._write_plan()
            self.log("INFO", "Multi-stage PLT time resets per stage; sign factor +1 (WP10 pending)")
            if self.interrupted:
                raise RunError("E801", f"Interrupted by signal {self.signal_number} during input gates")
            if self.options.check:
                if self.arithmetic_input is None:
                    self.log("INFO", "Check-only gates passed; runner not invoked")
                else:
                    self._check_arithmetic_input(runner, sidecar_path, constants_profile, material_report,
                                                 slot_report, mesh["doping"]["source_fields"])
            else:
                if runner is None:
                    self.manifest["arithmetic"]["capability_check"] = "failed"
                    raise RunError("E202", f"Runner '{runner_name}' for {precision} arithmetic cannot be resolved")
                if precision == "extended80":
                    self._require_extended80(runner)
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
                        exact = self.extended_states.get(index - 1)
                        if exact is not None:
                            current = exact.currents.get(contact)
                        else:
                            current = records[-1].get(f"contact_{contact}_total_A_per_um")
                            if current is None and records[-1].get("current_contact") == contact:
                                current = records[-1].get("current_total_A_per_um")
                            if current is not None:
                                current = float(current)
                        if current is None:
                            raise RunError("E407", f"Previous accepted state has no measured current for '{contact}'")
                        self.result = bind_current_start(self.result, index, current, voltages)
                        self.publisher.result = self.result
                        self._write_plan()
                        if precision == "extended80":
                            self._require_extended80(runner)
                        stage = self.result.stage_plan["stages"][index]
                        self.log("INFO", f"Current Time mapping bound to {current} A/um at stage {index}")
                    self.log("STAGE", f"Begin stage {index}: {stage['kind']}")
                    status, records, stderr, reported_precision = self._execute(index, runner)
                    reached, voltages = self._reached(index, records)
                    stopped = self._break_reached(index, records)
                    runtime = {"index": index, "runner_exit_code": status, "goal_reached": reached,
                               "accepted_points": len(records), "reached_voltages_V": voltages,
                               "requested_arithmetic_precision": precision,
                               "reported_arithmetic_precision": reported_precision}
                    if index in self.extended_states:
                        runtime["accepted_extended80_state"] = self.extended_states[index].record
                    if stopped:
                        runtime["stop_criterion_reached"] = stopped
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
                    if (status != 0 or not (reached or stopped)) and self.result.stage_plan["exit_on_failure"]:
                        raise RunError("E502", f"Stage {index} failed with ExitOnFailure", stderr=stderr)
                    if not records or not (self.work / stage["output_state"]).is_file():
                        raise RunError("E503", f"Stage {index} produced no accepted state for the next stage", stderr=stderr)
                    if index + 1 < count and (not reached or stage.get("fixed_contact_current")
                                               or stage.get("sweep", {}).get("mode") == "current"):
                        self.result = rebase_stage(self.result, index + 1, voltages, update_saved_state=True)
                        self.publisher.result = self.result
                        self.log("WARNING", f"W206: stage {index} reached {voltages}; rebasing later decks")
                        self._write_plan()
                        if precision == "extended80":
                            self._require_extended80(runner)
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

    def _require_extended80(self, runner):
        audit = self.manifest["arithmetic"]
        audit["capability_check"] = "failed"
        audit["executable"] = runner[-1]
        capabilities = arithmetic_capabilities(runner, self.work)
        audit["reported_capabilities"] = capabilities
        if (capabilities.get("arithmetic_precision"), capabilities.get("radix"),
                capabilities.get("significand_bits"), capabilities.get("max_exponent")) != (
                "extended80", 2, 64, 16384):
            raise RunError("E202", "Selected runner does not have the requested 80-bit arithmetic format")
        override = audit["linear_solver_override"]
        coverage = audit["frontend_input_check"]
        backends = capabilities.get("end_to_end_linear_backends")
        consumers = capabilities.get("end_to_end_input_consumers", [])
        scoped = False
        if capabilities.get("scope") != "end_to_end":
            profiles = capabilities.get("scoped_extended80_profiles", {})
            requested_profiles = {stage.get("profile") for stage in coverage["stages"]}
            if len(requested_profiles) == 1 and None not in requested_profiles and isinstance(profiles, dict):
                profile = profiles.get(next(iter(requested_profiles)))
                if isinstance(profile, dict):
                    scoped = True
                    backends = profile.get("linear_backends")
                    consumers = profile.get("input_consumers")
        if not isinstance(backends, list) or not all(isinstance(backend, str) for backend in backends):
            raise RunError("E202", "Selected runner did not identify its end-to-end linear backends")
        if override is not None and override not in backends:
            raise RunError("E202", f"VELA_LINEAR_SOLVER={override} is unsupported for end-to-end extended80")
        if capabilities.get("scope") != "end_to_end" and not scoped:
            raise RunError("E202", "ExtendedPrecision requires end-to-end extended80; "
                           f"the selected runner reports {capabilities.get('scope', 'unknown')} precision")
        if not isinstance(consumers, list) or not all(isinstance(item, str) for item in consumers):
            raise RunError("E202", "Selected runner did not identify its end-to-end input consumers")
        supported = set(consumers)
        for stage in coverage["stages"]:
            stage["binary64_consumers"] = [
                item for item in stage["required_consumers"] if item not in supported
            ]
        if any(stage["binary64_consumers"] or stage["missing_inputs"] for stage in coverage["stages"]):
            raise RunError("E202", "Frontend input conversion is binary64 on the selected path",
                           input_coverage=coverage)
        if not backends:
            raise RunError("E202", "Selected runner has no end-to-end extended80 linear backend")
        coverage["status"] = "passed"
        audit["capability_check"] = "passed"

    def _check_arithmetic_input(self, runner, sidecar_path, profile, materials, slot_binding, doping_fields):
        if runner is None:
            raise RunError("E202", "Extended80 runner is required for the arithmetic input check")
        source = collect_input(self.cmd, self.arithmetic_input,
            self.cmd.parent / self.ir["files"]["parameter"]["resolved"],
            self.cmd.parent / self.ir["files"]["grid"]["resolved"], sidecar_path, profile,
            tdr_export=self.work / "tdr_export", doping_fields=doping_fields)
        if slot_binding is None:
            source["materials"] = {row["name"]: row["arithmetic_input"] for row in materials["materials"]
                                   if "arithmetic_input" in row}
            source["sources"]["effective_doping"] = collect_effective_doping(self.work / "doping.csv")
        else:
            stage_decks = [self.work / stage["deck"] for stage in self.result.stage_plan["stages"]]
            reference = Path(slot_binding["inputs"]["deck"]["path"])
            sources, source["materials"] = collect_slot_input(self.work, slot_binding, reference, stage_decks)
            source["sources"].update(sources)
        input_path = self.work / "arithmetic_input.json"
        write_json(input_path, source)
        try:
            result = subprocess.run(runner + ["--check-arithmetic-input", input_path.name],
                                    cwd=self.work, capture_output=True, encoding="utf-8", timeout=90)
            if result.returncode:
                raise RunError("E202", "C++ arithmetic input check failed", stderr=result.stderr)
            report = json.loads(result.stdout)
            if (report["arithmetic_precision"], report["radix"], report["significand_bits"],
                    report["max_exponent"], report["simulation_started"]) != ("extended80", 2, 64, 16384, False):
                raise RunError("E202", "Arithmetic input check did not use the extended80 scalar in isolation")
            summary = {role: value["number_count"] for role, value in report["sources"].items()}
            if summary != {role: len(value["numbers"]) for role, value in source["sources"].items()} or any(
                    report["sources"][role]["sha256"] != value["sha256"] for role, value in source["sources"].items()):
                raise RunError("E202", "Arithmetic input report identity does not match its inputs")
            if set(report.get("materials", {})) != set(source["materials"]):
                raise RunError("E202", "Arithmetic input report does not cover the selected materials")
        except (OSError, subprocess.SubprocessError, ValueError, TypeError, KeyError) as error:
            raise RunError("E202", f"Cannot verify C++ arithmetic input: {error}") from error
        if slot_binding is not None:
            report["slot_binding"] = {
                "report": "slot_ialmob_binding.json",
                "materials_recomputed": sorted(source["materials"]),
                "scope": "input_audit_only",
            }
        write_json(self.work / "arithmetic_input_report.json", report)
        self.manifest["arithmetic"]["input_check"] = {"format": "extended80", "counts": summary,
            "report": "arithmetic_input_report.json", "simulation_started": False}
        self.log("INFO", f"Extended80 input check completed for {summary}")

    def _write_plan(self):
        self.manifest["arithmetic"]["planned_stages"] = [
            {"index": index, "requested": deck["solver"].get("arithmetic_precision", "double"),
             "require_end_to_end": deck["solver"].get("require_end_to_end_precision", False)}
            for index, deck in enumerate(self.result.decks)
        ]
        if self.manifest["arithmetic"]["requested"] == "extended80":
            stages = []
            selected_decks = self.result.decks[:1] if self.options.initial_only else self.result.decks
            for index, deck in enumerate(selected_decks):
                solver = deck["solver"]
                required = ["mesh", "doping", "materials", "temperature", "contacts", "sweep"]
                if solver.get("mobility", {}).get("model") == "ialmob":
                    required.append("ialmob")
                elif solver.get("extended80_stage_profile") == "slot_population_phumob_fixed_bias_v1":
                    required.append("phumob")
                if deck["sweep"].get("initial_state_file") or deck["sweep"].get("write_state_file"):
                    required.append("state_archive")
                missing = []
                control = deck.get("decimal_stage_controls")
                if not control or not all(key in control for key in ("contacts", "sweep", "solver")):
                    missing.append("decimal_stage_controls")
                elif "phumob" in required and not control["solver"].get("slot_source_coefficients"):
                    missing.append("slot_source_coefficients")
                stages.append({"index": index, "profile": solver.get("extended80_stage_profile"), "required_consumers": required,
                               "binary64_consumers": list(required), "missing_inputs": missing})
            self.manifest["arithmetic"]["frontend_input_check"] = {"status": "blocked", "stages": stages}
        for index, deck in enumerate(self.result.decks):
            stage = self.result.stage_plan["stages"][index]
            metadata = {"CmdSHA256": self.ir["source"]["sha256"], "Stage": index}
            swept = deck["sweep"]["contact"]
            held = stage.get("fixed_contact_current", {}).get("contact")
            for name, bias in stage["electrode_voltages_start_V"].items():
                if name not in (swept, held):
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
        strict = self.result.decks[index]["solver"].get("require_end_to_end_precision", False)
        if strict and stage["input_state"]:
            archive = stage["input_state"]
            for prior in self.result.stage_plan["stages"][:index]:
                if prior["output_state"] == archive:
                    verify_extended_stage_parent(self.work, archive, self.extended_states.get(prior["index"]))
                    break
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
        report_progress = not self.options.quiet or self.invocation is not None

        def read_output():
            for line in self.proc.stdout:
                lines.put(("RUNNER", line.rstrip()))

        def publish():
            seen = 0
            while not stop.wait(.05):
                try:
                    records = self.publisher.update(index)
                    if report_progress:
                        for row in records[seen:]:
                            lines.put(("PROGRESS", f"stage={index} bias={row['bias_V']} V accepted"))
                    seen = len(records)
                    progress_count[0] = seen
                except (FileNotFoundError, ValueError, PermissionError):
                    continue  # wait for complete CSV/file publication; final read is strict
                except Exception as error:
                    errors.put(error)
                    return
        reader = threading.Thread(target=read_output, daemon=True)
        writer = None if strict else threading.Thread(target=publish, daemon=True)
        reader.start()
        if writer:
            writer.start()
        captured = []

        def record(category, line):
            if category == "PROGRESS":
                if not self.options.quiet:
                    print(line, flush=True)
                if self.invocation is not None:
                    self.log(category, line)
            else:
                captured.append(line)
                self.log(category, line)

        try:
            while self.proc.poll() is None:
                if self.interrupted or not errors.empty():
                    self._stop_child()
                    break
                try:
                    record(*lines.get(timeout=.05))
                except queue.Empty:
                    pass
        finally:
            stop.set()
            # A slow filesystem can keep update() busy beyond five seconds.
            # Final publication must never overlap this publisher's mutation
            # of PLT/PVD temporary files and its accepted-record bookkeeping.
            if writer:
                writer.join()
            reader.join()
            self.proc.stdout.close()
        while not lines.empty():
            record(*lines.get_nowait())
        reported_precision = None
        accepted_summary = {}
        for line in captured:
            if "vtk_fields_skipped" in line:
                self.log("WARNING", "W202: " + line)
            if line.startswith("{"):
                try:
                    summary = json.loads(line)
                except ValueError:
                    continue
                if isinstance(summary, dict) and summary.get("config") == stage["deck"]:
                    reported_precision = summary.get("arithmetic_precision")
                    accepted_summary = summary
        if runner_log.exists():
            for line in runner_log.read_text(encoding="utf-8", errors="replace").splitlines():
                self.log("RUNTIME", line)
                if "vtk_fields_skipped" in line:
                    self.log("WARNING", "W202: " + line)
        if not errors.empty():
            raise RunError("E602", f"Output writer failed: {errors.get()}")
        if strict and not self.interrupted:
            records = accepted_records(self.work / stage["csv"])
            if records:
                self.extended_states[index] = read_accepted_extended_stage(
                    self.work, stage, accepted_summary, records, self.result.decks[index]["sweep"]["contact"],
                    runner)
            elif self.proc.returncode == 0:
                raise RunError("E202", "Extended80 stage reported success without accepted state data")
        if strict and index not in self.extended_states:
            return self.proc.returncode, [], "\n".join(captured), reported_precision
        try:
            records = self.publisher.update(index, final=True)
        except (OSError, ValueError) as error:
            raise RunError("E602", f"Output writer failed: {error}") from error
        if report_progress:
            for row in records[progress_count[0]:]:
                record("PROGRESS", f"stage={index} bias={row['bias_V']} V accepted")
        return self.proc.returncode, records, "\n".join(captured), reported_precision

    def _reached(self, index, records):
        stage = self.result.stage_plan["stages"][index]
        voltages = dict(stage["electrode_voltages_start_V"])
        if not records:
            return False, voltages
        sweep = stage.get("sweep", {})
        contact = self.result.decks[index]["sweep"]["contact"]
        exact = self.extended_states.get(index)
        if exact is not None:
            if sweep.get("mode") == "current":
                reached = abs(control_decimal(exact.currents[contact]) -
                              control_decimal(sweep["stop_current_A_per_um"])) <= Decimal("1e-10")
            else:
                reached = abs(control_decimal(exact.outer) -
                              control_decimal(stage["electrode_voltages_goal_V"][contact])) <= Decimal("1e-12")
            return reached, exact.voltages.copy()
        last = records[-1]
        voltages[contact] = float(last.get("inner_voltage_V") or last["bias_V"])
        if fixed := stage.get("fixed_contact_current"):
            held = fixed["contact"]
            voltages[held] = float(last[f"contact_{held}_voltage_V"])
        if sweep.get("mode") == "current":
            return abs(float(last["current_total_A_per_um"]) - sweep["stop_current_A_per_um"]) <= 1e-10, voltages
        outer = float(last.get("outer_voltage_V") or last["bias_V"])
        return abs(outer - stage["electrode_voltages_goal_V"][contact]) <= 1e-12, voltages

    def _break_reached(self, index, records):
        """A configured physical stop completes a stage before the sweep goal."""
        if not records:
            return None
        last = records[-1]
        criteria = self.result.stage_plan["stages"][index].get("break", {})
        if ionization := criteria.get("ionization_integral"):
            event = last.get("ionization_stop_criterion", "")
            if event.startswith(f"path_ionization_integral_rank_{ionization['rank']}_"):
                return "ionization_integral"
        for criterion in criteria.get("current", []):
            contact = criterion["contact"]
            exact = self.extended_states.get(index)
            if exact is not None:
                value = control_decimal(exact.currents[contact])
                threshold = lambda key: control_decimal(criterion[key])
            else:
                value = last.get(f"contact_{contact}_total_A_per_um")
                if value in (None, "") and last.get("current_contact") == contact:
                    value = last.get("current_total_A_per_um")
                if value in (None, "") or not math.isfinite(float(value)):
                    continue
                value = float(value)
                threshold = criterion.__getitem__
            if ("abs_val_A_per_um" in criterion and abs(value) >= threshold("abs_val_A_per_um") or
                    "min_val_A_per_um" in criterion and value <= threshold("min_val_A_per_um") or
                    "max_val_A_per_um" in criterion and value >= threshold("max_val_A_per_um")):
                return "current"
        return None
