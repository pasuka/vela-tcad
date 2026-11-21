"""Prepare frozen workbench cases and start the packaged solver."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

from vela_sdevice.log_levels import WORKBENCH_COMMAND_ENV


ACCEPTED_POINT = re.compile(r"stage=(\d+) bias=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?) V accepted")

@dataclass(frozen=True)
class Case:
    name: str
    source: Path
    arguments: list[str]
    hashes: dict[str, str]
    work: str
    log: str
    curve: str | None = None
    field: str | None = None
    original_script: str | None = None
    linear_solver: str = "sparselu"
    stop_criterion: str | None = None
    ctest: str | None = None
    name_zh: str | None = None


def _hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def prepare_cases(destination):
    """Freeze the inputs prepared by the numerical sd-sim regressions."""
    from vela_sdevice.cmd_parser import parse_cmd_file
    from vela_workbench_catalog import stage_catalog

    if destination.exists():
        raise FileExistsError(f"Case catalog already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temporary:
        staging = Path(temporary)
        catalog = staging / "cases"
        catalog.mkdir()
        cases = stage_catalog(catalog, staging)
        for key, case in cases.items():
            directory = catalog / key
            ir = parse_cmd_file(directory / case["arguments"][-1])
            log = Path(ir["files"]["output"]["resolved"])
            if log.anchor or ".." in log.parts:
                raise ValueError(f"Case {key} must use a relative log path within the run directory")
            case["log"] = log.as_posix()
            case["hashes"] = {
                path.name: _hash(path) for path in sorted(directory.iterdir())
                if path.name != case["original_script"]
            }
        (catalog / "cases.json").write_text(
            json.dumps(cases, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        catalog.rename(destination)


def read_cases(directory):
    manifest = json.loads((directory / "cases.json").read_text(encoding="utf-8"))
    return {key: Case(source=directory / key, **value) for key, value in manifest.items()}


def verify_inputs(case, directory):
    for name, checksum in case.hashes.items():
        if _hash(directory / name) != checksum:
            raise ValueError(f"Input checksum mismatch: {name}")


def write_run_failure(path, invocation, error, output=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as log:
        if log.tell() == 0:
            log.write(f"[INVOCATION] {invocation}\n")
        log.write(f"[FATAL] {error}\n")
        for line in output:
            log.write(f"[CONSOLE] {line}\n")


def run_case(launcher, case, directory, *, log_level="debug", on_progress=None):
    """Stream solver progress while sd-sim writes the case's main log."""
    env = dict(os.environ, VELA_LINEAR_SOLVER=case.linear_solver, VELA_LINEAR_THREADS="1")
    if case.linear_solver == "pardiso_mkl":
        env.update(MKL_THREADING_LAYER="SEQUENTIAL", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
        env.pop("MKL_CBWR", None)
    command = [str(launcher)]
    if log_level != "debug":
        command.extend(("--log-level", log_level))
    command.extend(case.arguments)
    invocation = subprocess.list2cmdline(command) if os.name == "nt" else shlex.join(command)
    env[WORKBENCH_COMMAND_ENV] = invocation
    log_path = directory / case.log
    console_tail = deque(maxlen=20)
    try:
        for name in case.hashes:
            shutil.copyfile(case.source / name, directory / name)
        verify_inputs(case, directory)
        with subprocess.Popen(
            command, cwd=directory, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        ) as process:
            for line in process.stdout:
                line = line.rstrip("\r\n")
                accepted = ACCEPTED_POINT.fullmatch(line)
                if accepted:
                    if on_progress is not None:
                        stage, bias = accepted.groups()
                        on_progress(int(stage) + 1, float(bias))
                elif line:
                    console_tail.append(line)
        if process.returncode == 0:
            return read_run_status(case, directory)
    except (OSError, ValueError, KeyError) as error:
        write_run_failure(log_path, invocation, error, console_tail)
        raise
    if not log_path.is_file() or not (directory / case.work / "run_manifest.json").is_file():
        write_run_failure(log_path, invocation, f"sd-sim exited with code {process.returncode}", console_tail)
    return ("Solver exited (code {})", process.returncode)


def read_run_status(case, directory):
    work = directory / case.work
    manifest = json.loads((work / "run_manifest.json").read_text(encoding="utf-8"))
    plan = json.loads((work / "stage_plan.json").read_text(encoding="utf-8"))
    stages = manifest["stages"]
    complete = (
        manifest["exit_code"] == 0 and stages and len(stages) == len(plan["stages"])
        and all(stage["runner_exit_code"] == 0 and stage["accepted_points"] > 0 for stage in stages)
        and all(stage["goal_reached"] for stage in stages[:-1])
    )
    if complete and case.stop_criterion:
        if stages[-1].get("stop_criterion_reached") == case.stop_criterion:
            return ("Current stop reached",)
    elif complete and stages[-1]["goal_reached"]:
        return ("Target reached",)
    return ("Some stage targets were not reached",)
