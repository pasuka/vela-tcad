"""Prepare frozen workbench cases and start the packaged solver."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from dataclasses import dataclass


@dataclass(frozen=True)
class Case:
    name: str
    source: Path
    arguments: list[str]
    hashes: dict[str, str]
    work: str
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
            case["hashes"] = {
                path.name: _hash(path) for path in sorted((catalog / key).iterdir())
                if path.name != case["original_script"]
            }
        (catalog / "cases.json").write_text(
            json.dumps(cases, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        catalog.rename(destination)


def read_cases(directory):
    manifest = json.loads((directory / "cases.json").read_text(encoding="utf-8"))
    return {key: Case(source=directory / key, **value) for key, value in manifest.items()}


def run_case(launcher, case, directory):
    """Use the frozen regression arguments, leaving published artifacts in place."""
    for name, checksum in case.hashes.items():
        copy = directory / name
        shutil.copyfile(case.source / name, copy)
        if _hash(copy) != checksum:
            raise ValueError(f"Input checksum mismatch: {name}")

    env = dict(os.environ, VELA_LINEAR_SOLVER=case.linear_solver, VELA_LINEAR_THREADS="1")
    if case.linear_solver == "pardiso_mkl":
        env.update(MKL_THREADING_LAYER="SEQUENTIAL", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
        env.pop("MKL_CBWR", None)
    command = [str(launcher), *case.arguments]
    with (directory / "sd-sim.log").open("wb") as log:
        log.write((subprocess.list2cmdline(command) + "\n\n").encode("utf-8"))
        log.flush()
        result = subprocess.run(
            command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        )
    if result.returncode:
        return ("Solver exited (code {})", result.returncode)
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
