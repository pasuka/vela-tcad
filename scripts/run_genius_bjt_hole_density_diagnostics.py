#!/usr/bin/env python3
"""Run non-baseline Genius NPN BJT hole-density solver diagnostics."""

from __future__ import annotations

import argparse
import copy
import json
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_OUTPUT = BUILD_ROOT / "hole_density_diagnosis"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
TRACE_NODES = [589, 144, 150, 7, 12]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase", choices=("probes", "ab", "enforce", "all"), default="all"
    )
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--fixture-root", type=Path, default=FIXTURE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def absolute(path: Path) -> str:
    return str(path.resolve())


def load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def run_config(
    runner: Path,
    config: Path,
    output_dir: Path,
    stage: str,
    require_success: bool = True,
) -> dict[str, object]:
    process = subprocess.run(
        [absolute(runner), "--config", absolute(config)],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    (output_dir / f"{stage}.stdout.log").write_text(
        process.stdout, encoding="utf-8", newline="\n"
    )
    (output_dir / f"{stage}.stderr.log").write_text(
        process.stderr, encoding="utf-8", newline="\n"
    )
    status: dict[str, object] = {
        "stage": stage,
        "config": absolute(config),
        "return_code": process.returncode,
    }
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if lines:
        try:
            status["runner_status"] = json.loads(lines[-1])
        except json.JSONDecodeError:
            status["runner_status"] = None
    write_json(output_dir / f"{stage}.status.json", status)
    if require_success and process.returncode != 0:
        raise RuntimeError(
            f"{stage} failed with exit code {process.returncode}; "
            f"see {output_dir / f'{stage}.stderr.log'}"
        )
    return status


def make_common_absolute(config: dict[str, object], fixture: Path) -> None:
    config["mesh_file"] = absolute(fixture / "vela" / "input" / "mesh.json")
    config["node_doping_file"] = absolute(fixture / "vela" / "input" / "doping.csv")
    config["materials_file"] = absolute(fixture / "vela" / "materials_sentaurus2022.json")


def run_probes(args: argparse.Namespace) -> list[dict[str, object]]:
    output = args.output_root / "probes"
    output.mkdir(parents=True, exist_ok=True)
    base = load(args.fixture_root / "vela" / "configs" / "m1_spatial_vce3.json")
    make_common_absolute(base, args.fixture_root)
    base.pop("sweep", None)
    base["state_file"] = absolute(BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv")
    statuses = []
    for simulation_type, name in [
        ("newton_carrier_term_probe", "carrier_terms"),
        ("newton_carrier_row_probe", "carrier_rows"),
    ]:
        config = copy.deepcopy(base)
        config["simulation_type"] = simulation_type
        config["output_csv"] = absolute(output / f"{name}.csv")
        if simulation_type == "newton_carrier_term_probe":
            config["carrier_term_probe"] = {"solved_equation_terms": True}
        path = output / f"{name}.json"
        write_json(path, config)
        statuses.append(run_config(args.runner, path, output, name))
    return statuses


def add_solver_diagnostics(
    config: dict[str, object], output: Path, stage: str
) -> None:
    solver = config["solver"]
    carrier = solver.setdefault("carrier_row_convergence", {})
    carrier["diagnostic_csv"] = absolute(output / f"{stage}_carrier_violations.csv")
    carrier["trace_csv"] = absolute(output / f"{stage}_carrier_trace.csv")
    carrier["trace_nodes"] = TRACE_NODES
    carrier["trace_first_iterations"] = 100
    carrier["trace_every_iterations"] = 1
    solver["local_update_diagnostics"] = {
        "enabled": True,
        "csv_file": absolute(output / f"{stage}_local_updates.csv"),
        "nodes": TRACE_NODES,
        "first_iterations": 100,
        "every_iterations": 1,
    }
    for name in (
        f"{stage}_carrier_violations.csv",
        f"{stage}_carrier_trace.csv",
        f"{stage}_local_updates.csv",
    ):
        path = output / name
        if path.exists():
            path.unlink()


def run_cap_step_probe(args: argparse.Namespace, cap_V: float) -> dict[str, object]:
    name = "cap_disabled" if cap_V == 0.0 else f"cap_{cap_V:g}V"
    output = args.output_root / "ab_step" / name
    output.mkdir(parents=True, exist_ok=True)
    config = load(args.fixture_root / "vela" / "configs" / "m1_spatial_vce3.json")
    make_common_absolute(config, args.fixture_root)
    config.pop("sweep", None)
    config["simulation_type"] = "newton_step_probe"
    config["state_file"] = absolute(
        BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv"
    )
    config["output_csv"] = absolute(output / "one_step.csv")
    config["solver"]["quasi_fermi_update_limit_V"] = cap_V
    path = output / "one_step_config.json"
    write_json(path, config)
    status = run_config(args.runner, path, output, "one_step")
    result = {
        "name": name,
        "quasi_fermi_update_limit_V": cap_V,
        "starting_state": config["state_file"],
        "status": status,
    }
    write_json(output / "variant_manifest.json", result)
    return result


def run_enforced_state(args: argparse.Namespace) -> dict[str, object]:
    output = args.output_root / "enforce"
    output.mkdir(parents=True, exist_ok=True)
    config = load(args.fixture_root / "vela" / "configs" / "m1_spatial_vce3.json")
    make_common_absolute(config, args.fixture_root)
    config.pop("sweep", None)
    config["simulation_type"] = "newton_solve_from_state"
    config["state_file"] = absolute(
        BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv"
    )
    config["output_state_file"] = absolute(output / "enforced_state.csv")
    carrier = config["solver"]["carrier_row_convergence"]
    carrier["mode"] = "enforce"
    carrier["eps_row"] = 1.0e-4
    carrier["min_newton_max_iter"] = 80
    add_solver_diagnostics(config, output, "enforced")
    path = output / "enforced_config.json"
    write_json(path, config)
    status = run_config(
        args.runner, path, output, "enforced", require_success=False
    )
    result = {
        "name": "cap_0.1V_carrier_rows_enforced",
        "quasi_fermi_update_limit_V": 0.1,
        "carrier_row_mode": "enforce",
        "carrier_row_eps": 1.0e-4,
        "starting_state": config["state_file"],
        "status": status,
    }
    write_json(output / "enforce_manifest.json", result)
    return result


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "baseline_unchanged": True,
        "trace_nodes": TRACE_NODES,
        "runs": [],
    }
    if args.phase in ("probes", "all"):
        manifest["runs"].append({"phase": "probes", "statuses": run_probes(args)})
    if args.phase in ("ab", "all"):
        for cap in (0.1, 0.05, 0.025, 0.0):
            manifest["runs"].append(run_cap_step_probe(args, cap))
    if args.phase in ("enforce", "all"):
        manifest["runs"].append(run_enforced_state(args))
    write_json(args.output_root / "run_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
