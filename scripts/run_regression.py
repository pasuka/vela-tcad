#!/usr/bin/env python3
"""Run the Vela engineering examples as a lightweight regression suite."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

EXAMPLES = [
    {
        "name": "pn_diode",
        "config": Path("examples/pn_diode/simulation.json"),
        "expected": [Path("outputs/pn_iv.csv"), Path("outputs/pn_sweep_0000_0V.vtk")],
        "checks": ["csv_converged", "finite_outputs", "iv_trend"],
    },
    {
        "name": "moscap",
        "config": Path("examples/moscap/simulation.json"),
        "expected": [Path("outputs/moscap_poisson.vtk")],
        "checks": ["finite_outputs", "moscap_interface"],
    },
    {
        "name": "nmos2d",
        "config": Path("examples/nmos2d/simulation.json"),
        "expected": [Path("outputs/nmos_poisson.vtk")],
        "checks": ["finite_outputs"],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--runner", type=Path, default=None, help="Path to vela_example_runner")
    parser.add_argument("--workdir", type=Path, default=None, help="Directory for copied examples and outputs")
    return parser.parse_args()


def default_runner(repo: Path) -> Path:
    candidates = [repo / "build" / "vela_example_runner", repo / "build" / "src" / "tools" / "vela_example_runner"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def copy_example(repo: Path, workdir: Path, name: str) -> Path:
    src = repo / "examples" / name
    dst = workdir / name
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    (dst / "outputs").mkdir(parents=True, exist_ok=True)
    return dst


def has_nan_or_inf(path: Path) -> bool:
    text = path.read_text(errors="ignore")
    for token in text.replace(",", " ").split():
        try:
            value = float(token)
        except ValueError:
            lowered = token.lower()
            if lowered in {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}:
                return True
            continue
        if not math.isfinite(value):
            return True
    return False


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def read_vtk_scalar(path: Path, field: str) -> list[float]:
    lines = path.read_text().splitlines()
    marker = f"SCALARS {field} double 1"
    for idx, line in enumerate(lines):
        if line.strip() == marker:
            start = idx + 2
            values = []
            for raw in lines[start:]:
                stripped = raw.strip()
                if not stripped or stripped.startswith("SCALARS") or stripped.startswith("POINT_DATA"):
                    break
                try:
                    values.append(float(stripped))
                except ValueError:
                    break
            return values
    raise AssertionError(f"VTK field '{field}' not found in {path}")


def node_index(mesh: dict[str, Any], x: float, y: float) -> int:
    best = min(mesh["nodes"], key=lambda n: abs(n["x"] - x) + abs(n["y"] - y))
    if abs(best["x"] - x) > 1e-15 or abs(best["y"] - y) > 1e-15:
        raise AssertionError(f"No mesh node near ({x}, {y})")
    return int(best["id"])


def check_csv_converged(example_dir: Path) -> None:
    rows = read_csv(example_dir / "outputs" / "pn_iv.csv")
    if not rows:
        raise AssertionError("PN IV CSV contains no sweep rows")
    bad = [row for row in rows if row.get("converged") != "1"]
    if bad:
        raise AssertionError(f"Non-converged PN sweep rows: {bad}")


def check_iv_trend(example_dir: Path) -> None:
    rows = read_csv(example_dir / "outputs" / "pn_iv.csv")
    currents = [abs(float(row["total_current"])) for row in rows]
    if currents[-1] + 1e-40 < currents[0]:
        raise AssertionError(f"PN diode current trend is not forward-increasing: {currents}")


def check_moscap_interface(example_dir: Path) -> dict[str, float]:
    cfg = json.loads((example_dir / "simulation.json").read_text())
    reg = cfg["regression"]
    mesh = json.loads((example_dir / cfg["mesh_file"]).read_text())
    psi = read_vtk_scalar(example_dir / cfg["output_vtk"], "potential_V")
    x = float(reg["x_probe"])
    y_if = float(reg["interface_y"])
    y_si = float(reg["silicon_probe_y"])
    y_ox = float(reg["oxide_probe_y"])
    i_if_a = node_index(mesh, x, y_if)
    i_if_b = node_index(mesh, x, y_if)
    jump = abs(psi[i_if_a] - psi[i_if_b])
    if jump > float(reg["potential_jump_tol"]):
        raise AssertionError(f"MOSCAP interface potential jump {jump} exceeds tolerance")
    eps_si = 11.7
    eps_ox = 3.9
    field_si = abs((psi[i_if_a] - psi[node_index(mesh, x, y_si)]) / (y_if - y_si))
    field_ox = abs((psi[node_index(mesh, x, y_ox)] - psi[i_if_b]) / (y_ox - y_if))
    disp_si = eps_si * field_si
    disp_ox = eps_ox * field_ox
    denom = max(abs(disp_si), abs(disp_ox), 1e-300)
    rel = abs(disp_si - disp_ox) / denom
    if rel > float(reg["displacement_rel_tol"]):
        raise AssertionError(f"MOSCAP displacement mismatch {rel} exceeds tolerance")
    return {"potential_jump": jump, "displacement_rel_error": rel}


def run_example(runner: Path, repo: Path, workdir: Path, spec: dict[str, Any]) -> dict[str, Any]:
    name = spec["name"]
    example_dir = copy_example(repo, workdir, name)
    config = example_dir / "simulation.json"
    proc = subprocess.run([str(runner), "--config", str(config)], cwd=example_dir, text=True, capture_output=True)
    result: dict[str, Any] = {
        "name": name,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
        "passed": False,
        "checks": {},
    }
    try:
        if proc.returncode != 0:
            raise AssertionError(f"runner exited with {proc.returncode}: {proc.stderr.strip()}")
        for rel in spec["expected"]:
            out = example_dir / rel
            if not out.exists() or out.stat().st_size == 0:
                raise AssertionError(f"missing or empty output file: {rel}")
        if "finite_outputs" in spec["checks"]:
            for rel in spec["expected"]:
                if has_nan_or_inf(example_dir / rel):
                    raise AssertionError(f"NaN/Inf detected in {rel}")
            result["checks"]["finite_outputs"] = True
        if "csv_converged" in spec["checks"]:
            check_csv_converged(example_dir)
            result["checks"]["csv_converged"] = True
        if "iv_trend" in spec["checks"]:
            check_iv_trend(example_dir)
            result["checks"]["iv_trend"] = True
        if "moscap_interface" in spec["checks"]:
            result["checks"]["moscap_interface"] = check_moscap_interface(example_dir)
        result["passed"] = True
    except Exception as ex:  # noqa: BLE001 - regression summary should capture all failures.
        result["error"] = str(ex)
    return result


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    runner = (args.runner or default_runner(repo)).resolve()
    workdir = (args.workdir or (repo / "build" / "regression_output")).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    if not runner.exists():
        print(f"Regression runner not found: {runner}", file=sys.stderr)
        return 2
    results = [run_example(runner, repo, workdir, spec) for spec in EXAMPLES]
    summary = {"runner": str(runner), "workdir": str(workdir), "passed": all(r["passed"] for r in results), "examples": results}
    (workdir / "regression_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
