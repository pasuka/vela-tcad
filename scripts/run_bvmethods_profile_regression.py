"""Run the autonomous opt-in BVMethods operator on physical QF checkpoints.

Native SG caches, source arrays, matrices and constitutive corrections are
not valid inputs. A separate comparison may inspect the exported operators.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import subprocess


def validate_manifest(data):
    if set(data) != {"config", "frames", "outer_targets_V"}:
        raise ValueError("manifest accepts only config, frames and outer_targets_V")
    frames = data["frames"]
    if not frames or frames[0].get("label") != "base":
        raise ValueError("first checkpoint must be base")
    labels = set()
    count = None
    for f in frames:
        if set(f) != {"label", "inner_V", "outer_V", "psi_V", "fn_V", "fp_V"}:
            raise ValueError("checkpoint accepts only potentials and contact voltages")
        if f["label"] in labels:
            raise ValueError("duplicate checkpoint label")
        labels.add(f["label"])
        n = len(f["psi_V"])
        if n <= 2587 or len(f["fn_V"]) != n or len(f["fp_V"]) != n:
            raise ValueError("BVMethods checkpoints require the same full mesh potentials")
        if count is not None and n != count:
            raise ValueError("checkpoint mesh dimensions differ")
        count = n
        numbers = f["psi_V"] + f["fn_V"] + f["fp_V"] + [f["inner_V"], f["outer_V"]]
        if not all(math.isfinite(x) for x in numbers):
            raise ValueError("nonfinite checkpoint")
    targets = data["outer_targets_V"]
    if len(targets) < 3 or targets[0] != 0 or not all(math.isfinite(v) for v in targets):
        raise ValueError("require at least three finite targets beginning at zero")
    if any(b <= a for a, b in zip(targets, targets[1:])):
        raise ValueError("targets must increase strictly")


def validate_summary(data, manifest):
    if data.get("operator") != "bvmethods_tri3_qf":
        raise ValueError("wrong operator")
    frames = data.get("frames", [])
    if [f["label"] for f in frames] != [f["label"] for f in manifest["frames"]]:
        raise ValueError("missing checkpoint exports")
    expected = 3 * len(manifest["frames"][0]["psi_V"]) + 1
    if any(f["dimension"] != expected or f["ledger_residual_difference_Linf"] != 0 for f in frames):
        raise ValueError("operator dimension or ledger mismatch")
    points = data.get("sweep", [])
    if [p["outer_V"] for p in points] != manifest["outer_targets_V"]:
        raise ValueError("incomplete multistep sweep")
    for p in points:
        if not p["converged"] or abs(p["load_V"]) > 1e-6:
            raise ValueError("coupled load-line did not converge")
        if not all(math.isfinite(p[k]) for k in ["inner_V", "current_A_per_um", "equation_inf", "load_V"]):
            raise ValueError("nonfinite accepted point")
    accepted = sum(t["alpha"] > 0 for p in data.get("continuation", points) for t in p["trace"])
    if accepted < 5:
        raise ValueError("regression did not exercise at least five accepted Newton updates")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text(encoding="utf-8"))
    validate_manifest(data)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "console.log").open("w", encoding="utf-8") as stream:
        result = subprocess.run([str(args.probe.resolve()), str(args.manifest.resolve()), str(args.output.resolve())],
                                stdout=stream, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        raise SystemExit(f"probe failed with exit {result.returncode}; see console.log")
    summary = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))
    validate_summary(summary, data)
    print(f"{len(summary['frames'])} checkpoints and {len(summary['sweep'])} outer targets passed")


if __name__ == "__main__":
    main()
