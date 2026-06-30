#!/usr/bin/env python3
"""Characterize the pn2d BV reverse-bias continuation/avalanche robustness gap.

The generated pn2d Sentaurus2018 BV deck is a focused low-bias smoke gate. The
avalanche physics and its analytic Jacobian already exist and are unit-tested;
what is unknown is the *empirical* continuation reach on the real pn2d mesh.

This probe is read-only with respect to the committed deck and the solver. It
copies the generated BV deck and runs progressively harder variants, each with
adaptive step control (retries + shrink/grow), and reports, per variant:
- the deepest reverse bias reached before non-convergence,
- the number of accepted sweep points,
- whether the full requested stop was reached.

Variants (all extend the sweep to a deeper reverse stop with adaptive steps):
- ``cont_base``     : current models (SRH, no avalanche) — continuation-only.
- ``cont_srh``      : same as base (kept for symmetry / explicit label).
- ``cont_avalanche``: add Selberherr impact ionization (the real BV target).

Nothing here changes ``reference_tcad/pn2d_sentaurus2018`` or any solver source;
outputs land under ``build/pn2d_bv_continuation_probe``.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        type=Path,
        default=REPO / "build" / "reference_tcad" / "pn2d_sentaurus2018",
        help="Generated reference workspace containing vela/simulation_bv.json.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "build" / "pn2d_bv_continuation_probe",
        help="Output directory for probe decks and curves.",
    )
    parser.add_argument(
        "--runner",
        type=Path,
        default=REPO / "build" / "vela_example_runner.exe",
        help="Path to vela_example_runner executable.",
    )
    parser.add_argument(
        "--stop",
        type=float,
        default=-5.0,
        help="Reverse-bias Anode voltage goal to attempt.",
    )
    parser.add_argument(
        "--step",
        type=float,
        default=-0.25,
        help="Nominal reverse-bias step.",
    )
    return parser.parse_args()


def resolve_input_paths(deck: dict[str, Any], base_dir: Path) -> dict[str, Any]:
    resolved = json.loads(json.dumps(deck))
    for key in [
        "mesh_file",
        "node_doping_file",
        "node_mobility_file",
        "node_lifetime_file",
        "interface_charge_file",
    ]:
        raw = resolved.get(key)
        if not isinstance(raw, str):
            continue
        candidate = Path(raw)
        if candidate.is_absolute():
            continue
        from_base = (base_dir / candidate).resolve()
        if from_base.exists():
            resolved[key] = str(from_base)
    return resolved


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def read_curve(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def build_variant(
    base_deck: dict[str, Any],
    base_dir: Path,
    tag: str,
    stop: float,
    step: float,
    impact_model: str,
) -> dict[str, Any]:
    deck = resolve_input_paths(json.loads(json.dumps(base_deck)), base_dir)
    solver = deck.setdefault("solver", {})
    solver["impact_ionization"] = {"model": impact_model}

    sweep = deck.setdefault("sweep", {})
    sweep["stop"] = stop
    sweep["step"] = step
    sweep["max_retries"] = 8
    sweep["shrink_factor"] = 0.5
    sweep["growth_factor"] = 1.2
    sweep["stop_on_failure"] = True
    deck["output_csv"] = f"pn2d_bv_{tag}.csv"
    return deck


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    runner = args.runner.resolve()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    base_bv_deck = json.loads((workspace / "vela" / "simulation_bv.json").read_text())
    base_deck_dir = workspace / "vela"

    variants = [
        ("cont_base", "none"),
        ("cont_avalanche", "selberherr"),
    ]

    rows: list[dict[str, Any]] = []
    for tag, impact_model in variants:
        deck = build_variant(
            base_bv_deck, base_deck_dir, tag, args.stop, args.step, impact_model
        )
        cfg = out_dir / f"simulation_bv_{tag}.json"
        write_json(cfg, deck)
        proc = subprocess.run(
            [str(runner), "--config", str(cfg)],
            cwd=out_dir,
            capture_output=True,
            text=True,
        )
        curve = read_curve(out_dir / deck["output_csv"])
        # bv_reverse mode records bias_V as the positive reverse-bias magnitude,
        # so the continuation depth is the largest-magnitude accepted bias.
        biases = [abs(float(r["bias_V"])) for r in curve] if curve else []
        deepest = max(biases) if biases else 0.0
        reached_stop = bool(biases) and deepest >= abs(args.stop) - 1.0e-9
        failure_reason = ""
        for r in curve:
            if r.get("converged") in ("0", "false", "False") and r.get("failure_reason"):
                failure_reason = r["failure_reason"]
                break
        rows.append(
            {
                "tag": tag,
                "impact_model": impact_model,
                "requested_stop_V": -abs(args.stop),
                "accepted_points": len(curve),
                "deepest_reverse_bias_V": deepest,
                "reached_requested_stop": reached_stop,
                "failure_reason": failure_reason,
                "runner_exit_code": proc.returncode,
            }
        )

    summary_csv = out_dir / "pn2d_bv_continuation_summary.csv"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("\n=== pn2d BV reverse-bias continuation probe ===")
    for row in rows:
        print(
            f"{row['tag']} (impact={row['impact_model']}): "
            f"points={row['accepted_points']}, "
            f"deepest_reverse={row['deepest_reverse_bias_V']:.3f} V "
            f"(requested {abs(row['requested_stop_V']):.3f} V), "
            f"reached_stop={row['reached_requested_stop']}, "
            f"failure={row['failure_reason'] or 'none'}, "
            f"exit={row['runner_exit_code']}"
        )
    print(f"\nSummary CSV: {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
