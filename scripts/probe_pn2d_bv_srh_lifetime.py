#!/usr/bin/env python3
"""P1 research probe: BV SRH constant-lifetime parity sweep.

The Sentaurus pn2d BV deck declares plain ``Recombination(SRH ...)`` with no
``DopingDependence`` keyword and logs ``no Lifetime file`` /
``no ModelParameters file``. In Sentaurus Device that means SRH uses a
*constant* lifetime equal to the silicon default ``taumax`` from the shipped
``Silicon.par`` parameter set (documented default ``1.0e-5 s`` for both
carriers), not a doping-dependent Scharfetter relation.

This probe quantifies how the BV 0.05 V parity orders move as the Vela BV SRH
lifetime is swept across:
- ``1.0e-6 s`` — the value currently used by the reference BV deck,
- ``3.0e-6 s`` — the prior tau-grid scan best for BV,
- ``1.0e-5 s`` — the documented Sentaurus silicon default ``taumax``.

It only overrides the SRH lifetime of an already-generated BV deck; nothing in
the solver or the committed reference deck is changed. Read-only research.
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

TAU_CANDIDATES_S = [
    (1.0e-6, "current reference BV deck"),
    (3.0e-6, "prior tau-grid scan best for BV"),
    (1.0e-5, "documented Sentaurus silicon default taumax"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace",
        type=Path,
        default=REPO / "build" / "pn2d_iv_promote_check",
        help="Generated reference workspace containing vela/simulation_bv.json.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "build" / "pn2d_bv_srh_probe",
        help="Output directory for probe decks, curves, and reports.",
    )
    parser.add_argument(
        "--runner",
        type=Path,
        default=REPO / "build" / "vela_example_runner.exe",
        help="Path to vela_example_runner executable.",
    )
    return parser.parse_args()


def read_curve(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


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


def main() -> int:
    args = parse_args()
    workspace = args.workspace.resolve()
    runner = args.runner.resolve()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    base_bv_deck = json.loads((workspace / "vela" / "simulation_bv.json").read_text())
    base_deck_dir = workspace / "vela"
    ref_bv = workspace / "reference_curves" / "pn2d_bv_reference.csv"
    if not ref_bv.exists():
        raise FileNotFoundError(f"Missing BV reference curve: {ref_bv}")

    rows: list[dict[str, Any]] = []
    for tau_s, note in TAU_CANDIDATES_S:
        tag = f"tau_{tau_s:.0e}".replace("-", "m").replace("+", "p")
        deck = resolve_input_paths(json.loads(json.dumps(base_bv_deck)), base_deck_dir)
        solver = deck.setdefault("solver", {})
        solver["taun"] = tau_s
        solver["taup"] = tau_s
        deck["output_csv"] = f"pn2d_bv_{tag}.csv"

        cfg = out_dir / f"simulation_bv_{tag}.json"
        write_json(cfg, deck)
        subprocess.run([str(runner), "--config", str(cfg)], cwd=out_dir, check=True)

        cmp_json = out_dir / f"pn2d_bv_0p05_{tag}.json"
        cmp_md = out_dir / f"pn2d_bv_0p05_{tag}.md"
        subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "compare_reference_curves.py"),
                "--reference",
                str(ref_bv),
                "--candidate",
                str(out_dir / deck["output_csv"]),
                "--output-json",
                str(cmp_json),
                "--output-md",
                str(cmp_md),
                "--kind",
                "iv",
                "--candidate-scale",
                "1",
                "--candidate-column",
                "current_total_A_per_um",
                "--bias-min",
                "0.05",
                "--bias-max",
                "0.05",
                "--min-points",
                "1",
            ],
            cwd=REPO,
            check=True,
        )

        bv_rows = read_curve(out_dir / deck["output_csv"])
        i_0p05 = None
        for row in bv_rows:
            if abs(float(row["bias_V"]) - 0.05) <= 1.0e-9:
                i_0p05 = float(row["current_total_A_per_um"])
                break

        cmp = json.loads(cmp_json.read_text())
        rows.append(
            {
                "tau_s": tau_s,
                "note": note,
                "bv_orders_at_0p05": cmp["iv"].get("orders_of_magnitude"),
                "i_0p05_A_per_um": i_0p05,
            }
        )

    summary_csv = out_dir / "pn2d_bv_srh_summary.csv"
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print("\n=== pn2d BV SRH constant-lifetime parity probe ===")
    for row in rows:
        print(
            f"tau={row['tau_s']:.1e} s ({row['note']}): "
            f"bv_orders_0p05={row['bv_orders_at_0p05']}, "
            f"I(0.05V)={row['i_0p05_A_per_um']:.3e} A/um"
        )
    print(f"\nSummary CSV: {summary_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
