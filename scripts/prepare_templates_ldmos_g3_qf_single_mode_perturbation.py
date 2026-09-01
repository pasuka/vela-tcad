#!/usr/bin/env python3
"""Prepare a symmetric one-node electron-QF perturbation for LDMOS G3."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--carrier-config", type=Path, required=True)
    parser.add_argument("--sentaurus-deck", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--node", type=int, default=3721)
    parser.add_argument("--amplitude-v", type=float, default=1.0e-6)
    args = parser.parse_args()

    if args.amplitude_v <= 0.0:
        raise ValueError("amplitude must be positive")
    state_path = args.state.resolve()
    config_path = args.carrier_config.resolve()
    output_dir = args.output_dir.resolve()
    fields, base_rows = read_csv(state_path)
    if "node_id" not in fields or "phin" not in fields:
        raise ValueError("state CSV must contain node_id and phin")

    matching = [row for row in base_rows if int(row["node_id"]) == args.node]
    if len(matching) != 1:
        raise ValueError(f"expected exactly one row for node {args.node}")
    base_phin = float(matching[0]["phin"])
    base_config = json.loads(config_path.read_text(encoding="utf-8"))
    if base_config.get("simulation_type") != "newton_carrier_term_probe":
        raise ValueError("carrier config must use newton_carrier_term_probe")
    deck_template = None
    if args.sentaurus_deck is not None:
        deck_template = args.sentaurus_deck.resolve().read_text(encoding="utf-8")
        for required in (
            'Load(FilePrefix="state_high_vsv")',
            'Output="probe_high_vsv_qf.log"',
            'NewtonPlot="probe_high_vsv_qf_newton_%d_%d.tdr"',
            'Plot(FilePrefix="probe_high_vsv_qf_loaded" NoOverWrite)',
        ):
            if required not in deck_template:
                raise ValueError(f"Sentaurus deck is missing expected token: {required}")

    variants: dict[str, dict[str, Any]] = {}
    for label, sign in (("minus", -1.0), ("plus", 1.0)):
        rows = [dict(row) for row in base_rows]
        perturbed = base_phin + sign * args.amplitude_v
        for row in rows:
            if int(row["node_id"]) == args.node:
                row["phin"] = format(perturbed, ".17g")
                break
        state_out = output_dir / f"state_node{args.node}_{label}.csv"
        carrier_out = output_dir / f"carrier_terms_node{args.node}_{label}.csv"
        config_out = output_dir / f"carrier_probe_node{args.node}_{label}.json"
        write_csv(state_out, fields, rows)

        config = dict(base_config)
        config["_comment"] = (
            f"Controlled electron-QF perturbation at node {args.node}: "
            f"{sign:+.0f} * {args.amplitude_v:.17g} V; all other state fields frozen."
        )
        config["state_file"] = str(state_out)
        config["output_csv"] = str(carrier_out)
        write_json(config_out, config)
        deck_out = None
        if deck_template is not None:
            prefix = f"probe_high_vsv_node{args.node}_{label}"
            deck = deck_template.replace(
                'Load(FilePrefix="state_high_vsv")',
                f'Load(FilePrefix="state_node{args.node}_{label}")',
            )
            deck = deck.replace("probe_high_vsv_qf", prefix)
            deck_out = output_dir / f"{prefix}.cmd"
            deck_out.write_text(deck, encoding="utf-8")
        variants[label] = {
            "sign": sign,
            "state_file": str(state_out),
            "carrier_config": str(config_out),
            "carrier_output": str(carrier_out),
            "sentaurus_deck": str(deck_out) if deck_out is not None else None,
            "target_phin_V": perturbed,
        }

    manifest = {
        "contract": {
            "state": "identical_high_endpoint_vsv_except_single_electron_qf_dof",
            "node_id": args.node,
            "amplitude_V": args.amplitude_v,
            "difference": "symmetric_central",
            "density_columns": "retained; both solvers reconstruct carrier density from QF unknowns",
            "production_defaults_changed": False,
        },
        "base_state": str(state_path),
        "base_carrier_config": str(config_path),
        "base_phin_V": base_phin,
        "variants": variants,
    }
    write_json(output_dir / "manifest.json", manifest)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
