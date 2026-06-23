#!/usr/bin/env python3
"""Prepare reproducible PN2D BV localization decks from an imported base deck."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--stop", required=True, type=float)
    parser.add_argument("--qf-limit", type=float)
    parser.add_argument("--max-update", type=float)
    parser.add_argument("--impact-model", choices=["keep", "none"], default="keep")
    parser.add_argument("--diagnostics", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_config = args.base_config.resolve()
    reference_root = args.reference_root.resolve()
    out_root = args.out_root.resolve()

    base = json.loads(base_config.read_text(encoding="utf-8-sig"))
    out_dir = out_root / args.case_name
    out_dir.mkdir(parents=True, exist_ok=True)

    base["mesh_file"] = str(reference_root / "vela" / "mesh.json")
    base["node_doping_file"] = str(reference_root / "doping.csv")
    base["materials_file"] = str(
        reference_root / "vela" / "pn2d_sentaurus2018_iv_materials.json"
    )
    base["output_csv"] = str(out_dir / "iv.csv")

    solver = base.setdefault("solver", {})
    if args.max_update is not None:
        solver["max_update"] = args.max_update
    if args.qf_limit is not None:
        solver["quasi_fermi_update_limit_V"] = args.qf_limit
    if args.impact_model == "none":
        solver["impact_ionization"] = {"model": "none"}
    if args.diagnostics:
        solver["diagnostics"] = True

    sweep = base.setdefault("sweep", {})
    sweep["stop"] = args.stop
    sweep["step"] = -0.05
    sweep["write_vtk"] = False
    sweep["max_step"] = 0.05
    sweep["min_step"] = 1.0e-10
    sweep["max_retries"] = 29
    if args.diagnostics:
        sweep["diagnostics"] = {
            "newton_history": {
                "enabled": True,
                "csv_file": str(out_dir / "newton_history.csv"),
            }
        }

    config_path = out_dir / "config.json"
    config_path.write_text(json.dumps(base, indent=2), encoding="utf-8")
    print(config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
