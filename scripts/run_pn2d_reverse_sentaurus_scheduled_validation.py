#!/usr/bin/env python3
"""Run PN2D BV validation with Sentaurus-aligned adaptive continuation."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

from generate_pn2d_config import render_named_template

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--base-config", type=Path)
    source.add_argument(
        "--template",
        action="store_true",
        help="Use the versioned pn2d_bv production template.",
    )
    parser.add_argument("--mesh-file", type=Path)
    parser.add_argument("--node-doping-file", type=Path)
    parser.add_argument("--materials-file", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    runner = args.runner.resolve()
    if args.template:
        missing = [
            name for name, value in (
                ("--mesh-file", args.mesh_file),
                ("--node-doping-file", args.node_doping_file),
                ("--materials-file", args.materials_file),
            ) if value is None
        ]
        if missing:
            parser.error(f"--template requires {', '.join(missing)}")
        base, _ = render_named_template(
            "pn2d_bv",
            {
                "solver_profile": "legacy",
                "mesh_file": str(args.mesh_file.resolve()),
                "node_doping_file": str(args.node_doping_file.resolve()),
                "materials_file": str(args.materials_file.resolve()),
            },
            allow_absolute_paths=True,
        )
        base_config_source = "configs/templates/pn2d_bv.template.json"
    else:
        base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
        base_config_source = str(args.base_config.resolve())
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, object] = {}

    for basis in ("net_doping", "cell_reconstructed_total_impurity"):
        for impact in ("on", "off"):
            tag = f"{basis}_{impact}"
            case_dir = out_dir / tag
            case_dir.mkdir(parents=True, exist_ok=True)
            cfg = json.loads(json.dumps(base))
            cfg["solver"]["diagnostics"] = False
            cfg["solver"]["max_iter"] = 80
            cfg["solver"]["handoff"]["newton_max_iter"] = 80
            cfg["solver"]["mobility"]["doping_concentration_basis"] = basis
            if impact == "off":
                cfg["solver"]["impact_ionization"] = {"model": "none"}

            sweep = cfg["sweep"]
            sweep.pop("bias_points", None)
            sweep.update(
                {
                    "start": 0.0,
                    "stop": -20.0,
                    "step": -0.05,
                    "initial_step": 1.0e-4,
                    "min_step": 1.0e-10,
                    "max_step": 0.05,
                    "growth_factor": 1.2,
                    "shrink_factor": 0.5,
                    "max_retries": 20,
                    "initialization": {"mode": "poisson_block"},
                    "write_vtk": False,
                    "diagnostics": {
                        "newton_history": {
                            "enabled": True,
                            "csv_file": str((case_dir / f"{tag}_newton.csv").resolve()),
                        }
                    },
                }
            )
            output = case_dir / f"{tag}.csv"
            cfg["output_csv"] = str(output.resolve())
            config_path = case_dir / f"simulation_{tag}.json"
            config_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
            log_path = case_dir / "runner.log"
            with log_path.open("w", encoding="utf-8") as log:
                result = subprocess.run(
                    [str(runner), "--config", str(config_path)],
                    cwd=case_dir,
                    env=os.environ.copy(),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            artifacts[tag] = {
                "returncode": result.returncode,
                "config": str(config_path),
                "output": str(output),
                "log": str(log_path),
            }
            if result.returncode:
                raise RuntimeError(f"{tag} failed; inspect {log_path}")

    (out_dir / "validation_manifest.json").write_text(
        json.dumps(
            {
                "base_config": base_config_source,
                "sentaurus_alignment": {
                    "initial_step": 1.0e-4,
                    "min_step": 1.0e-10,
                    "max_step": 0.05,
                    "growth_factor": 1.2,
                    "shrink_factor": 0.5,
                    "max_retries": 20,
                    "newton_max_iter": 80,
                    "initialization": "poisson_block",
                },
                "artifacts": artifacts,
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
