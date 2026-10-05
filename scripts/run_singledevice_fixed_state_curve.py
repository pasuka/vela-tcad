#!/usr/bin/env python3
"""Replay a 21-point Sentaurus SingleDevice state series through Vela current."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]


def run(command: list[str], cwd: Path | None = None) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def bias_for_index(index: int, intervals: int = 20) -> float:
    return -0.5 + (2.2 - (-0.5)) * index / intervals


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states-dir", type=Path, required=True)
    parser.add_argument("--vela-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--tdr-importer", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    args = parser.parse_args()

    args.work_dir.mkdir(parents=True, exist_ok=True)
    (args.work_dir / "exports").mkdir(exist_ok=True)
    (args.work_dir / "inventories").mkdir(exist_ok=True)
    for name in ("mesh.json", "doping.csv"):
        shutil.copy2(args.vela_dir / name, args.work_dir / name)

    summary_rows: list[dict[str, str | float]] = []
    for branch in ("lin", "sat"):
        base = json.loads(
            (args.vela_dir / f"simulation_idvg_{branch}.json").read_text(encoding="utf-8-sig"))
        material = Path(base["materials_file"])
        material_source = material if material.is_absolute() else args.vela_dir / material
        shutil.copy2(material_source, args.work_dir / material.name)
        base["materials_file"] = material.name
        for index in range(21):
            stem = f"{branch}_{index:04d}"
            tdr = args.states_dir / f"{branch}_state_{index:04d}_des.tdr"
            export = args.work_dir / "exports" / stem
            run([
                str(args.tdr_importer), "--tdr", str(tdr),
                "--inventory-json", str(args.work_dir / "inventories" / f"{stem}.json"),
                "--export-dir", str(export),
                "--compensated-doping-policy", "reported",
            ])
            state = args.work_dir / f"{stem}_state.csv"
            run([
                sys.executable, str(REPO / "scripts" / "sentaurus_fields_to_restart.py"),
                "--export-dir", str(export), "--output", str(state),
            ])
            bias = bias_for_index(index)
            deck = json.loads(json.dumps(base))
            deck["solver"]["method"] = "frozen_state"
            deck["sweep"].update({
                "start": bias,
                "stop": bias,
                "step": 0.135,
                "initial_state_file": state.name,
                "frozen_state_compute_current": True,
            })
            deck["output_csv"] = f"{stem}_current.csv"
            deck_path = args.work_dir / f"{stem}.json"
            deck_path.write_text(json.dumps(deck, indent=2) + "\n")
            run([str(args.runner), "--config", deck_path.name, "--log", "off"], cwd=args.work_dir)
            with (args.work_dir / deck["output_csv"]).open(newline="") as handle:
                row = next(csv.DictReader(handle))
            summary_rows.append({
                "branch": branch,
                "state_index": index,
                "bias_V": bias,
                "current_total_A_per_um": row["current_total_A_per_um"],
                "current_electron_A_per_um": row["current_electron_A_per_um"],
                "current_hole_A_per_um": row["current_hole_A_per_um"],
            })

    with (args.work_dir / "singledevice_fixed_state_curves.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
