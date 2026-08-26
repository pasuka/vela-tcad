#!/usr/bin/env python3
"""Create output-only derivative decks for Templates/LDMOS state capture.

The official materialized decks remain unchanged.  These derived decks add
only ``Plot`` statements at exact sweep coordinates and carry a manifest that
hashes both parent and derivative.  They must never be labelled as the
unaltered Applications Library oracle.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def threshold_neighborhood(curve: Path) -> float:
    with curve.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    points = [
        (float(row["bias_V"]), abs(float(row["current_total_A_per_um"])))
        for row in rows
    ]
    points = [(bias, current) for bias, current in points if math.isfinite(bias) and math.isfinite(current)]
    if len(points) < 3:
        raise ValueError(f"not enough finite Id-Vg points in {curve}")
    candidates: list[tuple[float, float]] = []
    for (v0, i0), (v1, i1) in zip(points, points[1:]):
        if v1 <= v0:
            continue
        slope = (math.log10(max(i1, 1.0e-300)) - math.log10(max(i0, 1.0e-300))) / (v1 - v0)
        candidates.append((slope, v1))
    if not candidates:
        raise ValueError(f"Id-Vg bias is not increasing in {curve}")
    return max(candidates)[1]


def add_after_once(text: str, marker: str, addition: str, label: str) -> str:
    count = text.count(marker)
    if count != 1:
        raise ValueError(f"{label}: expected one insertion marker, found {count}")
    return text.replace(marker, marker + addition, 1)


def build_idvg(text: str, threshold_v: float) -> tuple[str, dict[str, Any]]:
    text = add_after_once(
        text,
        "\tCoupled { Poisson Electron Hole }\n",
        '\tPlot(-Loadable FilePrefix="state_idvg_eq_0V")\n',
        "IdVg equilibrium",
    )
    times = sorted({0.0, threshold_v / 5.0, 0.5, 1.0})
    time_text = "; ".join(f"{value:.17g}" for value in times)
    marker = "\t\tCurrentPlot( Time= (Range=(0 1) Intervals= 30)  )\n"
    text = add_after_once(
        text,
        marker,
        f'\t\tPlot(-Loadable FilePrefix="state_idvg" NoOverWrite Time=({time_text}))\n',
        "IdVg sweep",
    )
    return text, {
        "threshold_neighborhood_V": threshold_v,
        "biases_V": [5.0 * value for value in times],
        "normalized_times": times,
    }


def build_idvd(text: str) -> tuple[str, dict[str, Any]]:
    biases = [0.0, 0.1, 1.0, 10.0, 20.0, 40.0]
    times = [value / 40.0 for value in biases]
    time_text = "; ".join(f"{value:.17g}" for value in times)
    markers = [
        "\t\tCurrentPlot( Time= (Range= (0 1) Intervals= 30) )\n",
        "\t\tCurrentPlot( Time= (Range= (0 1) Intervals= 30)  )\n",
    ]
    for gate, marker in zip((4, 8), markers, strict=True):
        text = add_after_once(
            text,
            marker,
            f'\t\tPlot(-Loadable FilePrefix="state_idvd_vg{gate}" NoOverWrite Time=({time_text}))\n',
            f"IdVd Vg={gate}",
        )
    return text, {"gate_biases_V": [4.0, 8.0], "drain_biases_V": biases, "normalized_times": times}


def build_bv(text: str) -> tuple[str, dict[str, Any]]:
    marker = '\t){ Coupled { Poisson Electron Hole Temperature } }\n'
    addition = (
        '\t\tPlot(-Loadable FilePrefix="state_bv_path" NoOverWrite '
        'Time=(Range=(0 1) Intervals=30))\n'
    )
    if text.count(marker) != 1:
        raise ValueError(f"BV continuation: expected one insertion marker, found {text.count(marker)}")
    replacement = (
        '\t){ Coupled { Poisson Electron Hole Temperature }\n'
        + addition +
        '\t}\n'
    )
    return text.replace(marker, replacement, 1), {
        "path_samples": 31,
        "selection_policy": "select exact saved continuation states after matching state terminal current/voltage to the original PLT",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--idvg-curve", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to overwrite derivative decks: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    threshold_v = threshold_neighborhood(args.idvg_curve)
    builders = {
        "IdVg.cmd": lambda value: build_idvg(value, threshold_v),
        "IdVd.cmd": build_idvd,
        "BVdss.cmd": build_bv,
    }
    manifest: dict[str, Any] = {
        "schema": "vela.templates_ldmos.state_deck_materialization.v1",
        "classification": "derived_output_only_control_not_official_oracle",
        "idvg_curve": {"path": str(args.idvg_curve), "sha256": sha256(args.idvg_curve)},
        "files": [],
    }
    for name, builder in builders.items():
        source = args.bundle_dir / name
        rendered, states = builder(source.read_text(encoding="utf-8"))
        output = args.output_dir / name
        output.write_text(rendered, encoding="utf-8")
        manifest["files"].append({
            "name": name,
            "parent_sha256": sha256(source),
            "derived_sha256": sha256(output),
            "only_intended_change": "additional non-loadable Plot state capture statements",
            "states": states,
        })
    source_parameter = args.bundle_dir / "sdevice.par"
    target_parameter = args.output_dir / "sdevice.par"
    target_parameter.write_bytes(source_parameter.read_bytes())
    manifest["files"].append({
        "name": "sdevice.par",
        "parent_sha256": sha256(source_parameter),
        "derived_sha256": sha256(target_parameter),
        "only_intended_change": "none",
        "states": {},
    })
    write_json(args.output_dir / "state_deck_manifest.json", manifest)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
