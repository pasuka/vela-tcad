#!/usr/bin/env python3
"""Refresh BJT VTK recombination fields and dependent comparisons from saved states."""

from __future__ import annotations

import argparse
import copy
import json
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
INDICES = (0, 10, 20, 30)


def run(command: list[str], *, allowed_returncodes: tuple[int, ...] = (0,)) -> None:
    completed = subprocess.run(
        command, cwd=REPO, text=True, encoding="utf-8", errors="replace",
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if completed.returncode not in allowed_returncodes:
        raise RuntimeError(completed.stderr or completed.stdout)


def export_state(
    *, runner: Path, base: dict[str, object], state: Path, output_vtk: Path,
    mesh: Path, doping: Path, bias: float, config_path: Path,
) -> None:
    config = copy.deepcopy(base)
    config.pop("sweep", None)
    config.update(
        {
            "simulation_type": "write_dd_state_vtk",
            "mesh_file": str(mesh.resolve()),
            "node_doping_file": str(doping.resolve()),
            "materials_file": str(
                (FIXTURE / "vela" / "materials_sentaurus2022.json").resolve()
            ),
            "state_file": str(state.resolve()),
            "output_vtk": str(output_vtk.resolve()),
        }
    )
    for contact in config["contacts"]:
        if contact["name"] == "collector":
            contact["bias"] = bias
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    run([str(runner.resolve()), "--config", str(config_path.resolve())])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runner", type=Path, default=REPO / "build-release" / "vela_example_runner.exe"
    )
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    args = parser.parse_args()

    base = json.loads(
        (FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json").read_text(
            encoding="utf-8"
        )
    )
    accepted = BUILD_ROOT / "m1_p0_acceptance"
    configs = accepted / "postprocess_refresh_configs"
    for index in INDICES:
        token = f"vce_{index:03d}"
        export_state(
            runner=args.runner,
            base=base,
            state=accepted / "states" / f"{token}.csv",
            output_vtk=accepted / "fields" / f"{token}.vtk",
            mesh=FIXTURE / "vela" / "input" / "mesh.json",
            doping=FIXTURE / "vela" / "input" / "doping.csv",
            bias=index / 10.0,
            config_path=configs / f"{token}.json",
        )

    refined = BUILD_ROOT / "mesh_sensitivity" / "local_refined"
    refined_base = json.loads(
        (refined / "vela" / "configs" / "m1_spatial_vce3.json").read_text(
            encoding="utf-8"
        )
    )
    export_state(
        runner=args.runner,
        base=refined_base,
        state=refined / "vela" / "m1_spatial_vce3_state.csv",
        output_vtk=refined / "vela" / "m1_spatial_vce3.vtk",
        mesh=refined / "vela" / "input" / "mesh.json",
        doping=refined / "vela" / "input" / "doping.csv",
        bias=3.0,
        config_path=refined / "vela" / "configs" / "m1_spatial_vce3_postprocess_refresh.json",
    )

    sentaurus_root = BUILD_ROOT / "m1_p0_multibias" / "sentaurus"
    transport_root = BUILD_ROOT / "m1_p0_multibias" / "transport_sources"
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "export_genius_bjt_accepted_transport_sources.py"),
            "--accepted-root", str(accepted),
            "--mesh-root", str(sentaurus_root / "vce_000"),
            "--output-root", str(transport_root),
            "--sentaurus-root", str(sentaurus_root),
        ],
        allowed_returncodes=(0, 1),
    )
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "compare_genius_bjt_transport_fields.py"),
            "--reference-root", str(FIXTURE),
            "--sentaurus-fields-root", str(sentaurus_root / "vce_030"),
            "--vela-vtk", str(accepted / "fields" / "vce_030.vtk"),
            "--output-dir", str(BUILD_ROOT / "m1_p0_relative_error" / "transport"),
        ],
        allowed_returncodes=(0, 1),
    )
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "compare_genius_bjt_transport_fields.py"),
            "--reference-root", str(FIXTURE),
            "--sentaurus-fields-root", str(refined / "sentaurus_exports" / "vce3"),
            "--vela-vtk", str(refined / "vela" / "m1_spatial_vce3.vtk"),
            "--output-dir", str(refined / "comparison" / "transport"),
        ]
    )
    run(
        [
            str(args.python.resolve()),
            str(REPO / "scripts" / "audit_genius_bjt_conservative_sections.py"),
        ]
    )
    print(json.dumps({"refreshed_biases_V": [0.0, 1.0, 2.0, 3.0], "refined_VCE_V": 3.0}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
