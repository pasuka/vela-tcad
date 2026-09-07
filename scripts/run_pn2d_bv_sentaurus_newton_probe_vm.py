#!/usr/bin/env python3
"""Extract fixed-transition Sentaurus residual and first-Newton-update TDRs."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_pn2d_bv_process_matrix_vm import (
    copy_remote_files,
    list_remote_files,
    run_checked,
)
from scripts.sentaurus_avalanche_controls import (
    validate_remote_root,
)


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
KNEE_BIASES = (-19.7, -19.8, -19.85, -19.9, -19.95, -20.0)
EXACT_TOLERANCE_V = 1.0e-10


def parse_branches(text: str) -> tuple[str, ...]:
    result = tuple(item.strip() for item in text.split(",") if item.strip())
    if not result or len(set(result)) != len(result):
        raise ValueError("--branches must contain unique branch names")
    unknown = sorted(set(result) - set(BRANCHES))
    if unknown:
        raise ValueError(f"unknown branches: {', '.join(unknown)}")
    return result


def bias_token(value: float) -> str:
    return ("m" if value < 0.0 else "p") + f"{abs(value):.6f}".replace(".", "p")


def _replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    result, count = re.subn(
        pattern, replacement, text, count=1, flags=re.S | re.I | re.M
    )
    if count != 1:
        raise ValueError(f"{label} was not found exactly once")
    return result


def make_probe_deck(
    template: str,
    *,
    source_prefix: str,
    target_bias: float,
    stem: str,
) -> str:
    """Build a one-iteration fixed-transition probe from a loadable exact state."""

    deck = template
    replacements = {
        "Plot": f'  Plot      = "{stem}.tdr"',
        "Current": f'  Current   = "{stem}.plt"',
        "Output": f'  Output    = "{stem}"',
    }
    for key, line in replacements.items():
        deck = _replace_once(
            deck,
            rf"^\s*{key}\s*=\s*\"[^\"]*\"\s*$",
            line,
            f"File.{key}",
        )
    newton_line = f'  NewtonPlot = "{stem}_newton_%d_%d_des.tdr"'
    if re.search(r"^\s*NewtonPlot\s*=", deck, flags=re.M | re.I):
        deck = _replace_once(
            deck,
            r"^\s*NewtonPlot\s*=\s*\"[^\"]*\"\s*$",
            newton_line,
            "File.NewtonPlot",
        )
    else:
        output_line = f'  Output    = "{stem}"'
        deck = deck.replace(output_line, output_line + "\n" + newton_line, 1)

    deck = re.sub(
        r"^\s*(?:CNormPrint|NewtonPlot\s*\([^)]*\)|Extrapolate|"
        r"AutoNPMinStepFactor\s*=\s*\S+|AutoCNPMinStepFactor\s*=\s*\S+)\s*$",
        "",
        deck,
        flags=re.M | re.I,
    )
    math_anchor = re.search(r"Math\s*\{", deck, flags=re.I)
    if math_anchor is None:
        raise ValueError("Math block is missing")
    controls = (
        "\n  CNormPrint"
        "\n  NewtonPlot(Error Residual Update)"
        "\n  AutoNPMinStepFactor=0"
        "\n  AutoCNPMinStepFactor=0"
    )
    deck = deck[: math_anchor.end()] + controls + deck[math_anchor.end() :]

    solve = f"""\
Solve {{
  Load(FilePrefix="{source_prefix}")
  Quasistationary(
    InitialStep=1 MinStep=0.9 MaxStep=1
    Increment=1 Decrement=2 NewtonPlotStep=2
    Goal {{ Name="Anode" Voltage={target_bias:.17g} }}
  ) {{
    Coupled(Iterations=1) {{ Poisson Electron Hole }}
  }}
}}
"""
    deck = _replace_once(deck, r"Solve\s*\{.*\}\s*$", solve.rstrip(), "Solve block")
    return deck.rstrip() + "\n"


def branch_bias_records(manifest: dict[str, Any], branch: str) -> list[dict[str, Any]]:
    record = next(
        (item for item in manifest["branch_records"] if item["branch"] == branch),
        None,
    )
    if record is None:
        raise ValueError(f"manifest is missing {branch}")
    return list(record["bias_records"])


def transition_sources(
    manifest: dict[str, Any],
    branch: str,
    targets: tuple[float, ...],
) -> list[tuple[float, dict[str, Any]]]:
    records = branch_bias_records(manifest, branch)
    requested = [float(record["requested_bias_V"]) for record in records]
    result: list[tuple[float, dict[str, Any]]] = []
    for target in targets:
        matches = [
            index
            for index, value in enumerate(requested)
            if abs(value - target) <= EXACT_TOLERANCE_V
        ]
        if len(matches) != 1 or matches[0] == 0:
            raise ValueError(f"{branch}: target {target:g} has no unique predecessor")
        target_index = matches[0]
        result.append((target, records[target_index - 1]))
    return result


def source_prefix_from_artifact(remote_source_root: str, artifact: dict[str, str]) -> str:
    name = Path(artifact["path"]).name
    suffix = "_des.tdr"
    if not name.endswith(suffix):
        raise ValueError(f"unexpected snapshot name: {name}")
    return f"{remote_source_root}/{artifact['path'][:-len(suffix)]}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--remote-source-root", required=True)
    parser.add_argument("--remote-root", required=True)
    parser.add_argument("--branches", default=",".join(BRANCHES))
    parser.add_argument("--biases", nargs="+", type=float, default=list(KNEE_BIASES))
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument(
        "--ssh-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\ssh.exe"),
    )
    parser.add_argument(
        "--scp-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\scp.exe"),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    branches = parse_branches(args.branches)
    targets = tuple(float(value) for value in args.biases)
    if len(set(targets)) != len(targets):
        raise ValueError("--biases contains duplicates")
    remote_source_root = validate_remote_root(args.remote_source_root)
    remote_root = validate_remote_root(args.remote_root)
    source_manifest_path = args.source_manifest.resolve()
    source_root = source_manifest_path.parent
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))

    cases: list[dict[str, Any]] = []
    for branch in branches:
        template_path = source_root / branch / "raw" / f"pn2d_bv_process_{branch}.cmd"
        tcl_path = source_root / branch / "raw" / "runtime_element_avalanche_probe.tcl"
        grid_path = source_root / branch / "raw" / "pn2d_msh.tdr"
        parameter_path = source_root / branch / "raw" / "models.par"
        template = template_path.read_text(encoding="ascii")
        for target, predecessor in transition_sources(manifest, branch, targets):
            token = bias_token(target)
            case_name = f"{branch}_{token}"
            case = output / branch / token
            bundle = case / "bundle"
            raw = case / "raw"
            bundle.mkdir(parents=True, exist_ok=True)
            stem = f"newton_probe_{branch}_{token}"
            source_prefix = source_prefix_from_artifact(
                remote_source_root, predecessor["snapshot_tdr"]
            )
            deck_path = bundle / f"{stem}.cmd"
            deck_path.write_text(
                make_probe_deck(
                    template,
                    source_prefix=source_prefix,
                    target_bias=target,
                    stem=stem,
                ),
                encoding="ascii",
                newline="\n",
            )
            for path in (tcl_path, grid_path, parameter_path):
                shutil.copy2(path, bundle / path.name)
            cases.append(
                {
                    "case_name": case_name,
                    "branch": branch,
                    "target_bias_V": target,
                    "predecessor_bias_V": float(predecessor["requested_bias_V"]),
                    "source_prefix": source_prefix,
                    "stem": stem,
                    "case": case,
                    "bundle": bundle,
                    "raw": raw,
                }
            )

    if args.dry_run:
        print(json.dumps({"cases": len(cases), "root": str(output)}, indent=2))
        return 0

    staging = output / "_remote_raw"
    local_complete = all(
        len(list(item["raw"].glob(f"{item['stem']}_newton_1_*_des.tdr"))) == 1
        for item in cases
    )
    if args.resume and not local_complete:
        names = list_remote_files(args.ssh_bin, args.ssh_target, remote_root)
        retained = [
            name
            for name in names
            if any(name.startswith(item["stem"]) for item in cases)
            and name.endswith((".tdr", ".log", ".cmd", ".plt", ".out"))
        ]
        copy_remote_files(
            args.scp_bin, args.ssh_target, remote_root, retained, staging
        )
        for item in cases:
            raw = item["raw"]
            raw.mkdir(parents=True, exist_ok=True)
            for path in staging.glob(f"{item['stem']}*"):
                shutil.copy2(path, raw / path.name)

    pending = [
        item
        for item in cases
        if not (
            args.resume
            and len(list(item["raw"].glob(f"{item['stem']}_newton_1_*_des.tdr")))
            == 1
        )
    ]
    if pending:
        run_checked([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote_root}"])
        common_names = {
            "pn2d_msh.tdr",
            "models.par",
            "runtime_element_avalanche_probe.tcl",
        }
        upload = [
            path
            for path in sorted(pending[0]["bundle"].iterdir())
            if path.name in common_names
        ]
        upload.extend(item["bundle"] / f"{item['stem']}.cmd" for item in pending)
        run_checked(
            [
                str(args.scp_bin),
                *[str(path) for path in upload],
                f"{args.ssh_target}:{remote_root}/",
            ]
        )
        commands = [
            (
                f"sdevice {shlex.quote(item['stem'] + '.cmd')} "
                f"> {shlex.quote(item['stem'] + '.run.out')} 2>&1"
            )
            for item in pending
        ]
        batch = f"cd {shlex.quote(remote_root)} && " + "; ".join(commands)
        subprocess.run(
            [str(args.ssh_bin), args.ssh_target, batch],
            check=False,
        )
        names = list_remote_files(args.ssh_bin, args.ssh_target, remote_root)
        retained = [
            name
            for name in names
            if any(name.startswith(item["stem"]) for item in pending)
            and name.endswith((".tdr", ".log", ".cmd", ".plt", ".out"))
        ]
        copy_remote_files(
            args.scp_bin, args.ssh_target, remote_root, retained, staging
        )
        for item in pending:
            raw = item["raw"]
            raw.mkdir(parents=True, exist_ok=True)
            for path in staging.glob(f"{item['stem']}*"):
                shutil.copy2(path, raw / path.name)

    for item in cases:
        first_updates = list(item["raw"].glob(f"{item['stem']}_newton_1_*_des.tdr"))
        if len(first_updates) != 1:
            raise RuntimeError(
                f"{item['case_name']}: expected one first-update TDR, "
                f"got {len(first_updates)}"
            )

    execution = {
        "schema": "vela.pn2d_bv_sentaurus_newton_probe_execution.v1",
        "status": "passed",
        "outcome": "fixed_transition_first_newton_observations_available",
        "source_manifest": str(source_manifest_path),
        "remote_source_root": remote_source_root,
        "remote_root": remote_root,
        "cases": [
            {
                key: value
                for key, value in item.items()
                if key not in {"case", "bundle", "raw"}
            }
            | {
                "raw": str(item["raw"]),
                "first_update_tdr": str(
                    next(item["raw"].glob(f"{item['stem']}_newton_1_*_des.tdr"))
                ),
            }
            for item in cases
        ],
    }
    (output / "execution.json").write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
