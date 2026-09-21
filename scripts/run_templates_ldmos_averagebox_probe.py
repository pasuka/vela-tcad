#!/usr/bin/env python3
"""Prepare and optionally run an LDMOS AverageBox diagnostic oracle.

The probe preserves the archived exact SProcess grid and G3-no-IALMob physics,
adds the documented box debug switches, and performs either one Poisson
initialization or one coupled Newton iteration from an archived state.  It
never edits the Applications Library or the archived parent run in place.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


REMOTE_ROOT = "/root/sentaurus_runs/vela_oracle_2022/templates_ldmos"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32" / "OpenSSH" / f"{name}.exe"
        )
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], *, check: bool = True) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv), check=check, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )


def render_deck(
    source: str,
    *,
    initial_state_prefix: str | None = None,
    gate_voltage_V: float = 0.0,
) -> str:
    if source.count("Math {") != 1 or source.count("Solve {") != 1:
        raise ValueError("expected one Math and one Solve section")
    if "AverageBoxMethod" in source or "BoxMeasureFromFile" in source:
        raise ValueError("source deck already contains a box debug policy")
    result, count = re.subn(
        r'(?m)^(\s*\{\s*name=\s*"gate".*?Voltage=\s*)[-+0-9.eE]+',
        rf"\g<1>{gate_voltage_V:.17g}",
        source,
        count=1,
    )
    if count != 1 and gate_voltage_V != 0.0:
        raise ValueError("could not set the gate voltage in the source deck")
    result = result.replace(
        "Math {",
        "Math {\n"
        "\tAverageBoxMethod\n"
        "\tBoxMeasureFromFile(GrdNumbering)\n"
        "\tBoxCoefficientsFromFile(GrdNumbering)",
        1,
    )
    result, count = re.subn(
        r"(?m)^Plot\s*\{",
        "Plot {\n"
        "\tBM_AngleElements\n"
        "\tBM_CoeffIntersectionNonDelaunayElements\n"
        "\tBM_ElementVolume\n"
        "\tBM_IntersectionNonDelaunayElements\n"
        "\tBM_VolumeIntersectionNonDelaunayElements",
        result,
        count=1,
    )
    if count != 1:
        raise ValueError("expected one top-level Plot section")
    solve = result.index("Solve {")
    load = (
        f'\tLoad(FilePrefix="{initial_state_prefix}")\n'
        if initial_state_prefix else ""
    )
    equations = "Poisson Electron Hole" if initial_state_prefix else "Poisson"
    return result[:solve] + (
        "Solve {\n" + load
        + f"\tCoupled(Iterations=1) {{ {equations} }}\n"
        '\tPlot(FilePrefix="ldmos_averagebox_probe")\n'
        "}\n"
    )


def render_state_capture_deck(source: str, gate_voltage_V: float) -> str:
    """Render the source G3 sweep to a loadable state at one gate voltage."""
    result, count = re.subn(
        r'(Goal\s*\{\s*name=\s*"gate"\s+voltage=\s*)[-+0-9.eE]+',
        rf"\g<1>{gate_voltage_V:.17g}",
        source,
        count=1,
        flags=re.IGNORECASE,
    )
    if count != 1:
        raise ValueError("expected one gate Goal in the source deck")
    solve = result.index("Solve {")
    closing = result.rfind("}")
    if closing <= solve:
        raise ValueError("could not locate the closing Solve brace")
    return (
        result[:closing]
        + '\n\tSave(FilePrefix="initial_state")\n'
        + result[closing:]
    )


def safe_extract(archive: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    root = output.resolve()
    with tarfile.open(archive, "r:gz") as stream:
        for member in stream.getmembers():
            target = (output / member.name).resolve()
            if target != root and root not in target.parents:
                raise ValueError(f"archive member escapes output: {member.name}")
        stream.extractall(output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--source-deck", default="IdVg.cmd")
    parser.add_argument(
        "--initial-tdr", type=Path,
        help="optional loadable Sentaurus state for a one-iteration coupled probe",
    )
    parser.add_argument(
        "--capture-loadable-state", action="store_true",
        help="first run the source sweep to create a loadable state at the gate voltage",
    )
    parser.add_argument("--gate-voltage-V", type=float, default=0.0)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite probe: {output}")
    bundle = output / "bundle"
    bundle.mkdir(parents=True)
    inputs = []
    for name in ("n1_fps.tdr", "sdevice.par"):
        src = source / name
        if not src.is_file():
            raise FileNotFoundError(src)
        dst = bundle / name
        shutil.copy2(src, dst)
        inputs.append({"path": name, "sha256": sha256(dst)})
    source_deck = source / args.source_deck
    if args.initial_tdr is not None and args.capture_loadable_state:
        raise ValueError("--initial-tdr and --capture-loadable-state are exclusive")
    initial_prefix = None
    if args.initial_tdr is not None:
        initial_tdr = args.initial_tdr.resolve()
        if not initial_tdr.is_file():
            raise FileNotFoundError(initial_tdr)
        initial_prefix = "initial_state"
        initial_copy = bundle / f"{initial_prefix}_des.tdr"
        shutil.copy2(initial_tdr, initial_copy)
        inputs.append({"path": initial_copy.name, "sha256": sha256(initial_copy)})
    elif args.capture_loadable_state:
        initial_prefix = "initial_state"
        capture = bundle / "capture_state.cmd"
        capture.write_text(
            render_state_capture_deck(
                source_deck.read_text(encoding="utf-8"), args.gate_voltage_V
            ),
            encoding="utf-8",
        )
        inputs.append({"path": capture.name, "sha256": sha256(capture)})
    deck = bundle / "averagebox_probe.cmd"
    deck.write_text(
        render_deck(
            source_deck.read_text(encoding="utf-8"),
            initial_state_prefix=initial_prefix,
            gate_voltage_V=args.gate_voltage_V,
        ),
        encoding="utf-8",
    )
    inputs.append({"path": deck.name, "sha256": sha256(deck)})

    run_id = "averagebox_node4492_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest: dict[str, object] = {
        "schema": "vela.templates_ldmos.averagebox_probe.v1",
        "classification": "derived_geometry_oracle_not_official_curve",
        "source_deck": str(source_deck),
        "source_deck_sha256": sha256(source_deck),
        "inputs": inputs,
        "box_policy": "AverageBoxMethod",
        "debug_contract": [
            "BoxMeasureFromFile(GrdNumbering)",
            "BoxCoefficientsFromFile(GrdNumbering)",
        ],
        "solve_scope": (
            "one_coupled_newton_iteration_from_archived_state"
            if initial_prefix else "one_poisson_initialization"
        ),
        "gate_voltage_V": args.gate_voltage_V,
        "state_source": (
            "captured_in_same_remote_run" if args.capture_loadable_state
            else ("provided_loadable_tdr" if args.initial_tdr else "none")
        ),
        "live": args.live,
        "run_id": run_id,
    }
    if args.live:
        remote = f"{REMOTE_ROOT}/{run_id}"
        quoted = shlex.quote(remote)
        run([args.ssh_bin, args.ssh_target,
             f"umask 077; test ! -e {quoted}; mkdir -p {quoted}"])
        run([
            args.scp_bin, *(str(path) for path in sorted(bundle.iterdir())),
            f"{args.ssh_target}:{remote}/",
        ])
        capture_command = (
            "sdevice capture_state.cmd > capture.out 2>&1; capture_rc=$?; "
            "printf '%s\\n' \"$capture_rc\" > capture_exitcode; "
            "if test \"$capture_rc\" -ne 0; then exit \"$capture_rc\"; fi; "
            if args.capture_loadable_state else ""
        )
        command = (
            f"cd {quoted}; rm -f MeasureCoefficients.debug; "
            + capture_command
            + "sdevice averagebox_probe.cmd > run.out 2>&1; rc=$?; "
            "printf '%s\\n' \"$rc\" > exitcode; "
            "tar -czf results.tgz averagebox_probe.cmd MeasureCoefficients.debug "
            "run.out exitcode capture.out capture_exitcode initial_state_des.sav "
            "initial_state_circuit_des.sav "
            "ldmos_averagebox_probe_des.tdr *_1_des.tdr n2_des.log 2>/dev/null; "
            "exit $rc"
        )
        completed = run([args.ssh_bin, args.ssh_target, command], check=False)
        archive = output / "results.tgz"
        run([args.scp_bin, f"{args.ssh_target}:{remote}/results.tgz", str(archive)])
        raw = output / "raw"
        safe_extract(archive, raw)
        manifest.update({
            "remote_directory": remote,
            "return_code": completed.returncode,
            "ssh_output": completed.stdout,
            "artifacts": [
                {
                    "path": path.relative_to(output).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in sorted(raw.iterdir()) if path.is_file()
            ],
        })
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    # A one-iteration diagnostic normally reaches the configured iteration
    # cap.  The run is a successful acquisition when the requested box debug
    # oracle was generated; retain the raw solver return code in the manifest.
    if args.live and not (output / "raw" / "MeasureCoefficients.debug").is_file():
        return int(manifest.get("return_code", 1)) or 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
