#!/usr/bin/env python3
"""Run the exact coarse7x3 Sentaurus high-bias process-variable probe."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.pn2d_high_bias_process_contract import SENTAURUS_RELEASE
from scripts.pn2d_sentaurus_process_run_contract import build_run_manifest
from scripts.run_pn2d_general_tri3_sentaurus_avalanche_controls_vm import (
    make_general_tcl,
    make_variant_deck,
)
from scripts.sentaurus_avalanche_controls import (
    validate_remote_root,
)


BIAS_V = -19.95
VARIANT = "implicit_default"
PROCESS_FIELDS = (
    "eVelocity",
    "hVelocity",
    "TotalCurrentDensity/Vector",
    "eIonIntegral",
    "hIonIntegral",
    "MeanIonIntegral",
    "Doping",
    "DonorConcentration",
    "AcceptorConcentration",
    "SpaceCharge",
    "SRHRecombination",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
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
    parser.add_argument(
        "--source-root",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--template-root",
        type=Path,
        default=Path("build-release/pn2d-minimal6-element-avalanche-replay-20260725"),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--remote-root", required=True)
    return parser.parse_args()


def write_ascii(path: Path, text: str) -> None:
    path.write_text(text, encoding="ascii", newline="\n")


def process_deck(template: str) -> str:
    deck = make_variant_deck(template, VARIANT, (BIAS_V,))
    anchor = "  hAlphaAvalanche\n}"
    if deck.count(anchor) != 1:
        raise ValueError("process Plot anchor must occur exactly once")
    replacement = "\n".join(("  hAlphaAvalanche", *("  " + name for name in PROCESS_FIELDS), "}"))
    return deck.replace(anchor, replacement, 1)


def process_tcl(template: str) -> str:
    tcl = make_general_tcl(template, (BIAS_V,))
    return tcl.replace("bias_V=%d", "bias_V=%.17g")


def sentaurus_release(ssh_bin: Path, ssh_target: str) -> str:
    result = subprocess.run(
        [str(ssh_bin), ssh_target, "sdevice -v 2>&1"],
        check=False,
        capture_output=True,
        text=True,
    )
    match = re.search(r"Version\s+([^\s*]+)", result.stdout + result.stderr)
    if match is None:
        raise RuntimeError("Sentaurus release was not found")
    return match.group(1)


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()
    remote_root = validate_remote_root(args.remote_root)
    output = args.output_root.resolve()
    bundle = output / "bundle"
    fetched = output / "fetched"
    bundle.mkdir(parents=True, exist_ok=True)
    fetched.mkdir(parents=True, exist_ok=True)

    release = sentaurus_release(args.ssh_bin, args.ssh_target)
    if release != SENTAURUS_RELEASE:
        raise RuntimeError(f"expected Sentaurus {SENTAURUS_RELEASE}, got {release}")

    source_root = args.source_root.resolve()
    template_root = args.template_root.resolve()
    sources = {
        "pn2d_msh.tdr": source_root / "pn2d_msh.tdr",
        "models.par": source_root / "models.par",
    }
    for name, source in sources.items():
        if not source.is_file():
            raise FileNotFoundError(f"missing probe input {name}: {source}")
        shutil.copy2(source, bundle / name)

    tcl_name = "runtime_element_avalanche_probe.tcl"
    deck_name = "runtime_general_tri3_avalanche_probe_implicit_default.cmd"
    write_ascii(
        bundle / tcl_name,
        process_tcl((template_root / tcl_name).read_text(encoding="ascii")),
    )
    write_ascii(
        bundle / deck_name,
        process_deck(
            (template_root / "runtime_element_avalanche_probe_default.cmd").read_text(
                encoding="ascii"
            )
        ),
    )

    remote = f"{remote_root}/coarse7x3/{VARIANT}"
    run([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote}"])
    for name in (*sources, tcl_name, deck_name):
        run([str(args.scp_bin), str(bundle / name), f"{args.ssh_target}:{remote}/"])
    run_name = "run_implicit_default.out"
    run(
        [
            str(args.ssh_bin),
            args.ssh_target,
            f"cd {remote} && sdevice {deck_name} > {run_name} 2>&1",
        ]
    )

    stem = "runtime_general_tri3_avalanche_probe_implicit_default"
    fetched_names = (
        run_name,
        f"{stem}.plt",
        f"{stem}.tdr",
        f"{stem}_des.log",
    )
    for name in fetched_names:
        run(
            [
                str(args.scp_bin),
                f"{args.ssh_target}:{remote}/{name}",
                str(fetched / name),
            ]
        )

    run_text = (fetched / run_name).read_text(encoding="ascii")
    observed = [
        float(value)
        for value in re.findall(
            r"AVAL_PROBE_BEGIN bias_V=(-?\d+(?:\.\d+)?)",
            run_text,
        )
    ]
    exact = len(observed) == 1 and abs(observed[0] - BIAS_V) <= 1.0e-12
    manifest = build_run_manifest(
        status="passed" if exact else "failed",
        experiment="pn2d_high_bias_process_probe",
        variant=VARIANT,
        exact_biases=(BIAS_V,),
        observed_biases=observed,
        remote_root=remote_root,
        bundle=bundle,
        fetched=fetched,
    )
    write_ascii(
        output / "manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    if not exact:
        raise RuntimeError(f"exact {BIAS_V:g} V callback was not observed: {observed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
