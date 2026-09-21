#!/usr/bin/env python3
"""Run general-Tri3 Sentaurus avalanche controls on the coarse7x3 mesh."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.pn2d_general_tri3_contract import (
    EXACT_BIASES_V,
    SCHEMA_ID,
    SENTAURUS_RELEASE,
)
from scripts.sentaurus_avalanche_controls import (
    make_solve_block,
    make_tcl,
    validate_biases,
    validate_remote_root,
)


VARIANTS = {
    "implicit_default": {
        "avalanche": "Avalanche(VanOverstraeten)",
    },
    "explicit_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
    },
    "explicit_electric_field": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
    },
    "grad_qf_use_qf_contacts": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "use_qf_contacts": True,
    },
    "electric_field_use_qf_contacts": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
        "use_qf_contacts": True,
    },
    "grad_qf_aval_dens_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": True,
    },
    "lowfield_mobility_avalanche_electric_field": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
        "use_qf_contacts": True,
        "disable_hfs": True,
    },
    "lowfield_mobility_avalanche_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "use_qf_contacts": True,
        "disable_hfs": True,
    },
}


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
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-replay-20260725"
        ),
    )
    parser.add_argument("--case-name", default="coarse7x3")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--remote-root", required=True)
    parser.add_argument(
        "--biases",
        nargs="+",
        type=int,
        default=[int(value) for value in EXACT_BIASES_V],
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_ascii(path: Path, text: str) -> None:
    path.write_text(text, encoding="ascii", newline="\n")


def write_manifest(path: Path, value: dict[str, object]) -> None:
    write_ascii(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def validate_case_name(value: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9._-]+", value) is None:
        raise ValueError(f"unsafe case name: {value!r}")
    return value


def make_general_tcl(template: str, biases: tuple[int, ...]) -> str:
    result = make_tcl(template, biases)
    old = """\
    set eqfp [$data ReadScalar $::des_data_vertex "eQuasiFermiPotential"]
    set qfp0 [tcl_cp_get_double $eqfp 0]

    set target ""
"""
    new = """\
    set eqfp [$data ReadScalar $::des_data_vertex "eQuasiFermiPotential"]
    set qfp_min 1.0e100
    set qfp_vertex_count [$mesh size_vertex]
    for {set qfp_index 0} {$qfp_index < $qfp_vertex_count} {incr qfp_index} {
        set qfp_value [tcl_cp_get_double $eqfp $qfp_index]
        if {$qfp_value < $qfp_min} {
            set qfp_min $qfp_value
        }
    }

    set target ""
"""
    if old not in result:
        raise ValueError("QFP target probe block was not found")
    result = result.replace(old, new, 1)
    result = result.replace(
        "abs($qfp0-$candidate)",
        "abs($qfp_min-$candidate)",
        1,
    )
    return result


def make_variant_deck(
    template: str,
    variant: str,
    biases: tuple[int, ...],
) -> str:
    if variant not in VARIANTS:
        raise ValueError(f"unknown variant: {variant}")
    validate_biases(biases)
    spec = VARIANTS[variant]
    result = template.replace("pn2d_minimal6.tdr", "pn2d_msh.tdr")
    result = result.replace(
        "runtime_element_avalanche_probe_default",
        f"runtime_general_tri3_avalanche_probe_{variant}",
    )
    default_avalanche = "Avalanche(VanOverstraeten)"
    if result.count(default_avalanche) != 1:
        raise ValueError("default avalanche selector must occur exactly once")
    result = result.replace(
        default_avalanche,
        str(spec["avalanche"]),
        1,
    )
    math_options: list[str] = []
    if spec.get("use_qf_contacts"):
        math_options.append(
            "ComputeGradQuasiFermiAtContacts=UseQuasiFermi"
        )
    if spec.get("aval_dens_grad_qf"):
        math_options.append("AvalDensGradQF")
    if math_options:
        result, count = re.subn(
            r"Math\s*\{",
            "Math {\n  " + "\n  ".join(math_options),
            result,
            count=1,
        )
        if count != 1:
            raise ValueError("Math block was not found exactly once")
    if spec.get("disable_hfs"):
        result, count = re.subn(
            r"Mobility\s*\(\s*DopingDependence\s*"
            r"HighFieldSaturation\s*\)",
            "Mobility(\n    DopingDependence\n  )",
            result,
            count=1,
            flags=re.S,
        )
        if count != 1:
            raise ValueError("HighFieldSaturation block was not removed once")
    result, count = re.subn(
        r"Solve\s*\{.*\}\s*$",
        make_solve_block(biases),
        result,
        count=1,
        flags=re.S,
    )
    if count != 1:
        raise ValueError("Solve block was not replaced exactly once")
    return result.rstrip() + "\n"


def sentaurus_release(ssh_bin: Path, ssh_target: str) -> str:
    completed = subprocess.run(
        [str(ssh_bin), ssh_target, "sdevice -v 2>&1"],
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    match = re.search(r"Version\s+([^\s*]+)", output)
    if match is None:
        raise RuntimeError("Sentaurus release was not found")
    return match.group(1)


def main() -> int:
    args = parse_args()
    case_name = validate_case_name(args.case_name)
    biases = tuple(args.biases)
    validate_biases(biases)
    remote_root = validate_remote_root(args.remote_root)
    source_root = args.source_root.resolve()
    template_root = args.template_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "manifest.json"

    source_files = {
        "pn2d_msh.tdr": source_root / "pn2d_msh.tdr",
        "models.par": source_root / "models.par",
    }
    for name, path in source_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing {case_name} source {name}: {path}")

    release = sentaurus_release(args.ssh_bin, args.ssh_target)
    if release != SENTAURUS_RELEASE:
        raise RuntimeError(
            f"expected Sentaurus {SENTAURUS_RELEASE}, got {release}"
        )

    tcl_template = (
        template_root / "runtime_element_avalanche_probe.tcl"
    ).read_text(encoding="ascii")
    deck_template = (
        template_root / "runtime_element_avalanche_probe_default.cmd"
    ).read_text(encoding="ascii")
    tcl_text = make_general_tcl(tcl_template, biases)

    static_hashes = {
        "tdr": sha256(source_files["pn2d_msh.tdr"]),
        "models.par": sha256(source_files["models.par"]),
    }
    manifest: dict[str, object] = {
        "schema": SCHEMA_ID,
        "status": "running",
        "experiment": "pn2d_general_tri3_sentaurus_avalanche_controls",
        "sentaurus_release": release,
        "ssh_target": args.ssh_target,
        "remote_root": remote_root,
        "exact_biases_V": list(biases),
        "variants": list(VARIANTS),
        "case_hashes": {case_name: static_hashes},
        "cases": {case_name: {}},
    }
    write_manifest(manifest_path, manifest)

    case_results: dict[str, object] = {}
    for variant in VARIANTS:
        local = output_root / case_name / variant
        bundle = local / "bundle"
        fetched = local / "fetched"
        bundle.mkdir(parents=True, exist_ok=True)
        fetched.mkdir(parents=True, exist_ok=True)
        for name, source in source_files.items():
            shutil.copy2(source, bundle / name)
        tcl_name = "runtime_element_avalanche_probe.tcl"
        write_ascii(bundle / tcl_name, tcl_text)
        deck_name = f"runtime_general_tri3_avalanche_probe_{variant}.cmd"
        write_ascii(
            bundle / deck_name,
            make_variant_deck(deck_template, variant, biases),
        )

        remote = f"{remote_root}/{case_name}/{variant}"
        run([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote}"])
        for name in ("pn2d_msh.tdr", "models.par", tcl_name, deck_name):
            run(
                [
                    str(args.scp_bin),
                    str(bundle / name),
                    f"{args.ssh_target}:{remote}/",
                ]
            )
        run_name = f"run_{variant}.out"
        run(
            [
                str(args.ssh_bin),
                args.ssh_target,
                f"cd {remote} && sdevice {deck_name} > {run_name} 2>&1",
            ]
        )
        stem = f"runtime_general_tri3_avalanche_probe_{variant}"
        fetched_names = (run_name, f"{stem}.plt", f"{stem}_des.log")
        for name in fetched_names:
            run(
                [
                    str(args.scp_bin),
                    f"{args.ssh_target}:{remote}/{name}",
                    str(fetched / name),
                ]
            )
        text = (fetched / run_name).read_text(
            encoding="ascii",
            errors="strict",
        )
        observed = sorted(
            {
                int(value)
                for value in re.findall(
                    r"AVAL_PROBE_BEGIN bias_V=(-?\d+)",
                    text,
                )
            },
            reverse=True,
        )
        status = "passed" if tuple(observed) == biases else "failed"
        case_results[variant] = {
            "status": status,
            "observed_biases_V": observed,
            "bundle_sha256": {
                name: sha256(bundle / name)
                for name in ("pn2d_msh.tdr", "models.par", tcl_name, deck_name)
            },
            "output_sha256": {
                name: sha256(fetched / name) for name in fetched_names
            },
        }
        manifest["cases"][case_name] = case_results
        write_manifest(manifest_path, manifest)
        if status != "passed":
            raise RuntimeError(
                f"{variant}: expected {biases}, got {tuple(observed)}"
            )

    manifest["status"] = "passed"
    write_manifest(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
