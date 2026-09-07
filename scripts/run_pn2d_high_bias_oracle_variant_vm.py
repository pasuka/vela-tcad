#!/usr/bin/env python3
"""Run one exact-lattice coarse7x3 Sentaurus high-bias oracle variant."""

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

from scripts.pn2d_high_bias_process_contract import (
    EXACT_HIGH_BIAS_V,
    SENTAURUS_RELEASE,
)
from scripts.pn2d_sentaurus_process_run_contract import build_run_manifest
from scripts.run_pn2d_general_tri3_sentaurus_avalanche_controls_vm import (
    VARIANTS,
    make_general_tcl,
    make_variant_deck,
)
from scripts.run_pn2d_high_bias_process_probe_vm import (
    PROCESS_FIELDS,
    run,
    sentaurus_release,
    write_ascii,
)
from scripts.sentaurus_avalanche_controls import (
    validate_biases,
    validate_remote_root,
)


AVALANCHE_DISABLED = "avalanche_disabled"
ORACLE_VARIANTS = (*VARIANTS, AVALANCHE_DISABLED)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", choices=ORACLE_VARIANTS, required=True)
    parser.add_argument(
        "--biases",
        nargs="+",
        type=float,
        default=list(EXACT_HIGH_BIAS_V),
    )
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


def oracle_deck(template: str, variant: str, biases: tuple[float, ...]) -> str:
    builder_variant = "implicit_default" if variant == AVALANCHE_DISABLED else variant
    deck = make_variant_deck(template, builder_variant, biases)
    if variant == AVALANCHE_DISABLED:
        old_stem = "runtime_general_tri3_avalanche_probe_implicit_default"
        deck = deck.replace(
            old_stem,
            "runtime_general_tri3_avalanche_probe_avalanche_disabled",
        )
        deck, count = re.subn(
            r"\n\s*Avalanche\(VanOverstraeten\)\n",
            "\n",
            deck,
            count=1,
        )
        if count != 1:
            raise ValueError("avalanche term was not removed exactly once")
    anchor = "  hAlphaAvalanche\n}"
    if deck.count(anchor) != 1:
        raise ValueError("process Plot anchor must occur exactly once")
    fields = "\n".join(
        ("  hAlphaAvalanche", *("  " + name for name in PROCESS_FIELDS), "}")
    )
    return deck.replace(anchor, fields, 1)


def oracle_tcl(template: str, biases: tuple[float, ...]) -> str:
    tcl = make_general_tcl(template, biases).replace(
        "bias_V=%d",
        "bias_V=%.17g",
    )
    read_anchor = (
        '    set generation_total [$data ReadScalar $::des_data_vertex '
        '"AvalancheGeneration"]\n'
    )
    reads = read_anchor + """\
    set velocity_n [$data ReadScalar $::des_data_vertex "eVelocity"]
    set velocity_p [$data ReadScalar $::des_data_vertex "hVelocity"]
    set total_current [$data ReadVector $::des_data_vertex "TotalCurrentDensity"]
    set ion_n [$data ReadScalar $::des_data_vertex "eIonIntegral"]
    set ion_p [$data ReadScalar $::des_data_vertex "hIonIntegral"]
    set ion_mean [$data ReadScalar $::des_data_vertex "MeanIonIntegral"]
    set doping [$data ReadScalar $::des_data_vertex "DopingConcentration"]
    set space_charge [$data ReadScalar $::des_data_vertex "SpaceCharge"]
    set srh [$data ReadScalar $::des_data_vertex "srhRecombination"]
"""
    if tcl.count(read_anchor) != 1:
        raise ValueError("runtime process-read anchor must occur exactly once")
    tcl = tcl.replace(read_anchor, reads, 1)
    vertex_anchor = """\
            [tcl_cp_get_double $generation_total $vertex_index]]
    }
"""
    process_put = """\
            [tcl_cp_get_double $generation_total $vertex_index]]
        puts [format "AVAL_PROBE_PROCESS bias_V=%.17g vertex=%d velocity_n_cm_s=%.17g velocity_p_cm_s=%.17g total_current_x_A_cm2=%.17g total_current_y_A_cm2=%.17g ion_n=%.17g ion_p=%.17g ion_mean=%.17g doping_cm3=%.17g space_charge_cm3=%.17g srh_cm3_s=%.17g" \
            $target $vertex_index \
            [tcl_cp_get_double $velocity_n $vertex_index] \
            [tcl_cp_get_double $velocity_p $vertex_index] \
            [tcl_cp_get_double2 $total_current 0 $vertex_index] \
            [tcl_cp_get_double2 $total_current 1 $vertex_index] \
            [tcl_cp_get_double $ion_n $vertex_index] \
            [tcl_cp_get_double $ion_p $vertex_index] \
            [tcl_cp_get_double $ion_mean $vertex_index] \
            [tcl_cp_get_double $doping $vertex_index] \
            [tcl_cp_get_double $space_charge $vertex_index] \
            [tcl_cp_get_double $srh $vertex_index]]
    }
"""
    if tcl.count(vertex_anchor) != 1:
        raise ValueError("runtime process-output anchor must occur exactly once")
    return tcl.replace(vertex_anchor, process_put, 1)


def main() -> int:
    args = parse_args()
    variant = args.variant
    biases = tuple(args.biases)
    validate_biases(biases)
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
            raise FileNotFoundError(f"missing oracle input {name}: {source}")
        shutil.copy2(source, bundle / name)

    tcl_name = "runtime_element_avalanche_probe.tcl"
    stem = f"runtime_general_tri3_avalanche_probe_{variant}"
    deck_name = f"{stem}.cmd"
    write_ascii(
        bundle / tcl_name,
        oracle_tcl(
            (template_root / tcl_name).read_text(encoding="ascii"),
            biases,
        ),
    )
    write_ascii(
        bundle / deck_name,
        oracle_deck(
            (template_root / "runtime_element_avalanche_probe_default.cmd").read_text(
                encoding="ascii"
            ),
            variant,
            biases,
        ),
    )

    remote = f"{remote_root}/coarse7x3/{variant}"
    run([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote}"])
    for name in (*sources, tcl_name, deck_name):
        run([str(args.scp_bin), str(bundle / name), f"{args.ssh_target}:{remote}/"])
    run_name = f"run_{variant}.out"
    run(
        [
            str(args.ssh_bin),
            args.ssh_target,
            f"cd {remote} && sdevice {deck_name} > {run_name} 2>&1",
        ]
    )
    for name in (run_name, f"{stem}.plt", f"{stem}.tdr", f"{stem}_des.log"):
        run(
            [
                str(args.scp_bin),
                f"{args.ssh_target}:{remote}/{name}",
                str(fetched / name),
            ]
        )

    text = (fetched / run_name).read_text(encoding="ascii")
    observed = tuple(
        float(value)
        for value in re.findall(
            r"AVAL_PROBE_BEGIN bias_V=(-?\d+(?:\.\d+)?)",
            text,
        )
    )
    exact = observed == biases
    manifest = build_run_manifest(
        status="passed" if exact else "failed",
        experiment="pn2d_exact_high_bias_oracle_variant",
        variant=variant,
        exact_biases=biases,
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
        raise RuntimeError(f"exact lattice mismatch: expected {biases}, got {observed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
