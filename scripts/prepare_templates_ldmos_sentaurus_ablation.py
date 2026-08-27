#!/usr/bin/env python3
"""Generate fail-closed, one-factor Templates/LDMOS Id-Vg ablation decks."""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
from pathlib import Path
from typing import Callable, Sequence


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(text: str, pattern: str, replacement: str, label: str) -> str:
    result, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"{label}: expected exactly one match, found {count}")
    return result


def remove_hqp(text: str) -> str:
    return replace_once(text, r"\s*hQuantumPotential\(density\)", "", "remove hQP")


def remove_eqp(text: str) -> str:
    return replace_once(text, r"\s*eQuantumPotential\(density\)", "", "remove eQP")


def remove_ialmob(text: str) -> str:
    return replace_once(
        text, r"\s*Enormal\s*\(IALMob\(AutoOrientation\)\)", "", "remove IALMob")


def remove_high_field(text: str) -> str:
    return replace_once(text, r"\s*HighFieldSaturation", "", "remove high field")


def replace_polysi_with_barrier(text: str, barrier_V: float) -> str:
    return replace_once(
        text,
        r'Material\s*=\s*"PolySi"\(N\)',
        f"Barrier= {barrier_V:.17g}",
        "replace PolySi(N) gate",
    )


def equilibrium_capture(text: str) -> str:
    """Stop after the initial coupled equilibrium and save that state.

    This derivative changes no Physics statement.  It only removes the two
    bias ramps after the initial fixed point and adds one non-loadable Plot.
    The strict tail check prevents silently dropping commands outside Solve.
    """
    marker = "\tCoupled { Poisson Electron Hole }\n"
    if text.count(marker) != 1:
        raise ValueError(
            f"equilibrium capture: expected exactly one marker, found {text.count(marker)}")
    before, after = text.split(marker, 1)
    if not after.rstrip().endswith("}"):
        raise ValueError("equilibrium capture: source has content after the Solve block")
    return (
        before + marker +
        '\tPlot(-Loadable FilePrefix="G4_equilibrium_state")\n'
        "}\n"
    )


def prepare(source: Path, output: Path, gate_barrier_V: float) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=True)
    original = source.read_text(encoding="utf-8")
    transforms: list[tuple[str, str | None, Callable[[str], str]]] = [
        ("G0-original", None, lambda value: value),
        ("G1-no-hQP", "G0-original", remove_hqp),
        ("G2-no-eQP", "G1-no-hQP", remove_eqp),
        ("G3-no-IALMob", "G2-no-eQP", remove_ialmob),
        ("G4-no-highfield", "G3-no-IALMob", remove_high_field),
    ]
    texts: dict[str, str] = {}
    records: list[dict[str, object]] = []
    for name, parent, transform in transforms:
        parent_text = original if parent is None else texts[parent]
        text = transform(parent_text)
        texts[name] = text
        deck = output / f"{name}.cmd"
        deck.write_text(text, encoding="utf-8")
        diff = output / f"{name}.diff"
        diff.write_text("".join(difflib.unified_diff(
            parent_text.splitlines(keepends=True), text.splitlines(keepends=True),
            fromfile=parent or source.name, tofile=name)), encoding="utf-8")
        records.append({
            "id": name, "parent": parent, "deck": deck.name,
            "deck_sha256": sha256(deck), "diff": diff.name,
            "diff_sha256": sha256(diff),
        })

    equilibrium_name = "G4-equilibrium-capture"
    equilibrium_text = equilibrium_capture(texts["G4-no-highfield"])
    equilibrium = output / f"{equilibrium_name}.cmd"
    equilibrium.write_text(equilibrium_text, encoding="utf-8")
    equilibrium_diff = output / f"{equilibrium_name}.diff"
    equilibrium_diff.write_text("".join(difflib.unified_diff(
        texts["G4-no-highfield"].splitlines(keepends=True),
        equilibrium_text.splitlines(keepends=True),
        fromfile="G4-no-highfield", tofile=equilibrium_name)), encoding="utf-8")
    records.append({
        "id": equilibrium_name,
        "parent": "G4-no-highfield",
        "deck": equilibrium.name,
        "deck_sha256": sha256(equilibrium),
        "diff": equilibrium_diff.name,
        "diff_sha256": sha256(equilibrium_diff),
        "classification": "derived_output_and_path_control_not_curve_oracle",
        "only_intended_changes": [
            "stop after the initial coupled equilibrium",
            "write one non-loadable equilibrium state",
        ],
    })

    control_name = "G-contact-poly-barrier-control"
    control_text = replace_polysi_with_barrier(original, gate_barrier_V)
    control = output / f"{control_name}.cmd"
    control.write_text(control_text, encoding="utf-8")
    control_diff = output / f"{control_name}.diff"
    control_diff.write_text("".join(difflib.unified_diff(
        original.splitlines(keepends=True), control_text.splitlines(keepends=True),
        fromfile="G0-original", tofile=control_name)), encoding="utf-8")
    records.append({
        "id": control_name, "parent": "G0-original", "deck": control.name,
        "deck_sha256": sha256(control), "diff": control_diff.name,
        "diff_sha256": sha256(control_diff),
    })

    manifest: dict[str, object] = {
        "schema": "vela.templates_ldmos.sentaurus_idvg_ablation.v1",
        "source": str(source.resolve()),
        "source_sha256": sha256(source),
        "gate_barrier_V": gate_barrier_V,
        "single_factor_chain": records,
        "execution_status": "prepared_not_run",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gate-barrier-V", type=float, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(
        args.source.resolve(), args.output_dir.resolve(), args.gate_barrier_V), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
