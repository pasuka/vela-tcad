#!/usr/bin/env python3
"""Generate the fail-closed Templates/LDMOS Id-Vd D1--D5 chain."""

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


def replace_exact(text: str, old: str, new: str, count: int, label: str) -> str:
    observed = text.count(old)
    if observed != count:
        raise ValueError(f"{label}: expected {count} matches, found {observed}")
    return text.replace(old, new)


def remove_isothermal_coupling(text: str) -> str:
    result, count = re.subn(
        r"\nThermode\s*\{.*?\n\}\s*\n",
        "\n",
        text,
        count=1,
        flags=re.DOTALL,
    )
    if count != 1:
        raise ValueError(f"remove Thermode: expected 1 match, found {count}")
    return replace_exact(
        result,
        "Coupled { Poisson Electron Hole Temperature }",
        "Coupled { Poisson Electron Hole }",
        4,
        "remove Temperature equation",
    )


def remove_hrec_velocity(text: str) -> str:
    result, count = re.subn(r"\s+hRecVelocity\s*=\s*1\.93E6", "", text)
    if count != 2:
        raise ValueError(f"remove hRecVelocity: expected 2 matches, found {count}")
    return result


def remove_once(text: str, pattern: str, label: str) -> str:
    result, count = re.subn(pattern, "", text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"{label}: expected 1 match, found {count}")
    return result


def prepare(source: Path, parameter_file: Path, output: Path) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=False)
    original = source.read_text(encoding="utf-8")
    transforms: list[tuple[str, str | None, Callable[[str], str], str]] = [
        ("D0-original", None, lambda value: value, "official full-physics baseline"),
        ("D1-isothermal", "D0-original", remove_isothermal_coupling,
         "remove Thermode and Temperature equation coupling"),
        ("D2-no-hRecVelocity", "D1-isothermal", remove_hrec_velocity,
         "remove source/drain hRecVelocity"),
        ("D3-no-hQP", "D2-no-hRecVelocity",
         lambda value: remove_once(value, r"\s*hQuantumPotential\(density\)", "remove hQP"),
         "remove hQuantumPotential"),
        ("D4-classical", "D3-no-hQP",
         lambda value: remove_once(value, r"\s*eQuantumPotential\(density\)", "remove eQP"),
         "remove eQuantumPotential"),
        ("D5-no-IALMob", "D4-classical",
         lambda value: remove_once(
             value, r"\s*Enormal\s*\(IALMob\(AutoOrientation\)\)", "remove IALMob"),
         "remove IALMob; retain HighFieldSaturation"),
    ]
    texts: dict[str, str] = {}
    records: list[dict[str, object]] = []
    for name, parent, transform, factor in transforms:
        parent_text = original if parent is None else texts[parent]
        rendered = transform(parent_text)
        texts[name] = rendered
        variant = output / name
        variant.mkdir()
        deck = variant / "IdVd.cmd"
        deck.write_text(rendered, encoding="utf-8")
        parameter_destination = variant / "sdevice.par"
        parameter_destination.write_bytes(parameter_file.read_bytes())
        diff = variant / "parent.diff"
        diff.write_text("".join(difflib.unified_diff(
            parent_text.splitlines(keepends=True),
            rendered.splitlines(keepends=True),
            fromfile=parent or source.name,
            tofile=name,
        )), encoding="utf-8")
        records.append({
            "id": name,
            "parent": parent,
            "factor": factor,
            "deck": f"{name}/IdVd.cmd",
            "deck_sha256": sha256(deck),
            "parameter_sha256": sha256(parameter_destination),
            "diff": f"{name}/parent.diff",
            "diff_sha256": sha256(diff),
        })
    manifest: dict[str, object] = {
        "schema": "vela.templates_ldmos.sentaurus_idvd_ablation.v1",
        "source": str(source.resolve()),
        "source_sha256": sha256(source),
        "parameter_source": str(parameter_file.resolve()),
        "parameter_sha256": sha256(parameter_file),
        "single_factor_chain": records,
        "execution_status": "prepared_not_run",
        "common_invariants": [
            "same final SProcess n1_fps.tdr",
            "same Fermi/OldSlotboom/SRH/Auger and mobility parameters",
            "same gate/drain targets, CurrentPlot points and Math settings",
            "same two-dimensional width and electrical contacts",
        ],
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(
        args.source.resolve(), args.parameter_file.resolve(), args.output_dir.resolve()),
        indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
