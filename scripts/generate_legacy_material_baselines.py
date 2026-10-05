#!/usr/bin/env python3
"""Freeze inherited pre-WP13 scalars for historical calibrated TCAD decks.

These explicit files preserve old experiments; they are not default libraries.
Original partial files and sealed manifests remain unchanged. All numeric
concentrations and mobilities use unit_scaling (cm^-3 and cm^2/V/s), including
the historical Nc_m3/Nv_m3 field spelling.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LEGACY = {
    "Si": dict(name="Si", eps_r=11.7, ni=1e10, mun=1350.0, mup=480.0,
               bandgap_eV=1.12, electron_affinity_eV=4.05,
               Nc_m3=2.8e19, Nv_m3=1.04e19, temperature_K=300.0),
    "SiO2": dict(name="SiO2", eps_r=3.9, ni=0.0, mun=0.0, mup=0.0,
                 bandgap_eV=9.0, electron_affinity_eV=0.95, temperature_K=300.0),
}
SOURCES = [
    "reference_tcad/bvmethods_sentaurus2018/vela/materials_sentaurus2018.json",
    "reference_tcad/pn2d_sentaurus2018/source/pn2d_sentaurus2018_iv_materials.json",
    "reference_tcad/schottky_charon_sentaurus2018/vela/materials_mu_n_sentaurus2018.json",
    "reference_tcad/schottky_charon_sentaurus2018/vela/materials_mu_n_nc_eg_sentaurus2018.json",
    "reference_tcad/singledevice_sentaurus2018/vela/materials_sentaurus2018.json",
    "tests/fixtures/pn2d_bv_rejected_transition/materials.json",
]
BUILTIN_OUTPUT = "reference_tcad/schottky_charon_sentaurus2018/vela/materials_legacy_builtin_300k.json"


def complete(document: dict, source: str) -> dict:
    result = dict(document)
    result["_comment"] = ("WP13 frozen baseline: explicit pre-alignment defaults (9f9d2f5d) "
                          f"plus {source}; unit_scaling, 300 K. Generated; do not edit.")
    result["materials"] = [dict(LEGACY.get(entry["name"], {}), **entry)
                           for entry in document["materials"]]
    return result


def outputs() -> dict[Path, dict]:
    result = {REPO / BUILTIN_OUTPUT: complete({"materials": list(LEGACY.values())}, "built-ins")}
    for source in SOURCES:
        path = REPO / source
        document = complete(json.loads(path.read_text(encoding="utf-8-sig")), source)
        result[path.with_name(path.stem + "_legacy_300k.json")] = document
        if "singledevice" in source:
            result[path.parent.parent / "source" / (path.stem + "_legacy_300k.json")] = document
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    stale = []
    for path, document in outputs().items():
        text = json.dumps(document, indent=2) + "\n"
        if args.check:
            if not path.is_file() or path.read_text(encoding="utf-8") != text:
                stale.append(str(path.relative_to(REPO)))
        else:
            path.write_text(text, encoding="utf-8", newline="\n")
    if stale:
        print("Stale frozen material baselines: " + ", ".join(stale))
        return 1
    print("Frozen legacy material baselines are current" if args.check else "Wrote frozen legacy material baselines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
