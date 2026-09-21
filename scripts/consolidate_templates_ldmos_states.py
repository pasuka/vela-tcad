#!/usr/bin/env python3
"""Consolidate validated IdVg/IdVd states with a corrected BV capture."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from run_templates_ldmos_sentaurus_vm import records, sha256_file, utc_now, write_json


def copy_one(source: Path, destination: Path, provenance: list[dict[str, str]]) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, destination)
    provenance.append({
        "output": destination.name,
        "source": str(source),
        "source_sha256": sha256_file(source),
        "output_sha256": sha256_file(destination),
    })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-states", type=Path, required=True)
    parser.add_argument("--bv-states", type=Path, required=True)
    parser.add_argument("--state-decks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite consolidated evidence: {output}")
    raw, decks = output / "raw", output / "decks"
    raw.mkdir(parents=True)
    decks.mkdir()
    provenance: list[dict[str, str]] = []

    base_raw = args.base_states.resolve() / "raw"
    for pattern in ("state_idvg_*.tdr", "state_idvd_*.tdr"):
        for source in sorted(base_raw.glob(pattern)):
            copy_one(source, raw / source.name, provenance)
    for name in ("IdVg_n2_des.plt", "IdVd_Vg1_n4_des.plt", "IdVd_Vg2_n4_des.plt"):
        copy_one(base_raw / name, raw / name, provenance)

    bv_raw = args.bv_states.resolve() / "raw"
    for source in sorted(bv_raw.glob("state_bv_path_*.tdr")):
        copy_one(source, raw / source.name, provenance)
    copy_one(bv_raw / "n6_des.tdr", raw / "state_bv_path_final_des.tdr", provenance)
    copy_one(bv_raw / "n6_des.plt", raw / "n6_des.plt", provenance)

    for source in sorted(args.state_decks.resolve().iterdir()):
        if source.is_file():
            copy_one(source, decks / source.name, provenance)
    manifest = {
        "schema": "vela.templates_ldmos.consolidated_state_capture.v1",
        "classification": "derived_output_only_control_not_official_oracle",
        "generated_at": utc_now(),
        "policy": (
            "IdVg/IdVd come from the first output-only capture; BV comes from the corrected "
            "explicit-continuation-time capture; the derived final device state supplies the "
            "criterion-post fallback when the requested original-run time exceeds the derived run end."
        ),
        "provenance": provenance,
        "records": records(raw),
    }
    write_json(output / "consolidation_manifest.json", manifest)
    print(json.dumps({
        "output": str(output),
        "state_tdr_count": len(list(raw.glob("state_*.tdr"))),
        "record_count": len(manifest["records"]),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
