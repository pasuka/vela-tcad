"""Freeze, run, analyze, and verify the post-main M77 requalification."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import run_simplemos_m75_m77_poisson_charge_volume_decomposition as base


SCRIPT = Path(__file__).resolve()
BASE_SCRIPT = Path(base.__file__).resolve()
REPO = SCRIPT.parents[1]

base.PROFILES["M77"] = {
    "slug": "combined_poisson_charge_volume_post_main",
    "title": "post-main combined electron-hole-dopant Poisson charge volume",
    "prefix": "combined_poisson_volume_post_main",
    "contract": "simplemos_m77_combined_poisson_charge_volume_post_main_contract_v2.json",
    "freeze": "simplemos_m77_combined_poisson_charge_volume_post_main_contract_freeze_v2.json",
    "evidence": "simplemos_m77_combined_poisson_charge_volume_post_main_evidence.json",
}
base.SCRIPT = SCRIPT

_base_paths = base.paths


def post_main_paths(milestone: str) -> dict[str, Path]:
    result = _base_paths(milestone)
    result["doc"] = (
        REPO / "docs/validation" /
        "simplemos_m77_combined_poisson_charge_volume_post_main_2026-09-05.md"
    )
    result["artifact"] = (
        REPO / "docs/validation/reports/simplemos_m77_post_main/artifact.json"
    )
    return result


base.paths = post_main_paths
_base_source_paths = base.source_paths


def post_main_source_paths(contract: dict) -> list[Path]:
    paths = [
        path for path in _base_source_paths(contract)
        if path.name != "test_mos_mixed_material.cpp"
    ]
    paths.extend([
        BASE_SCRIPT,
        REPO / "CMakeLists.txt",
        REPO / "tests/test_poisson_charge_volume.cpp",
    ])
    return sorted(set(path.resolve() for path in paths))


base.source_paths = post_main_source_paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    milestone = "M77"
    item = base.paths(milestone)

    if args.freeze_contract:
        base.freeze_contract(milestone)
        print(json.dumps({
            "status": "frozen_before_execution",
            "milestone": milestone,
            "contract_sha256": base.sha256(item["contract"]),
        }))
        return

    contract = base.validate_contract(milestone)
    if args.verify:
        print(json.dumps(base.verify(milestone)["summary"], indent=2))
        return

    manifest = base.run_candidate(milestone, args.jobs) if args.run else \
        base.read_json(item["run_manifest"])
    if args.run and not args.analyze:
        print(json.dumps({
            "status": manifest["status"],
            "workflow_count": len(manifest["workflows"]),
        }))
        return
    if not args.analyze:
        parser.error("choose --freeze-contract, --run, --analyze, or --verify")

    result = base.analyze(milestone, contract)
    base.freeze_results(milestone, result)
    report = result[0]
    print(json.dumps({
        "status": report["status"],
        "classification": report["classification"],
        "summary": report["summary"],
        "acceptance": report["acceptance"],
    }, indent=2))


if __name__ == "__main__":
    main()
