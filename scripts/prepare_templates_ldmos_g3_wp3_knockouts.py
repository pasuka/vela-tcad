#!/usr/bin/env python3
"""Prepare the frozen WP3 T4 G3 Sentaurus knockout matrix without running it.

The generator is intentionally offline.  It derives three one-factor decks
from an already classical G3 endpoint-capture deck and emits a reviewable VM
execution plan.  It never contacts the Sentaurus VM.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any, Sequence


VARIANTS: dict[str, dict[str, str]] = {
    "k1_hfs_off": {
        "label": "K1",
        "hypothesis": "H1",
        "sentaurus_delta": "remove Mobility(HighFieldSaturation)",
        "vela_delta": "solver.mobility.model=constant; no HFS fallback",
    },
    "k2_boltzmann": {
        "label": "K2",
        "hypothesis": "H2_degeneracy_support",
        "sentaurus_delta": "Fermi statistics to default Boltzmann statistics",
        "vela_delta": "carrier_statistics=boltzmann with matching contact reconstruction",
    },
    "k3_bgn_off": {
        "label": "K3",
        "hypothesis": "H2_band_edge_placement",
        "sentaurus_delta": "remove EffectiveIntrinsicDensity(OldSlotboom)",
        "vela_delta": "bandgap_narrowing disabled with matching contact reconstruction",
    },
}

BIAS_POINTS_V = (1.0, 1.166666666666667)
FORBIDDEN_SOURCE_TOKENS = ("IALMob", "eQuantumPotential", "hQuantumPotential", "Predictor")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def replace_exactly_once(text: str, pattern: str, replacement: str, label: str) -> str:
    result, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise ValueError(f"{label}: expected exactly one match, found {count}")
    return result


def validate_frozen_g3_source(text: str) -> None:
    forbidden = [token for token in FORBIDDEN_SOURCE_TOKENS if token.lower() in text.lower()]
    if forbidden:
        raise ValueError(f"source is not the frozen classical G3 contract: {forbidden}")
    required_patterns = {
        "Fermi": r"Physics\s*\{\s*Fermi\s*\}",
        "HighFieldSaturation": r"Mobility\s*\(\s*HighFieldSaturation\s*\)",
        "OldSlotboom": r"EffectiveIntrinsicDensity\s*\(\s*OldSlotboom\s*\)",
        "SRH": r"\bSRH\s*\(",
        "Auger": r"\bAuger\b",
        "drain goal": r"Goal\s*\{[^}]*Name\s*=\s*\"drain\"[^}]*Voltage\s*=\s*0\.1\b",
    }
    missing = [label for label, pattern in required_patterns.items() if not re.search(pattern, text, re.I | re.S)]
    if missing:
        raise ValueError(f"source is missing frozen G3 clauses: {missing}")
    if len(re.findall(r"\bHighFieldSaturation\b", text, re.I)) != 1:
        raise ValueError("source must contain exactly one HighFieldSaturation selector")
    if len(re.findall(r"EffectiveIntrinsicDensity\s*\(\s*OldSlotboom\s*\)", text, re.I)) != 1:
        raise ValueError("source must contain exactly one OldSlotboom selector")


def apply_knockout(text: str, variant: str) -> str:
    if variant == "k1_hfs_off":
        return replace_exactly_once(
            text, r"\s*Mobility\s*\(\s*HighFieldSaturation\s*\)",
            "\n\t* K1: HFS intentionally removed; constant low-field mobility",
            "K1 HFS removal",
        )
    if variant == "k2_boltzmann":
        return replace_exactly_once(
            text, r"Physics\s*\{\s*Fermi\s*\}",
            "* K2: Fermi intentionally removed; Sentaurus default Boltzmann statistics",
            "K2 Fermi-to-Boltzmann replacement",
        )
    if variant == "k3_bgn_off":
        return replace_exactly_once(
            text, r"\s*EffectiveIntrinsicDensity\s*\(\s*OldSlotboom\s*\)", "",
            "K3 OldSlotboom removal",
        )
    raise ValueError(f"unknown knockout variant: {variant}")


def replace_file_outputs(text: str, variant: str) -> str:
    replacements = {
        "Output": f'Output=     "{variant}.log"',
        "Current": f'Current=    "{variant}.plt"',
        "Plot": f'Plot=       "{variant}.tdr"',
    }
    result = text
    for key, replacement in replacements.items():
        result = replace_exactly_once(
            result, rf"^\s*{key}\s*=\s*\"[^\"]+\"\s*$", f"\t{replacement}",
            f"isolate {key} output",
        )
    return result


def endpoint_solve_block(variant: str) -> str:
    return f'''Solve {{
\tCoupled(Iterations=100 LineSearchDamping=1e-4){{ Poisson }}
\tCoupled {{ Poisson Electron Hole }}

\tQuasistationary(
\t\tInitialStep=1e-2 Increment=1.35
\t\tMinStep=1e-5 MaxStep=0.2
\t\tGoal {{ Name="drain" Voltage=0.1 }}
\t){{ Coupled {{ Poisson Electron Hole }} }}

\tNewCurrentPrefix="{variant}_vg1p0_"
\tQuasistationary(
\t\tDoZero InitialStep=0.01 Increment=1.35
\t\tMinStep=1e-5 MaxStep=0.1
\t\tGoal {{ Name="gate" Voltage=1.0 }}
\t){{ Coupled {{ Poisson Electron Hole }}
\t\tCurrentPlot(Time=(1))
\t}}
\tPlot(-Loadable FilePrefix="{variant}_vg1p0")

\tNewCurrentPrefix="{variant}_vg1p166667_"
\tQuasistationary(
\t\tInitialStep=0.01 Increment=1.35
\t\tMinStep=1e-5 MaxStep=0.1
\t\tGoal {{ Name="gate" Voltage=1.166666666666667 }}
\t){{ Coupled {{ Poisson Electron Hole }}
\t\tCurrentPlot(Time=(1))
\t}}
\tPlot(-Loadable FilePrefix="{variant}_vg1p166667")
}}
'''


def replace_solve_block(text: str, variant: str) -> str:
    match = re.search(r"(?ms)^Solve\s*\{.*\}\s*$", text)
    if not match:
        raise ValueError("source must end with exactly one Solve block")
    if re.search(r"(?m)^Solve\s*\{", text[:match.start()]):
        raise ValueError("source contains multiple Solve blocks")
    return text[:match.start()] + endpoint_solve_block(variant)


def assert_isolated_variant(text: str, variant: str) -> None:
    lowered = text.lower()
    if "ialmob" in lowered or "predictor" in lowered:
        raise ValueError(f"{variant}: forbidden IALMob/predictor selector")
    counts = {
        "Fermi": len(re.findall(r"Physics\s*\{\s*Fermi\s*\}", text, re.I)),
        "HFS": len(re.findall(r"\bHighFieldSaturation\b", text, re.I)),
        "BGN": len(re.findall(r"EffectiveIntrinsicDensity\s*\(\s*OldSlotboom\s*\)", text, re.I)),
    }
    expected = {
        "k1_hfs_off": {"Fermi": 1, "HFS": 0, "BGN": 1},
        "k2_boltzmann": {"Fermi": 0, "HFS": 1, "BGN": 1},
        "k3_bgn_off": {"Fermi": 1, "HFS": 1, "BGN": 0},
    }[variant]
    if counts != expected:
        raise ValueError(f"{variant}: physics isolation failed: {counts} != {expected}")
    if len(re.findall(r"Plot\s*\(\s*-Loadable\s+FilePrefix=", text, re.I)) != 2:
        raise ValueError(f"{variant}: exactly two endpoint state plots are required")
    if "Intervals=30" in text or "Voltage=5.0" in text.replace(" ", ""):
        raise ValueError(f"{variant}: accidental 31-point/full IdVg sweep")


def prepare(source_deck: Path, parameter_file: Path, grid_file: Path, output: Path) -> dict[str, Any]:
    source_deck = source_deck.resolve()
    parameter_file = parameter_file.resolve()
    grid_file = grid_file.resolve()
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty T4 preparation: {output}")
    source_text = source_deck.read_text(encoding="utf-8")
    validate_frozen_g3_source(source_text)
    if not parameter_file.is_file():
        raise FileNotFoundError(parameter_file)
    if not grid_file.is_file():
        raise FileNotFoundError(grid_file)
    output.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    for variant, metadata in VARIANTS.items():
        variant_dir = output / variant
        variant_dir.mkdir()
        transformed = apply_knockout(source_text, variant)
        transformed = replace_file_outputs(transformed, variant)
        transformed = replace_solve_block(transformed, variant)
        assert_isolated_variant(transformed, variant)
        deck = variant_dir / "IdVg.cmd"
        deck.write_text(transformed, encoding="utf-8")
        diff = variant_dir / "IdVg.diff"
        diff.write_text("".join(difflib.unified_diff(
            source_text.splitlines(keepends=True), transformed.splitlines(keepends=True),
            fromfile=source_deck.name, tofile=f"{variant}/IdVg.cmd",
        )), encoding="utf-8")
        replay_contract = {
            "schema": "vela.templates_ldmos.g3_wp3_knockout_replay_contract.v1",
            "variant": variant,
            "bias_points_V": list(BIAS_POINTS_V),
            "sentaurus_physics_delta": metadata["sentaurus_delta"],
            "vela_physics_delta": metadata["vela_delta"],
            "state_pairing": "variant_physics_on_same_variant_converged_state_only",
            "fixed_constraints": {
                "drain_voltage_V": 0.1,
                "ialmob_enabled": False,
                "vela_predictor_enabled": False,
                "idvg_31_point_run": False,
            },
            "required_endpoint_artifacts": [
                f"{variant}_vg1p0*.tdr", f"{variant}_vg1p166667*.tdr",
            ],
        }
        contract_path = variant_dir / "replay_contract.json"
        write_json(contract_path, replay_contract)
        sealed_parameter = variant_dir / "sdevice.par"
        shutil.copyfile(parameter_file, sealed_parameter)
        state_manifest = {
            "schema": "vela.templates_ldmos.conditional_state_deck.v1",
            "classification": f"WP3 T4 {metadata['label']} isolated endpoint state regeneration",
            "variant": variant,
            "physics_delta_from_frozen_g3": metadata["sentaurus_delta"],
            "bias_points_V": list(BIAS_POINTS_V),
            "deck": "IdVg.cmd",
            "deck_sha256": sha256_file(deck),
            "parameter": "sdevice.par",
            "parameter_sha256": sha256_file(sealed_parameter),
            "expected_grid_name": "n1_fps.tdr",
            "expected_grid_sha256": sha256_file(grid_file),
            "replay_contract": "replay_contract.json",
            "replay_contract_sha256": sha256_file(contract_path),
            "execution_status": "prepared_not_run",
        }
        state_manifest_path = variant_dir / "state_deck_manifest.json"
        write_json(state_manifest_path, state_manifest)
        records.append({
            "variant": variant,
            **metadata,
            "deck": f"{variant}/IdVg.cmd",
            "deck_sha256": sha256_file(deck),
            "diff": f"{variant}/IdVg.diff",
            "diff_sha256": sha256_file(diff),
            "replay_contract": f"{variant}/replay_contract.json",
            "replay_contract_sha256": sha256_file(contract_path),
            "parameter_copy": f"{variant}/sdevice.par",
            "parameter_copy_sha256": sha256_file(sealed_parameter),
            "state_deck_manifest": f"{variant}/state_deck_manifest.json",
            "state_deck_manifest_sha256": sha256_file(state_manifest_path),
            "endpoint_prefixes": [f"{variant}_vg1p0", f"{variant}_vg1p166667"],
        })

    execution_plan = {
        "schema": "vela.templates_ldmos.g3_wp3_knockout_execution_plan.v1",
        "status": "prepared_not_submitted",
        "sentaurus_release": "T-2022.03-SP2",
        "ssh_target": "sentaurus",
        "remote_root": "~/sentaurus_runs/vela_oracle_2022/templates_ldmos_g3_wp3_t4",
        "parent_grid_requirement": "copy the sealed phase01 n1_fps.tdr into each independent variant directory",
        "expected_parent_grid_sha256": sha256_file(grid_file),
        "per_variant_sequence": [
            "create a unique remote directory and refuse overwrite",
            "copy the same sealed n1_fps.tdr; verify its SHA-256",
            "upload only IdVg.cmd and the sealed sdevice.par",
            "run sdevice IdVg.cmd once; this converges and plots only the two registered endpoints",
            "archive log, PLT, both endpoint TDRs, timing, exit code and sha256sum output",
            "download to ignored reference_staging and record local SHA-256",
            "import each endpoint TDR independently with sentaurus_import",
            "replay only with the matching variant physics contract",
            "write one analysis/summary.json per variant per bias",
        ],
        "host_runner": {
            "script": "scripts/run_templates_ldmos_state_capture.py",
            "arguments": [
                "--run-dir <sealed-phase01-run-dir>",
                "--state-decks <prepared-root>/<variant>",
                "--stages idvg",
                "--output-name g3_wp3_t4_<variant>_<unique-stamp>"
            ],
            "sealing": "runner hashes the state-deck manifest and every downloaded raw artifact, including endpoint TDRs"
        },
        "import_command": "build-release/sentaurus_import.exe --tdr <endpoint.tdr> --export-dir <endpoint-import-dir> --coordinate-unit um --compensated-doping-policy reported",
        "prohibited": ["IALMob", "Vela predictor", "31-point IdVg", "production-default changes"],
        "expected_runs": 3,
        "expected_endpoint_states": 6,
    }
    write_json(output / "execution_plan.json", execution_plan)
    manifest = {
        "schema": "vela.templates_ldmos.g3_wp3_knockout_manifest.v1",
        "frozen_plan_commit": "01f20ace88dc2b68a7bc3788534244dadf0d18f2",
        "status": "prepared_not_run",
        "source_deck": str(source_deck),
        "source_deck_sha256": sha256_file(source_deck),
        "parameter_file": str(parameter_file),
        "parameter_file_sha256": sha256_file(parameter_file),
        "grid_file": str(grid_file),
        "grid_file_sha256": sha256_file(grid_file),
        "bias_points_V": list(BIAS_POINTS_V),
        "variants": records,
        "execution_plan_sha256": sha256_file(output / "execution_plan.json"),
    }
    write_json(output / "manifest.json", manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-deck", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--grid-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(
        args.source_deck, args.parameter_file, args.grid_file, args.output_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
