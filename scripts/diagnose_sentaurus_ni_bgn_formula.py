#!/usr/bin/env python3
"""Check Sentaurus PN2D ni_eff/BGN formula hypotheses against exported data."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


VT_300 = 8.617333262145e-5 * 300.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-par", type=Path, required=True)
    parser.add_argument("--materials-json", type=Path, required=True)
    parser.add_argument("--inferred-ni-csv", type=Path, required=True)
    parser.add_argument("--contact-inferred-ni-csv", type=Path)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--reference-ni-cm3", type=float, default=1.0e10)
    parser.add_argument("--thermal-voltage", type=float, default=VT_300)
    return parser.parse_args()


def parse_parameter(text: str, name: str, section: str | None = None) -> float:
    scope = text
    if section is not None:
        match = re.search(rf"^{re.escape(section)}\s*\{{(?P<body>.*?)^\}}", text, re.S | re.M)
        if not match:
            raise ValueError(f"missing section {section}")
        scope = match.group("body")
    pattern = rf"^\s*{re.escape(name)}\s*=\s*([-+0-9.eE]+)"
    match = re.search(pattern, scope, re.M)
    if not match:
        raise ValueError(f"missing parameter {name}")
    return float(match.group(1))


def parse_models(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    return {
        "bandgap_dEg0_old_slotboom_eV": parse_parameter(text, "dEg0(OldSlotboom)", "Bandgap"),
        "old_slotboom_Ebgn_eV": parse_parameter(text, "Ebgn", "OldSlotboom"),
        "old_slotboom_Nref_cm3": parse_parameter(text, "Nref", "OldSlotboom"),
        "old_slotboom_C": parse_parameter(text, "C", "OldSlotboom"),
        "bandgap_Bgn2Chi": parse_parameter(text, "Bgn2Chi", "Bandgap"),
        "bandgap_Eg0_eV": parse_parameter(text, "Eg0", "Bandgap"),
        "bandgap_alpha_eV_per_K": parse_parameter(text, "alpha", "Bandgap"),
        "bandgap_beta_K": parse_parameter(text, "beta", "Bandgap"),
    }


def read_material_ni(path: Path) -> float:
    data = json.loads(path.read_text(encoding="utf-8"))
    for item in data.get("materials", []):
        if item.get("name") == "Si":
            return float(item["ni"])
    raise ValueError(f"no Si material in {path}")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def first_numeric(row: dict[str, str], candidates: list[str]) -> float | None:
    for key in candidates:
        raw = row.get(key)
        if raw in (None, ""):
            continue
        try:
            value = float(raw)
        except ValueError:
            continue
        if math.isfinite(value):
            return value
    return None


def median(values: list[float]) -> float:
    clean = sorted(values)
    if not clean:
        raise ValueError("empty median input")
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def sentaurus_global_ni(path: Path) -> float:
    values = []
    for row in read_rows(path):
        bias = first_numeric(row, ["bias_V"])
        if bias is None or abs(bias) > 1.0e-9:
            continue
        value = first_numeric(row, ["median_ni_eff_cm3"])
        if value is not None:
            values.append(value)
    return median(values)


def sentaurus_contact_ni(path: Path | None) -> float | None:
    if path is None or not path.exists():
        return None
    values = []
    for row in read_rows(path):
        bias = first_numeric(row, ["bias_V"])
        location = row.get("location")
        value = first_numeric(row, ["median_ni_eff_cm3"])
        if bias is not None and abs(bias) <= 1.0e-9 and location == "contact" and value is not None:
            values.append(value)
    return median(values) if values else None


def doping_stats(path: Path) -> dict[str, float]:
    rows = read_rows(path)
    values = []
    for row in rows:
        value = first_numeric(
            row,
            [
                "component0",
                "DopingConcentration",
                "doping_cm3",
                "net_doping_cm3",
                "donors_cm3",
                "acceptors_cm3",
                "value",
            ],
        )
        if value is not None:
            values.append(abs(value))
    nonzero = [value for value in values if value > 0.0]
    if not nonzero:
        raise ValueError(f"no nonzero doping values in {path}")
    return {
        "points": float(len(values)),
        "abs_min_nonzero_cm3": min(nonzero),
        "abs_median_nonzero_cm3": median(nonzero),
        "abs_max_cm3": max(nonzero),
    }


def old_slotboom_term(params: dict[str, float], doping_cm3: float) -> float:
    x = math.log(doping_cm3 / params["old_slotboom_Nref_cm3"])
    return params["old_slotboom_Ebgn_eV"] * (
        x + math.sqrt(x * x + params["old_slotboom_C"])
    )


def ni_eff(ni_cm3: float, delta_eV: float, vt: float) -> float:
    return ni_cm3 * math.exp(delta_eV / (2.0 * vt))


def hypothesis_rows(
    params: dict[str, float],
    reference_ni_cm3: float,
    material_ni_cm3: float,
    target_ni_cm3: float,
    doping_cm3: float,
    vt: float,
) -> list[dict[str, Any]]:
    dEg0 = params["bandgap_dEg0_old_slotboom_eV"]
    term = old_slotboom_term(params, doping_cm3)
    raw_delta = dEg0 + term
    clamped_delta = max(raw_delta, 0.0)
    target_delta_from_reference = 2.0 * vt * math.log(target_ni_cm3 / reference_ni_cm3)
    offset_hack = target_delta_from_reference - params["old_slotboom_Ebgn_eV"] * math.sqrt(params["old_slotboom_C"])
    rows = [
        {
            "hypothesis": "vela_default_ni_plus_offset_clamped_delta",
            "fitted_to_target": False,
            "base_ni_cm3": reference_ni_cm3,
            "deltaEg_eV": clamped_delta,
            "predicted_ni_eff_cm3": ni_eff(reference_ni_cm3, clamped_delta, vt),
            "notes": "Current Vela old_slotboom semantics if material ni is not overridden.",
        },
        {
            "hypothesis": "models_par_material_ni_plus_offset_clamped_delta",
            "fitted_to_target": False,
            "base_ni_cm3": material_ni_cm3,
            "deltaEg_eV": clamped_delta,
            "predicted_ni_eff_cm3": ni_eff(material_ni_cm3, clamped_delta, vt),
            "notes": "Material Eg/dEg0 imported, but Vela still includes negative dEg0 in delta and clamps it.",
        },
        {
            "hypothesis": "models_par_material_ni_plus_raw_delta",
            "fitted_to_target": False,
            "base_ni_cm3": material_ni_cm3,
            "deltaEg_eV": raw_delta,
            "predicted_ni_eff_cm3": ni_eff(material_ni_cm3, raw_delta, vt),
            "notes": "No clamp, with dEg0 counted in both material ni and BGN delta.",
        },
        {
            "hypothesis": "sentaurus_split_dEg0_into_material_ni_then_positive_slotboom_term",
            "fitted_to_target": False,
            "base_ni_cm3": material_ni_cm3,
            "deltaEg_eV": term,
            "predicted_ni_eff_cm3": ni_eff(material_ni_cm3, term, vt),
            "notes": "dEg0 affects Bandgap/material ni; OldSlotboom positive term affects ni_eff.",
        },
        {
            "hypothesis": "legacy_offset_hack_with_reference_ni",
            "fitted_to_target": True,
            "base_ni_cm3": reference_ni_cm3,
            "deltaEg_eV": target_delta_from_reference,
            "predicted_ni_eff_cm3": target_ni_cm3,
            "notes": f"Equivalent one-parameter fit at Nref; offset_eV={offset_hack:.12g}.",
        },
    ]
    for row in rows:
        predicted = float(row["predicted_ni_eff_cm3"])
        row["target_ni_eff_cm3"] = target_ni_cm3
        row["ratio_predicted_over_target"] = predicted / target_ni_cm3
        row["log10_predicted_over_target"] = math.log10(predicted / target_ni_cm3)
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    params = parse_models(args.models_par)
    material_ni_cm3 = read_material_ni(args.materials_json)
    target_ni_cm3 = sentaurus_contact_ni(args.contact_inferred_ni_csv) or sentaurus_global_ni(args.inferred_ni_csv)
    doping = doping_stats(args.doping_csv)
    rows = hypothesis_rows(
        params,
        args.reference_ni_cm3,
        material_ni_cm3,
        target_ni_cm3,
        doping["abs_median_nonzero_cm3"],
        args.thermal_voltage,
    )
    write_csv(args.out_dir / "ni_bgn_formula_hypotheses.csv", rows)
    physical_rows = [row for row in rows if not bool(row.get("fitted_to_target"))]
    summary = {
        "models_par": str(args.models_par),
        "materials_json": str(args.materials_json),
        "inferred_ni_csv": str(args.inferred_ni_csv),
        "contact_inferred_ni_csv": str(args.contact_inferred_ni_csv) if args.contact_inferred_ni_csv else None,
        "doping_csv": str(args.doping_csv),
        "thermal_voltage_eV": args.thermal_voltage,
        "parameters": params,
        "doping_stats_cm3": doping,
        "reference_ni_cm3": args.reference_ni_cm3,
        "models_par_material_ni_cm3": material_ni_cm3,
        "sentaurus_target_ni_eff_cm3": target_ni_cm3,
        "old_slotboom_positive_term_eV_at_median_doping": old_slotboom_term(
            params,
            doping["abs_median_nonzero_cm3"],
        ),
        "rows": rows,
        "best_physical_hypothesis": min(
            physical_rows,
            key=lambda row: abs(float(row["log10_predicted_over_target"])),
        )["hypothesis"],
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    best = summary["best_physical_hypothesis"]
    lines = [
        "# Sentaurus ni_eff/BGN Formula Verification",
        "",
        f"- Target Sentaurus contact/global ni_eff: `{target_ni_cm3:.12g} cm^-3`.",
        f"- Exported nonzero doping median/max: `{doping['abs_median_nonzero_cm3']:.12g}` / `{doping['abs_max_cm3']:.12g} cm^-3`.",
        f"- models.par material ni: `{material_ni_cm3:.12g} cm^-3`.",
        f"- OldSlotboom positive term at median doping: `{summary['old_slotboom_positive_term_eV_at_median_doping']:.12g} eV`.",
        f"- Best non-fitted physical hypothesis: `{best}`.",
        "",
        "## Hypotheses",
        "",
        "| hypothesis | fitted | base ni cm^-3 | deltaEg eV | predicted ni_eff cm^-3 | ratio vs target |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| {hypothesis} | {fitted} | {base_ni_cm3:.6g} | {deltaEg_eV:.6g} | "
            "{predicted_ni_eff_cm3:.6g} | {ratio_predicted_over_target:.6g} |".format(
                fitted=str(row["fitted_to_target"]).lower(),
                **row,
            )
        )
    lines.append("")
    lines.append(
        "Interpretation: Sentaurus data are matched by applying `dEg0(OldSlotboom)` "
        "through the material Bandgap/intrinsic-density path, then applying only "
        "the positive OldSlotboom term to `ni_eff`."
    )
    (args.out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
