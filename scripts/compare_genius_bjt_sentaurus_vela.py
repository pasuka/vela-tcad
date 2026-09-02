#!/usr/bin/env python3
"""Compare the Genius NPN BJT Sentaurus oracle with Vela DC sweeps."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


TERMINALS = ("collector", "base", "emitter")
TARGET_VCE = tuple(index / 10.0 for index in range(31))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def select_sentaurus(rows: list[dict[str, str]], target: float) -> dict[str, str]:
    row = min(rows, key=lambda item: abs(float(item["VCE_V"]) - target))
    if abs(float(row["VCE_V"]) - target) > 1.0e-8:
        raise ValueError(f"Sentaurus curve has no point at VCE={target:g} V")
    return row


def select_vela_terminals(
    rows: list[dict[str, str]], target: float
) -> tuple[dict[str, dict[str, str]], float]:
    available = sorted({float(row["bias_V"]) for row in rows})
    bias = min(available, key=lambda value: abs(value - target))
    error = abs(bias - target)
    if error > 1.0e-8:
        raise ValueError(f"Vela curve has no point at VCE={target:g} V")
    selected = {
        row["contact"].lower(): row
        for row in rows
        if abs(float(row["bias_V"]) - bias) <= 1.0e-12
    }
    missing = set(TERMINALS) - set(selected)
    if missing:
        raise ValueError(f"Vela VCE={target:g} V is missing contacts: {sorted(missing)}")
    return selected, error


def beta_abs(ic: float, ib: float) -> float:
    return abs(ic / ib) if ib else math.inf


def magnitude_ratio(actual: float, reference: float) -> float:
    return abs(actual) / abs(reference) if reference else math.inf


def log_magnitude_error(actual: float, reference: float) -> float:
    if not actual or not reference:
        return math.inf
    return abs(math.log10(abs(actual) / abs(reference)))


def finite_stats(values: list[float]) -> dict[str, float | None]:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return {"median": None, "maximum": None}
    middle = len(finite) // 2
    median = (
        finite[middle]
        if len(finite) % 2
        else 0.5 * (finite[middle - 1] + finite[middle])
    )
    return {"median": median, "maximum": finite[-1]}


def compare_model(
    model: str,
    sentaurus_path: Path,
    vela_balance_path: Path,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    sentaurus = read_csv(sentaurus_path)
    vela = read_csv(vela_balance_path)
    output: list[dict[str, object]] = []
    max_voltage_error = 0.0
    all_converged = True

    for target in TARGET_VCE:
        srow = select_sentaurus(sentaurus, target)
        terminals, voltage_error = select_vela_terminals(vela, target)
        max_voltage_error = max(max_voltage_error, voltage_error)
        currents = {
            terminal: float(terminals[terminal]["current_total_A_per_um"])
            for terminal in TERMINALS
        }
        converged = all(terminals[name]["converged"] == "1" for name in TERMINALS)
        all_converged = all_converged and converged
        kcl = sum(currents.values())
        scale = max(abs(value) for value in currents.values())
        relative_kcl = abs(kcl) / max(scale, 1.0e-30)
        sic = float(srow["Ic_A_per_um"])
        sib = float(srow["Ib_A_per_um"])
        sie = float(srow["Ie_A_per_um"])
        vic = currents["collector"]
        vib = currents["base"]
        vie = currents["emitter"]
        sbeta = beta_abs(sic, sib)
        vbeta = beta_abs(vic, vib)
        output.append(
            {
                "model": model,
                "VBE_V": float(srow["VBE_V"]),
                "VCE_V": target,
                "sentaurus_Ic_A_per_um": sic,
                "vela_Ic_A_per_um": vic,
                "Ic_magnitude_ratio_vela_over_sentaurus": magnitude_ratio(vic, sic),
                "Ic_absolute_log10_error": log_magnitude_error(vic, sic),
                "sentaurus_Ib_A_per_um": sib,
                "vela_Ib_A_per_um": vib,
                "Ib_magnitude_ratio_vela_over_sentaurus": magnitude_ratio(vib, sib),
                "Ib_absolute_log10_error": log_magnitude_error(vib, sib),
                "sentaurus_Ie_A_per_um": sie,
                "vela_Ie_A_per_um": vie,
                "sentaurus_beta_abs": sbeta,
                "vela_beta_abs": vbeta,
                "beta_magnitude_ratio_vela_over_sentaurus": magnitude_ratio(vbeta, sbeta),
                "sentaurus_kcl_abs_A_per_um": float(srow["kcl_abs_A_per_um"]),
                "sentaurus_kcl_relative": float(srow["kcl_relative"]),
                "vela_kcl_abs_A_per_um": abs(kcl),
                "vela_kcl_relative": relative_kcl,
                "vela_converged": int(converged),
            }
        )

    active = [row for row in output if float(row["VCE_V"]) >= 0.5]
    final = output[-1]
    max_kcl_abs = max(float(row["vela_kcl_abs_A_per_um"]) for row in output)
    max_kcl_relative = max(float(row["vela_kcl_relative"]) for row in output)
    operational_pass = (
        len(output) == 31
        and max_voltage_error <= 1.0e-8
        and all_converged
        and (max_kcl_abs <= 1.0e-18 or max_kcl_relative <= 1.0e-6)
    )
    summary: dict[str, object] = {
        "model": model,
        "sentaurus_input": str(sentaurus_path),
        "vela_terminal_balance_input": str(vela_balance_path),
        "compared_point_count": len(output),
        "maximum_voltage_selection_error_V": max_voltage_error,
        "all_vela_points_converged": all_converged,
        "maximum_vela_absolute_kcl_residual_A_per_um": max_kcl_abs,
        "maximum_vela_relative_kcl_residual": max_kcl_relative,
        "active_region_VCE_range_V": [0.5, 3.0],
        "active_region_Ic_absolute_log10_error": finite_stats(
            [float(row["Ic_absolute_log10_error"]) for row in active]
        ),
        "active_region_Ib_absolute_log10_error": finite_stats(
            [float(row["Ib_absolute_log10_error"]) for row in active]
        ),
        "at_VCE_3V": {
            "sentaurus_Ic_A_per_um": final["sentaurus_Ic_A_per_um"],
            "vela_Ic_A_per_um": final["vela_Ic_A_per_um"],
            "Ic_magnitude_ratio_vela_over_sentaurus": final[
                "Ic_magnitude_ratio_vela_over_sentaurus"
            ],
            "sentaurus_Ib_A_per_um": final["sentaurus_Ib_A_per_um"],
            "vela_Ib_A_per_um": final["vela_Ib_A_per_um"],
            "Ib_magnitude_ratio_vela_over_sentaurus": final[
                "Ib_magnitude_ratio_vela_over_sentaurus"
            ],
            "sentaurus_beta_abs": final["sentaurus_beta_abs"],
            "vela_beta_abs": final["vela_beta_abs"],
        },
        "numerical_parity_gate": {
            "status": "characterization_only",
            "reason": "WP3-WP5 records the first common-input comparison; no numerical parity tolerance was pre-registered.",
        },
        "operational_pass": operational_pass,
    }
    return output, summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        handle.write(text)


def write_markdown(path: Path, summary: dict[str, object]) -> None:
    lines = [
        "# Genius NPN BJT Sentaurus/Vela comparison",
        "",
        "All 31 requested collector biases use directly reported collector, base, and emitter currents.",
        "Operational pass means exact bias selection, solver convergence, and terminal KCL passed.",
        "Numerical parity is intentionally characterization-only because WP0-WP2 did not pre-register a Vela error tolerance.",
        "",
        "| Model | Operational | Sentaurus Ic @ 3 V (A/um) | Vela Ic @ 3 V (A/um) | Ic ratio | Sentaurus beta | Vela beta |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for model in ("M0", "M1"):
        item = summary[model]
        point = item["at_VCE_3V"]
        lines.append(
            f"| {model} | {item['operational_pass']} | "
            f"{point['sentaurus_Ic_A_per_um']:.9e} | {point['vela_Ic_A_per_um']:.9e} | "
            f"{point['Ic_magnitude_ratio_vela_over_sentaurus']:.6g} | "
            f"{point['sentaurus_beta_abs']:.6g} | {point['vela_beta_abs']:.6g} |"
        )
    lines.extend(["", f"Overall operational pass: **{summary['operational_pass']}**", ""])
    write_text_lf(path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--vela-run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    summaries: dict[str, object] = {
        "schema_version": 1,
        "comparison_scope": "WP3-WP5 first common-input characterization",
    }
    for model, stem in (("M0", "m0"), ("M1", "m1")):
        rows, model_summary = compare_model(
            model,
            args.reference_root / "reference_curves" / f"bjt_{stem}_output.csv",
            args.vela_run_root / f"{stem}_collector_terminal_balance.csv",
        )
        write_csv(args.output_dir / f"{stem}_sentaurus_vela.csv", rows)
        summaries[model] = model_summary
    summaries["operational_pass"] = bool(
        summaries["M0"]["operational_pass"] and summaries["M1"]["operational_pass"]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_text_lf(
        args.output_dir / "comparison_summary.json",
        json.dumps(summaries, indent=2, allow_nan=False) + "\n",
    )
    write_markdown(args.output_dir / "comparison_summary.md", summaries)
    print(json.dumps(summaries, indent=2, allow_nan=False))
    return 0 if summaries["operational_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
