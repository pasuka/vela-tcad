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


def evaluate_numerical_parity(
    active_rows: list[dict[str, object]],
    contract: dict[str, object],
) -> dict[str, object]:
    status = str(contract.get("status", "characterization_only"))
    reason = str(contract.get("reason", ""))
    if status == "characterization_only":
        return {"status": status, "reason": reason, "pass": None}
    if status != "asserted":
        raise ValueError(f"unsupported numerical parity status: {status}")
    if not active_rows:
        raise ValueError("asserted numerical parity gate has no active-region rows")

    metric_columns = {
        "Ic": "Ic_absolute_log10_error",
        "Ib": "Ib_absolute_log10_error",
        "beta": "beta_absolute_log10_error",
    }
    thresholds = {
        name: float(contract[f"maximum_{name}_absolute_log10_error"])
        for name in metric_columns
    }
    if any(not math.isfinite(value) or value < 0.0 for value in thresholds.values()):
        raise ValueError("numerical parity thresholds must be finite and non-negative")
    observed = {
        name: max(float(row[column]) for row in active_rows)
        for name, column in metric_columns.items()
    }
    passed = all(observed[name] <= thresholds[name] for name in metric_columns)
    return {
        "status": "asserted",
        "reason": reason,
        "thresholds_maximum_absolute_log10_error": thresholds,
        "observed_maximum_absolute_log10_error": observed,
        "pass": passed,
    }


def compare_model(
    model: str,
    sentaurus_path: Path,
    vela_balance_path: Path,
    active_region: tuple[float, float],
    parity_contract: dict[str, object],
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
                "beta_absolute_log10_error": log_magnitude_error(vbeta, sbeta),
                "sentaurus_kcl_abs_A_per_um": float(srow["kcl_abs_A_per_um"]),
                "sentaurus_kcl_relative": float(srow["kcl_relative"]),
                "vela_kcl_abs_A_per_um": abs(kcl),
                "vela_kcl_relative": relative_kcl,
                "vela_converged": int(converged),
            }
        )

    active = [
        row
        for row in output
        if active_region[0] <= float(row["VCE_V"]) <= active_region[1]
    ]
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
        "active_region_VCE_range_V": list(active_region),
        "active_region_Ic_absolute_log10_error": finite_stats(
            [float(row["Ic_absolute_log10_error"]) for row in active]
        ),
        "active_region_Ib_absolute_log10_error": finite_stats(
            [float(row["Ib_absolute_log10_error"]) for row in active]
        ),
        "active_region_beta_absolute_log10_error": finite_stats(
            [float(row["beta_absolute_log10_error"]) for row in active]
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
        "numerical_parity_gate": evaluate_numerical_parity(active, parity_contract),
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
        "M1 numerical parity is evaluated over VCE=0.5-3.0 V against the pre-registered maximum log-error thresholds.",
        "",
        "| Model | Operational | Numerical parity | Sentaurus Ic @ 3 V (A/um) | Vela Ic @ 3 V (A/um) | Ic ratio | Ib ratio | Sentaurus beta | Vela beta |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for model in ("M0", "M1"):
        item = summary[model]
        point = item["at_VCE_3V"]
        lines.append(
            f"| {model} | {item['operational_pass']} | {item['numerical_parity_gate']['pass']} | "
            f"{point['sentaurus_Ic_A_per_um']:.9e} | {point['vela_Ic_A_per_um']:.9e} | "
            f"{point['Ic_magnitude_ratio_vela_over_sentaurus']:.6g} | "
            f"{point['Ib_magnitude_ratio_vela_over_sentaurus']:.6g} | "
            f"{point['sentaurus_beta_abs']:.6g} | {point['vela_beta_abs']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## M1 numerical parity gate",
            "",
            "| Observable | Maximum allowed absolute log10 error | Observed maximum | Pass |",
            "|---|---:|---:|---:|",
        ]
    )
    gate = summary["M1"]["numerical_parity_gate"]
    for name in ("Ic", "Ib", "beta"):
        threshold = gate["thresholds_maximum_absolute_log10_error"][name]
        observed = gate["observed_maximum_absolute_log10_error"][name]
        lines.append(
            f"| {name} | {threshold:.6g} | {observed:.9g} | {observed <= threshold} |"
        )
    lines.extend(
        [
            "",
            f"Overall operational pass: **{summary['operational_pass']}**",
            f"Asserted numerical parity pass: **{summary['numerical_parity_pass']}**",
            f"Overall pass: **{summary['overall_pass']}**",
            "",
        ]
    )
    write_text_lf(path, "\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--vela-run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    threshold_document = json.loads(
        (args.reference_root / "contracts" / "comparison_thresholds.json").read_text(
            encoding="utf-8"
        )
    )
    parity = threshold_document["wp3_wp5_vela_comparison"]["numerical_parity"]
    active_region = tuple(float(value) for value in parity["active_region_VCE_range_V"])
    if len(active_region) != 2 or active_region[0] > active_region[1]:
        raise ValueError("numerical parity active-region range must contain ordered bounds")

    summaries: dict[str, object] = {
        "schema_version": 2,
        "comparison_scope": "WP3-WP5 common-input comparison with pre-registered M1 numerical parity",
    }
    for model, stem in (("M0", "m0"), ("M1", "m1")):
        rows, model_summary = compare_model(
            model,
            args.reference_root / "reference_curves" / f"bjt_{stem}_output.csv",
            args.vela_run_root / f"{stem}_collector_terminal_balance.csv",
            active_region,
            parity["models"][model],
        )
        write_csv(args.output_dir / f"{stem}_sentaurus_vela.csv", rows)
        summaries[model] = model_summary
    summaries["operational_pass"] = bool(
        summaries["M0"]["operational_pass"] and summaries["M1"]["operational_pass"]
    )
    asserted_gates = [
        summaries[model]["numerical_parity_gate"]
        for model in ("M0", "M1")
        if summaries[model]["numerical_parity_gate"]["status"] == "asserted"
    ]
    summaries["numerical_parity_pass"] = bool(asserted_gates) and all(
        gate["pass"] for gate in asserted_gates
    )
    summaries["overall_pass"] = bool(
        summaries["operational_pass"] and summaries["numerical_parity_pass"]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_text_lf(
        args.output_dir / "comparison_summary.json",
        json.dumps(summaries, indent=2, allow_nan=False) + "\n",
    )
    write_markdown(args.output_dir / "comparison_summary.md", summaries)
    print(json.dumps(summaries, indent=2, allow_nan=False))
    return 0 if summaries["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
