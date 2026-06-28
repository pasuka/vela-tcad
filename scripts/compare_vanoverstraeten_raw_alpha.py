#!/usr/bin/env python3
"""Compare default and raw Van Overstraeten alpha against Sentaurus nodes."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path
from typing import Iterable


FIELDS = [
    "bias_V",
    "node_id",
    "x_um",
    "y_um",
    "Fn_sentaurus_V_per_cm",
    "Fn_self_default_V_per_cm",
    "Fn_self_raw_V_per_cm",
    "alpha_n_sentaurus_cm_inv",
    "alpha_n_self_default_cm_inv",
    "alpha_n_self_raw_cm_inv",
    "alpha_n_default_ratio",
    "alpha_n_raw_ratio",
    "Fp_sentaurus_V_per_cm",
    "Fp_self_default_V_per_cm",
    "Fp_self_raw_V_per_cm",
    "alpha_p_sentaurus_cm_inv",
    "alpha_p_self_default_cm_inv",
    "alpha_p_self_raw_cm_inv",
    "alpha_p_default_ratio",
    "alpha_p_raw_ratio",
]

def row_float(row: dict[str, float | str | None], name: str, default: float = 0.0) -> float:
    value = row.get(name)
    return value if isinstance(value, float) else default


REGIONS = [
    ("全部节点", lambda row: True),
    ("高场节点", lambda row: max(row_float(row, "Fn_self_raw_V_per_cm"),
                                  row_float(row, "Fp_self_raw_V_per_cm")) > 1.0e4),
    ("结中心区域", lambda row: 0.9 <= row_float(row, "x_um") <= 1.1),
    ("结左肩部", lambda row: 0.7 <= row_float(row, "x_um") <= 0.85),
    ("结右肩部", lambda row: 1.15 <= row_float(row, "x_um") <= 1.3),
]


def optional_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def gamma_factor(temperature_k: float, reference_temperature_k: float, phonon_energy_ev: float) -> float:
    k_boltzmann_ev_per_k = 8.617333262145e-5
    ref_arg = phonon_energy_ev / (2.0 * k_boltzmann_ev_per_k * reference_temperature_k)
    arg = phonon_energy_ev / (2.0 * k_boltzmann_ev_per_k * temperature_k)
    denom = math.tanh(arg)
    return math.tanh(ref_arg) / denom if abs(denom) > 0.0 else 1.0


def raw_vanoverstraeten_alpha_cm_inv(
    field_v_per_cm: float,
    *,
    low_a_cm_inv: float,
    high_a_cm_inv: float,
    low_b_v_per_cm: float,
    high_b_v_per_cm: float,
    switch_field_v_per_cm: float,
    gamma: float,
) -> float:
    field = abs(field_v_per_cm)
    if field <= 0.0 or gamma <= 0.0:
        return 0.0
    low_field = field < switch_field_v_per_cm
    prefactor = low_a_cm_inv if low_field else high_a_cm_inv
    critical = low_b_v_per_cm if low_field else high_b_v_per_cm
    if prefactor <= 0.0:
        return 0.0
    exponent = max(-critical * gamma / field, -700.0)
    return gamma * prefactor * math.exp(exponent)


def ratio(value: float | None, reference: float | None) -> float | None:
    if value is None or reference is None or reference == 0.0:
        return None
    return value / reference


def fmt(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def read_normalized_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [
            {str(key).strip(): str(value).strip() for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def load_rows(path: Path, args: argparse.Namespace) -> list[dict[str, float | str | None]]:
    source_rows = read_normalized_csv(path)
    if source_rows and {"quantity", "sentaurus_value", "vela_value_scaled_to_sentaurus_units"}.issubset(source_rows[0]):
        return load_long_rows(source_rows, args)
    return load_wide_rows(source_rows, args)


def add_raw_alpha(row: dict[str, float | str | None], args: argparse.Namespace, gamma: float) -> None:
    fn_raw = row.get("Fn_self_raw_V_per_cm")
    fp_raw = row.get("Fp_self_raw_V_per_cm")
    if row.get("alpha_n_self_raw_cm_inv") is None:
        row["alpha_n_self_raw_cm_inv"] = raw_vanoverstraeten_alpha_cm_inv(
            fn_raw if isinstance(fn_raw, float) else 0.0,
            low_a_cm_inv=args.electron_a_low_cm_inv,
            high_a_cm_inv=args.electron_a_high_cm_inv,
            low_b_v_per_cm=args.electron_b_low_v_per_cm,
            high_b_v_per_cm=args.electron_b_high_v_per_cm,
            switch_field_v_per_cm=args.switch_field_v_per_cm,
            gamma=gamma,
        )
    if row.get("alpha_p_self_raw_cm_inv") is None:
        row["alpha_p_self_raw_cm_inv"] = raw_vanoverstraeten_alpha_cm_inv(
            fp_raw if isinstance(fp_raw, float) else 0.0,
            low_a_cm_inv=args.hole_a_low_cm_inv,
            high_a_cm_inv=args.hole_a_high_cm_inv,
            low_b_v_per_cm=args.hole_b_low_v_per_cm,
            high_b_v_per_cm=args.hole_b_high_v_per_cm,
            switch_field_v_per_cm=args.switch_field_v_per_cm,
            gamma=gamma,
        )
    row["alpha_n_default_ratio"] = ratio(
        row.get("alpha_n_self_default_cm_inv"), row.get("alpha_n_sentaurus_cm_inv"))  # type: ignore[arg-type]
    row["alpha_n_raw_ratio"] = ratio(
        row.get("alpha_n_self_raw_cm_inv"), row.get("alpha_n_sentaurus_cm_inv"))  # type: ignore[arg-type]
    row["alpha_p_default_ratio"] = ratio(
        row.get("alpha_p_self_default_cm_inv"), row.get("alpha_p_sentaurus_cm_inv"))  # type: ignore[arg-type]
    row["alpha_p_raw_ratio"] = ratio(
        row.get("alpha_p_self_raw_cm_inv"), row.get("alpha_p_sentaurus_cm_inv"))  # type: ignore[arg-type]


def load_wide_rows(
    source_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict[str, float | str | None]]:
    gamma = gamma_factor(args.temperature_k, args.reference_temperature_k, args.phonon_energy_ev)
    rows: list[dict[str, float | str | None]] = []
    for source in source_rows:
        fn_default = optional_float(source.get("Fn_self_default_V_per_cm"))
        fp_default = optional_float(source.get("Fp_self_default_V_per_cm"))
        fn_raw = optional_float(source.get("Fn_self_raw_V_per_cm"))
        fp_raw = optional_float(source.get("Fp_self_raw_V_per_cm"))
        if fn_raw is None:
            fn_raw = fn_default
        if fp_raw is None:
            fp_raw = fp_default
        alpha_n_sent = optional_float(source.get("alpha_n_sentaurus_cm_inv"))
        alpha_p_sent = optional_float(source.get("alpha_p_sentaurus_cm_inv"))
        alpha_n_default = optional_float(source.get("alpha_n_self_default_cm_inv"))
        alpha_p_default = optional_float(source.get("alpha_p_self_default_cm_inv"))
        row = {
            "bias_V": source.get("bias_V", ""),
            "node_id": source.get("node_id", ""),
            "x_um": optional_float(source.get("x_um")) or 0.0,
            "y_um": optional_float(source.get("y_um")) or 0.0,
            "Fn_sentaurus_V_per_cm": optional_float(source.get("Fn_sentaurus_V_per_cm")),
            "Fn_self_default_V_per_cm": fn_default,
            "Fn_self_raw_V_per_cm": fn_raw,
            "alpha_n_sentaurus_cm_inv": alpha_n_sent,
            "alpha_n_self_default_cm_inv": alpha_n_default,
            "alpha_n_self_raw_cm_inv": optional_float(source.get("alpha_n_self_raw_cm_inv")),
            "alpha_n_default_ratio": None,
            "alpha_n_raw_ratio": None,
            "Fp_sentaurus_V_per_cm": optional_float(source.get("Fp_sentaurus_V_per_cm")),
            "Fp_self_default_V_per_cm": fp_default,
            "Fp_self_raw_V_per_cm": fp_raw,
            "alpha_p_sentaurus_cm_inv": alpha_p_sent,
            "alpha_p_self_default_cm_inv": alpha_p_default,
            "alpha_p_self_raw_cm_inv": optional_float(source.get("alpha_p_self_raw_cm_inv")),
            "alpha_p_default_ratio": None,
            "alpha_p_raw_ratio": None,
        }
        add_raw_alpha(row, args, gamma)
        rows.append(row)
    return rows


def load_long_rows(
    source_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> list[dict[str, float | str | None]]:
    gamma = gamma_factor(args.temperature_k, args.reference_temperature_k, args.phonon_energy_ev)
    grouped: dict[tuple[str, str, str, str], dict[str, float | str | None]] = {}
    for source in source_rows:
        key = (
            source.get("bias_V", ""),
            source.get("node_id", ""),
            source.get("x_um", ""),
            source.get("y_um", ""),
        )
        row = grouped.setdefault(key, {
            "bias_V": source.get("bias_V", ""),
            "node_id": source.get("node_id", ""),
            "x_um": optional_float(source.get("x_um")) or 0.0,
            "y_um": optional_float(source.get("y_um")) or 0.0,
            "Fn_sentaurus_V_per_cm": None,
            "Fn_self_default_V_per_cm": None,
            "Fn_self_raw_V_per_cm": None,
            "alpha_n_sentaurus_cm_inv": None,
            "alpha_n_self_default_cm_inv": None,
            "alpha_n_self_raw_cm_inv": None,
            "alpha_n_default_ratio": None,
            "alpha_n_raw_ratio": None,
            "Fp_sentaurus_V_per_cm": None,
            "Fp_self_default_V_per_cm": None,
            "Fp_self_raw_V_per_cm": None,
            "alpha_p_sentaurus_cm_inv": None,
            "alpha_p_self_default_cm_inv": None,
            "alpha_p_self_raw_cm_inv": None,
            "alpha_p_default_ratio": None,
            "alpha_p_raw_ratio": None,
        })
        quantity = source.get("quantity", "")
        sentaurus = optional_float(source.get("sentaurus_value"))
        vela = optional_float(source.get("vela_value_scaled_to_sentaurus_units"))
        if quantity == "electric_field":
            row["Fn_sentaurus_V_per_cm"] = sentaurus
            row["Fp_sentaurus_V_per_cm"] = sentaurus
            row["Fn_self_default_V_per_cm"] = vela
            row["Fp_self_default_V_per_cm"] = vela
            row["Fn_self_raw_V_per_cm"] = vela
            row["Fp_self_raw_V_per_cm"] = vela
        elif quantity in {"electron_impact_ionization_drive", "electron_avalanche_drive"}:
            row["Fn_sentaurus_V_per_cm"] = sentaurus
            row["Fn_self_default_V_per_cm"] = vela
            row["Fn_self_raw_V_per_cm"] = vela
        elif quantity in {"hole_impact_ionization_drive", "hole_avalanche_drive"}:
            row["Fp_sentaurus_V_per_cm"] = sentaurus
            row["Fp_self_default_V_per_cm"] = vela
            row["Fp_self_raw_V_per_cm"] = vela
        elif quantity == "electron_alpha_avalanche":
            row["alpha_n_sentaurus_cm_inv"] = sentaurus
            row["alpha_n_self_default_cm_inv"] = vela
        elif quantity == "hole_alpha_avalanche":
            row["alpha_p_sentaurus_cm_inv"] = sentaurus
            row["alpha_p_self_default_cm_inv"] = vela

    rows = list(grouped.values())
    for row in rows:
        add_raw_alpha(row, args, gamma)
    return rows


def write_csv(path: Path, rows: Iterable[dict[str, float | str | None]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: fmt(row[name]) if isinstance(row.get(name), float) else row.get(name, "")
                             for name in FIELDS})


def median(values: Iterable[float | None]) -> float | None:
    clean = [value for value in values if value is not None and math.isfinite(value)]
    return statistics.median(clean) if clean else None


def build_summary(rows: list[dict[str, float | str | None]]) -> str:
    lines = [
        "# Van Overstraeten Raw Alpha Compare",
        "",
        "| region | nodes | median alpha_n default/S | median alpha_n raw/S | median alpha_p default/S | median alpha_p raw/S | median raw/default uplift |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    shoulder_uplifts: list[float] = []
    shoulder_raw_ratios: list[float] = []
    for name, predicate in REGIONS:
        subset = [row for row in rows if predicate(row)]  # type: ignore[arg-type]
        n_default = median(row["alpha_n_default_ratio"] for row in subset)  # type: ignore[index]
        n_raw = median(row["alpha_n_raw_ratio"] for row in subset)  # type: ignore[index]
        p_default = median(row["alpha_p_default_ratio"] for row in subset)  # type: ignore[index]
        p_raw = median(row["alpha_p_raw_ratio"] for row in subset)  # type: ignore[index]
        uplifts = []
        for row in subset:
            for carrier in ("n", "p"):
                default = row[f"alpha_{carrier}_self_default_cm_inv"]
                raw = row[f"alpha_{carrier}_self_raw_cm_inv"]
                if isinstance(default, float) and isinstance(raw, float) and default > 0.0:
                    uplifts.append(raw / default)
        uplift = median(uplifts)
        if name in {"结左肩部", "结右肩部"}:
            if uplift is not None:
                shoulder_uplifts.append(uplift)
            for value in (n_raw, p_raw):
                if value is not None:
                    shoulder_raw_ratios.append(value)
        lines.append(
            f"| {name} | {len(subset)} | {fmt(n_default)} | {fmt(n_raw)} | "
            f"{fmt(p_default)} | {fmt(p_raw)} | {fmt(uplift)} |"
        )

    lines.append("")
    uplift = median(shoulder_uplifts)
    raw_ratio = median(shoulder_raw_ratios)
    if uplift is not None and uplift >= 2.0 and (raw_ratio is None or raw_ratio >= 0.25):
        lines.append("alpha 偏小主要来自 cutoff/smoothing/RefDens 抑制")
    else:
        lines.append("问题不在 cutoff/smoothing/RefDens，继续检查公式、参数、分段和单位")
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--elements-csv", type=Path, default=None,
                        help="Reserved for callers that reconstruct gradients before this comparison.")
    parser.add_argument("--output-csv", type=Path,
                        default=Path("build/diagnostics/vanoverstraeten_raw_alpha_compare.csv"))
    parser.add_argument("--summary-md", type=Path,
                        default=Path("build/diagnostics/vanoverstraeten_raw_alpha_summary.md"))
    parser.add_argument("--electron-a-low-cm-inv", type=float, default=7.03e5)
    parser.add_argument("--electron-a-high-cm-inv", type=float, default=7.03e5)
    parser.add_argument("--electron-b-low-v-per-cm", type=float, default=1.231e6)
    parser.add_argument("--electron-b-high-v-per-cm", type=float, default=1.231e6)
    parser.add_argument("--hole-a-low-cm-inv", type=float, default=1.582e6)
    parser.add_argument("--hole-a-high-cm-inv", type=float, default=6.71e5)
    parser.add_argument("--hole-b-low-v-per-cm", type=float, default=2.036e6)
    parser.add_argument("--hole-b-high-v-per-cm", type=float, default=1.693e6)
    parser.add_argument("--switch-field-v-per-cm", type=float, default=4.0e5)
    parser.add_argument("--phonon-energy-ev", type=float, default=0.063)
    parser.add_argument("--reference-temperature-k", type=float, default=300.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = load_rows(args.input_csv, args)
    write_csv(args.output_csv, rows)
    args.summary_md.parent.mkdir(parents=True, exist_ok=True)
    args.summary_md.write_text(build_summary(rows), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
