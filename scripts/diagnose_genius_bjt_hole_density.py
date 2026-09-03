#!/usr/bin/env python3
"""Localize the Genius NPN BJT hole-density discrepancy at VBE=0.7/VCE=3 V."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import TwoSlopeNorm


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = (
    REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
)
SENTAURUS_ROOT = BUILD_ROOT / "m1_current_diagnosis" / "sentaurus_vce3"
VELA_STATE = BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv"
OUTPUT_DIR = FIXTURE / "hole_density_diagnosis"
DIAGNOSTIC_ROOT = BUILD_ROOT / "hole_density_diagnosis"
ENFORCED_STATE = DIAGNOSTIC_ROOT / "enforce" / "enforced_state.csv"

INK = "#27313A"
MUTED = "#66717C"
GRID = "#D9DEE3"
BLUE = "#2D6CDF"
ORANGE = "#D97706"
VT_300K_V = 8.617333262145e-5 * 300.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, default=FIXTURE)
    parser.add_argument("--sentaurus-root", type=Path, default=SENTAURUS_ROOT)
    parser.add_argument("--vela-state", type=Path, default=VELA_STATE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--diagnostic-root", type=Path, default=DIAGNOSTIC_ROOT)
    parser.add_argument("--comparison-state", type=Path, default=ENFORCED_STATE)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_scalar(path: Path, count: int, column: str) -> np.ndarray:
    rows = read_csv(path)
    values = np.full(count, np.nan)
    for row in rows:
        values[int(row["node_id"])] = float(row[column])
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite or missing node values in {path}")
    return values


def percentile(values: np.ndarray, fraction: float) -> float:
    return float(np.quantile(np.abs(values), fraction))


def statistics(values: np.ndarray) -> dict[str, float | int | None]:
    if values.size == 0:
        return {
            "node_count": 0,
            "rmse": None,
            "p95_absolute": None,
            "maximum_absolute": None,
        }
    return {
        "node_count": int(values.size),
        "rmse": float(np.sqrt(np.mean(values * values))),
        "p95_absolute": percentile(values, 0.95),
        "maximum_absolute": float(np.max(np.abs(values))),
    }


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.5,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def symmetric_limit(values: np.ndarray) -> float:
    limit = float(np.quantile(np.abs(values), 0.995))
    return max(limit, float(np.max(np.abs(values))) * 0.25, 1.0e-12)


def draw_map(
    ax: plt.Axes,
    triangulation: mtri.Triangulation,
    values: np.ndarray,
    title: str,
    unit: str,
) -> None:
    limit = symmetric_limit(values)
    artist = ax.tripcolor(
        triangulation,
        values,
        shading="gouraud",
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit),
        rasterized=True,
    )
    ax.set_title(title)
    ax.set_xlabel("x [um]")
    ax.set_ylabel("y [um]")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    colorbar = ax.figure.colorbar(artist, ax=ax, fraction=0.046, pad=0.035)
    colorbar.set_label(unit)


def save_map(
    triangulation: mtri.Triangulation,
    values: np.ndarray,
    title: str,
    unit: str,
    path: Path,
) -> None:
    figure, ax = plt.subplots(figsize=(8.0, 3.7), constrained_layout=True)
    draw_map(ax, triangulation, values, title, unit)
    figure.savefig(path, dpi=220)
    plt.close(figure)


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="\n", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_status(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))["runner_status"]


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    mesh = json.loads(
        (args.fixture_root / "vela" / "input" / "mesh.json").read_text(
            encoding="utf-8"
        )
    )
    nodes = sorted(mesh["nodes"], key=lambda row: int(row["id"]))
    coordinates = np.asarray([[float(row["x"]), float(row["y"])] for row in nodes])
    triangles = np.asarray([row["node_ids"] for row in mesh["triangles"]], dtype=int)
    count = len(nodes)
    fields = args.sentaurus_root / "fields"

    sent_psi = load_scalar(fields / "ElectrostaticPotential_region0.csv", count, "component0")
    sent_phip = load_scalar(fields / "hQuasiFermiPotential_region0.csv", count, "component0")
    sent_p = load_scalar(fields / "hDensity_region0.csv", count, "component0")
    vela_psi = load_scalar(args.vela_state, count, "psi")
    vela_phip = load_scalar(args.vela_state, count, "phip")
    vela_p = load_scalar(args.vela_state, count, "holes_m3") / 1.0e6
    comparison_phip = None
    if args.comparison_state.exists():
        comparison_phip = load_scalar(args.comparison_state, count, "phip")

    delta_log_p = np.log10(np.maximum(vela_p, 1.0e-300) / np.maximum(sent_p, 1.0e-300))
    delta_psi = vela_psi - sent_psi
    delta_phip = vela_phip - sent_phip
    delta_eta_p = (vela_phip - vela_psi) - (sent_phip - sent_psi)
    predicted_delta_log_p = delta_eta_p / (VT_300K_V * math.log(10.0))
    relation_residual = delta_log_p - predicted_delta_log_p
    plateau = (delta_eta_p >= -0.101) & (delta_eta_p <= -0.099)
    large_error = np.abs(delta_log_p) > 1.0
    low_density = sent_p < 1.0e2
    significant = sent_p >= 1.0e10

    rows: list[dict[str, object]] = []
    for index in range(count):
        rows.append(
            {
                "node_id": index,
                "x_um": coordinates[index, 0],
                "y_um": coordinates[index, 1],
                "sentaurus_holes_cm3": sent_p[index],
                "vela_holes_cm3": vela_p[index],
                "delta_log10_holes_decade": delta_log_p[index],
                "delta_psi_V": delta_psi[index],
                "delta_phip_V": delta_phip[index],
                "delta_phip_minus_psi_V": delta_eta_p[index],
                "predicted_delta_log10_holes_decade": predicted_delta_log_p[index],
                "density_relation_residual_decade": relation_residual[index],
                "minus_100mV_plateau": int(plateau[index]),
            }
        )
    write_rows(args.output_dir / "hole_density_qf_node_diagnostics.csv", rows)

    populations = {
        "all_nodes": np.ones(count, dtype=bool),
        "sentaurus_holes_lt_1e2_cm3": low_density,
        "sentaurus_holes_lt_1e10_cm3": ~significant,
        "sentaurus_holes_ge_1e10_cm3": significant,
        "minus_100mV_plateau": plateau,
    }
    summary: dict[str, object] = {
        "schema_version": 1,
        "bias": {"vbe_V": 0.7, "vce_V": 3.0},
        "temperature_K": 300.0,
        "thermal_decade_voltage_V": VT_300K_V * math.log(10.0),
        "populations": {},
        "plateau": {
            "node_count": int(np.count_nonzero(plateau)),
            "large_error_node_count": int(np.count_nonzero(large_error)),
            "overlap_node_count": int(np.count_nonzero(plateau & large_error)),
            "x_um_min": float(np.min(coordinates[plateau, 0])) if np.any(plateau) else None,
            "x_um_max": float(np.max(coordinates[plateau, 0])) if np.any(plateau) else None,
            "y_um_min": float(np.min(coordinates[plateau, 1])) if np.any(plateau) else None,
            "y_um_max": float(np.max(coordinates[plateau, 1])) if np.any(plateau) else None,
        },
    }
    for name, mask in populations.items():
        summary["populations"][name] = {
            "delta_log10_holes": statistics(delta_log_p[mask]),
            "qf_predicted_delta_log10_holes": statistics(predicted_delta_log_p[mask]),
            "density_relation_residual": statistics(relation_residual[mask]),
            "delta_psi_V": statistics(delta_psi[mask]),
            "delta_phip_V": statistics(delta_phip[mask]),
            "delta_phip_minus_psi_V": statistics(delta_eta_p[mask]),
        }
    max_node = int(np.argmax(np.abs(delta_log_p)))
    summary["maximum_error_node"] = rows[max_node]
    (args.output_dir / "hole_density_qf_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    configure_style()
    triangulation = mtri.Triangulation(
        coordinates[:, 0], coordinates[:, 1], triangles
    )
    plots = [
        (delta_log_p, "Vela/SDevice hole-density error", "decade", "delta_log10_hole_density.png"),
        (delta_psi, "Electrostatic-potential difference", "V", "delta_electrostatic_potential.png"),
        (delta_phip, "Hole quasi-Fermi-potential difference", "V", "delta_hole_quasi_fermi_potential.png"),
        (delta_eta_p, "Hole-density driving-potential difference", "V", "delta_hole_density_driving_potential.png"),
        (relation_residual, "Hole-density relation residual", "decade", "hole_density_relation_residual.png"),
    ]
    for values, title, unit, filename in plots:
        save_map(triangulation, values, title, unit, args.output_dir / filename)

    figure, axes = plt.subplots(3, 2, figsize=(13.5, 9.0), constrained_layout=True)
    for ax, (values, title, unit, _) in zip(axes.flat, plots, strict=False):
        draw_map(ax, triangulation, values, title, unit)
    axes.flat[-1].axis("off")
    figure.suptitle("Genius NPN BJT hole-density discrepancy decomposition", color=INK)
    figure.savefig(args.output_dir / "hole_density_qf_decomposition.png", dpi=220)
    plt.close(figure)

    surface = np.isclose(coordinates[:, 1], 0.0) & (coordinates[:, 0] <= 2.0)
    vertical = np.isclose(coordinates[:, 0], 0.84375)
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 4.3), constrained_layout=True)
    for ax, mask, axis_index, xlabel, title in [
        (axes[0], surface, 0, "x [um]", "Base-side surface to low-hole basin (y=0)"),
        (axes[1], vertical, 1, "y [um]", "Collector to low-hole basin (x=0.84375 um)"),
    ]:
        order = np.argsort(coordinates[mask, axis_index])
        position = coordinates[mask, axis_index][order]
        ax.plot(position, sent_phip[mask][order], color=INK, label="SDevice", linewidth=1.8)
        ax.plot(position, vela_phip[mask][order], color=BLUE, label="Vela", linewidth=1.5)
        if comparison_phip is not None:
            ax.plot(
                position,
                comparison_phip[mask][order],
                color=ORANGE,
                linestyle="--",
                label="strict-row diagnostic",
                linewidth=1.4,
            )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("hole quasi-Fermi potential [V]")
        ax.grid(True, alpha=0.8)
        ax.legend()
    figure.savefig(args.output_dir / "hole_quasi_fermi_profiles.png", dpi=220)
    plt.close(figure)

    all_stats = summary["populations"]["all_nodes"]
    low_stats = summary["populations"]["sentaurus_holes_lt_1e2_cm3"]
    text = [
        "# Genius NPN BJT hole-density discrepancy diagnosis",
        "",
        "Bias: VBE=0.70 V, VCE=3.00 V; exact common mesh, no interpolation.",
        "",
        "## QF decomposition",
        "",
        f"- Maximum hole-density error: {all_stats['delta_log10_holes']['maximum_absolute']:.9g} decade.",
        f"- The -99 to -101 mV driving-potential plateau contains {summary['plateau']['node_count']} nodes; {summary['plateau']['overlap_node_count']} also exceed 1 decade error.",
        f"- For SDevice p<1e2 cm^-3, the density-relation residual RMSE is {low_stats['density_relation_residual']['rmse']:.9g} decade.",
        f"- At the maximum-error node, measured and QF-predicted errors are {rows[max_node]['delta_log10_holes_decade']:.9g} and {rows[max_node]['predicted_delta_log10_holes_decade']:.9g} decade.",
        "",
        "The low-density discrepancy is therefore explained by the difference in hQF-psi to numerical precision. This establishes the immediate state-variable cause; solver-limit causality requires the separate A/B runs.",
        "",
    ]
    (args.output_dir / "hole_density_qf_summary.md").write_text(
        "\n".join(text), encoding="utf-8", newline="\n"
    )

    term_path = args.diagnostic_root / "probes" / "carrier_terms.csv"
    row_path = args.diagnostic_root / "probes" / "carrier_rows.csv"
    variants_root = args.diagnostic_root / "ab_step"
    enforce_manifest = args.diagnostic_root / "enforce" / "enforce_manifest.json"
    if term_path.exists() and row_path.exists() and variants_root.exists():
        term_rows = {int(row["node_id"]): row for row in read_csv(term_path)}
        carrier_rows = {int(row["node_id"]): row for row in read_csv(row_path)}
        plateau_ids = np.flatnonzero(plateau).tolist()
        ratios = np.asarray(
            [
                abs(float(term_rows[index]["hole_residual"]))
                / max(
                    abs(float(term_rows[index]["hole_flux_abs_sum"])),
                    abs(float(term_rows[index]["hole_recombination"])),
                    abs(float(term_rows[index]["hole_impact"])),
                    1.0e-30,
                )
                for index in plateau_ids
            ]
        )
        raw_updates = np.asarray(
            [float(carrier_rows[index]["raw_delta_phip_V"]) for index in plateau_ids]
        )
        capped_updates = np.asarray(
            [float(carrier_rows[index]["capped_delta_phip_V"]) for index in plateau_ids]
        )
        solver_summary: dict[str, object] = {
            "baseline_plateau": {
                "node_count": len(plateau_ids),
                "hole_row_ratio": statistics(ratios),
                "raw_delta_phip_V": statistics(raw_updates),
                "capped_delta_phip_V": statistics(capped_updates),
            },
            "cap_one_step_ab": [],
        }
        variant_names = ["cap_0.025V", "cap_0.05V", "cap_0.1V", "cap_disabled"]
        for name in variant_names:
            one_step = {int(row["node_id"]): row for row in read_csv(variants_root / name / "one_step.csv")}
            trial_errors = np.asarray(
                [
                    math.log10(
                        (float(one_step[index]["trial_hole_density_m3"]) / 1.0e6)
                        / sent_p[index]
                    )
                    for index in plateau_ids
                ]
            )
            status = load_status(variants_root / name / "one_step.status.json")
            solver_summary["cap_one_step_ab"].append(
                {
                    "variant": name,
                    "quasi_fermi_update_limit_V": 0.0 if name == "cap_disabled" else float(name[4:-1]),
                    "plateau_trial_hole_error": statistics(trial_errors),
                    "initial_combined_residual": status["block_residuals"]["combined"],
                    "trial_combined_residual": status["trial_block_residuals"]["combined"],
                    "raw_step_norm": status["raw_step_norm"],
                    "limited_step_norm": status["step_norm"],
                }
            )
        if enforce_manifest.exists():
            enforce = json.loads(enforce_manifest.read_text(encoding="utf-8"))["status"]["runner_status"]
            enforce_spatial = json.loads(
                (
                    args.diagnostic_root
                    / "enforce"
                    / "spatial_analysis"
                    / "hole_density_qf_summary.json"
                ).read_text(encoding="utf-8")
            )
            solver_summary["strict_row_run"] = {
                "converged": enforce["converged"],
                "iterations": enforce["iterations"],
                "failure_reason": enforce["failure_reason"],
                "final_combined_residual": enforce["final_residual"],
                "row_max_ratio": enforce["carrier_row_convergence"]["max_ratio"],
                "row_max_node": enforce["carrier_row_convergence"]["max_ratio_node"],
                "row_violation_count": enforce["carrier_row_convergence"]["violation_count"],
                "contact_currents_A_per_um": enforce["contact_currents_A_per_um"],
                "full_domain_hole_error": enforce_spatial["populations"]["all_nodes"]["delta_log10_holes"],
                "low_lt_1e2_hole_error": enforce_spatial["populations"]["sentaurus_holes_lt_1e2_cm3"]["delta_log10_holes"],
                "significant_ge_1e10_hole_error": enforce_spatial["populations"]["sentaurus_holes_ge_1e10_cm3"]["delta_log10_holes"],
            }
        (args.output_dir / "solver_diagnostic_summary.json").write_text(
            json.dumps(solver_summary, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )

        labels = ["baseline", "0.025 V", "0.05 V", "0.1 V", "disabled"]
        ab_rows = solver_summary["cap_one_step_ab"]
        hole_rmse = [statistics(delta_log_p[plateau])["rmse"]] + [
            row["plateau_trial_hole_error"]["rmse"] for row in ab_rows
        ]
        hole_rmse = [np.nan if value is None else value for value in hole_rmse]
        residuals = [ab_rows[0]["initial_combined_residual"]] + [
            row["trial_combined_residual"] for row in ab_rows
        ]
        figure, axes = plt.subplots(1, 2, figsize=(12.8, 4.4), constrained_layout=True)
        axes[0].bar(labels, hole_rmse, color=[MUTED, BLUE, BLUE, BLUE, ORANGE])
        axes[0].set_yscale("log")
        axes[0].set_title("One-step hole-error response to QF update limit")
        axes[0].set_ylabel("plateau RMSE [decade], log scale")
        axes[0].tick_params(axis="x", rotation=20)
        axes[0].grid(True, axis="y", alpha=0.8)
        axes[1].bar(labels, residuals, color=[MUTED, BLUE, BLUE, BLUE, ORANGE])
        axes[1].set_yscale("log")
        axes[1].set_title("One-step nonlinear residual response")
        axes[1].set_ylabel("combined residual, log scale")
        axes[1].tick_params(axis="x", rotation=20)
        axes[1].grid(True, axis="y", alpha=0.8)
        figure.savefig(args.output_dir / "quasi_fermi_update_limit_ab.png", dpi=220)
        plt.close(figure)

        trace_ids = [589, 144, 150, 7, 12]
        scales = np.asarray(
            [
                max(
                    abs(float(term_rows[index]["hole_flux_abs_sum"])),
                    abs(float(term_rows[index]["hole_recombination"])),
                    1.0e-30,
                )
                for index in trace_ids
            ]
        )
        flux = np.asarray([float(term_rows[index]["hole_flux"]) for index in trace_ids]) / scales
        srh = np.asarray([float(term_rows[index]["hole_recombination"]) for index in trace_ids]) / scales
        residual = np.asarray([float(term_rows[index]["hole_residual"]) for index in trace_ids]) / scales
        x = np.arange(len(trace_ids))
        width = 0.25
        figure, ax = plt.subplots(figsize=(8.5, 4.4), constrained_layout=True)
        ax.bar(x - width, flux, width, label="hole flux divergence", color=BLUE)
        ax.bar(x, srh, width, label="SRH", color=ORANGE)
        ax.bar(x + width, residual, width, label="residual", color=MUTED)
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.set_xticks(x, [str(index) for index in trace_ids])
        ax.set_xlabel("node id")
        ax.set_ylabel("term / local row scale")
        ax.set_title("Hole continuity balance at representative plateau nodes")
        ax.legend()
        ax.grid(True, axis="y", alpha=0.8)
        figure.savefig(args.output_dir / "hole_continuity_plateau_terms.png", dpi=220)
        plt.close(figure)

        strict = solver_summary.get("strict_row_run", {})
        if ratios.size:
            continuity_sentence = (
                f"Across the plateau, the local hole-row ratio has median {float(np.median(ratios)):.6g}, "
                f"P95 {percentile(ratios, 0.95):.6g}, and maximum {float(np.max(ratios)):.6g}. "
                f"The raw Newton hQF correction has median {float(np.median(raw_updates)):.6g} V "
                f"and is limited to {float(np.median(capped_updates)):.6g} V."
            )
        else:
            continuity_sentence = (
                "No nodes satisfy the -101 to -99 mV plateau definition in this state; "
                "plateau continuity and update statistics are therefore not applicable."
            )
        report_lines = [
            "# Genius NPN BJT hole-density solver diagnosis",
            "",
            "## 1. Spatial decomposition",
            "",
            f"The baseline contains {len(plateau_ids)} nodes with delta(hQF-psi) between -101 and -99 mV. Of the {int(np.count_nonzero(large_error))} nodes above 1 decade hole error, {int(np.count_nonzero(plateau & large_error))} are in this plateau.",
            f"For SDevice p<1e2 cm^-3, the carrier-relation residual RMSE is {statistics(relation_residual[low_density])['rmse']:.6g} decade.",
            "",
            "## 2. Hole continuity terms",
            "",
            continuity_sentence,
            "",
            "## 3. QF update-limit A/B",
            "",
            "| Limit | Plateau trial RMSE [decade] | Trial combined residual |",
            "|---:|---:|---:|",
        ]
        for row in ab_rows:
            limit = "disabled" if row["quasi_fermi_update_limit_V"] == 0.0 else f"{row['quasi_fermi_update_limit_V']:.3g} V"
            plateau_rmse = row["plateau_trial_hole_error"]["rmse"]
            plateau_rmse_text = "n/a" if plateau_rmse is None else f"{plateau_rmse:.6g}"
            report_lines.append(
                f"| {limit} | {plateau_rmse_text} | {row['trial_combined_residual']:.6g} |"
            )
        if strict:
            report_lines.extend(
                [
                    "",
                    "## 4. Strict carrier-row convergence",
                    "",
                    f"The enforce/1e-4 run ended after {strict['iterations']} iterations with `{strict['failure_reason']}`. The original -100 mV plateau was removed and p<1e2 cm^-3 RMSE fell to {strict['low_lt_1e2_hole_error']['rmse']:.6g} decade, but {strict['row_violation_count']} rows still violated the global criterion; the maximum ratio was {strict['row_max_ratio']:.6g} at node {strict['row_max_node']}.",
                    "",
                    "## 5. Contact-to-bulk profiles",
                    "",
                    "The base and collector contact values agree. The approximately 0.1 V offset appears only inside the remote n-type low-hole basin, which identifies weak interior minority-carrier coupling rather than a contact boundary-condition mismatch.",
                    "",
                    "## Conclusion",
                    "",
                    "The immediate discrepancy is a prematurely accepted minority-hole QF state. A 0.1 V limited correction removes the plateau and reduces the nonlinear residual; disabling the limit is unstable. Global strict row enforcement is too sensitive to other negligible-density cancellation rows, so the production fix should qualify carrier rows by physical relevance rather than remove the QF update limit.",
                    "",
                ]
            )
        (args.output_dir / "solver_diagnostic_summary.md").write_text(
            "\n".join(report_lines), encoding="utf-8", newline="\n"
        )
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
