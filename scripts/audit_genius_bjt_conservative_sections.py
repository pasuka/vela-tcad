#!/usr/bin/env python3
"""Audit conservative BJT SG fluxes across base and collector-side sections."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_INDICES = (0, 10, 20, 30)
DEFAULT_CUTS_UM = (("base_side", 0.45), ("collector_side", 0.85))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def section_current(
    edge_rows: list[dict[str, str]], *, axis: str, cut_um: float
) -> dict[str, float | int]:
    """Return current toward increasing ``axis`` through a discrete dual-face cut."""
    coordinate = "x" if axis == "x" else "y"
    electron_particles = 0.0
    hole_particles = 0.0
    crossing_edges = 0
    for row in edge_rows:
        c0 = float(row[f"{coordinate}0"]) * 1.0e6
        c1 = float(row[f"{coordinate}1"]) * 1.0e6
        low0 = c0 <= cut_um
        low1 = c1 <= cut_um
        if low0 == low1:
            continue
        sign = 1.0 if low0 else -1.0
        electron_particles += sign * float(
            row["electron_particle_line_flux_per_m_s"]
        )
        hole_particles += sign * float(row["hole_particle_line_flux_per_m_s"])
        crossing_edges += 1
    electron = -Q_C * electron_particles * 1.0e-6
    hole = Q_C * hole_particles * 1.0e-6
    return {
        "crossing_edge_count": crossing_edges,
        "electron_A_per_um": electron,
        "hole_A_per_um": hole,
        "total_A_per_um": electron + hole,
    }


def relative_error(actual: float, reference: float, floor: float = 1.0e-30) -> float:
    return abs(actual - reference) / max(abs(reference), floor)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--edge-root", type=Path, default=BUILD_ROOT / "m1_p0_multibias" / "edge_audit"
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "transport_sources" / "sources",
    )
    parser.add_argument(
        "--terminal-csv",
        type=Path,
        default=BUILD_ROOT / "m1_p0_acceptance" / "terminal_currents.csv",
    )
    parser.add_argument(
        "--comparison-csv",
        type=Path,
        default=BUILD_ROOT / "m1_p0_relative_error" / "overall" / "m1_sentaurus_vela.csv",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "conservative_sections",
    )
    parser.add_argument(
        "--thresholds",
        type=Path,
        default=REPO / "reference_tcad" / "genius_bjt_sentaurus2022" / "contracts" / "comparison_thresholds.json",
    )
    parser.add_argument("--indices", type=int, nargs="+", default=list(DEFAULT_INDICES))
    args = parser.parse_args()

    terminals = {int(row["index"]): row for row in read_csv(args.terminal_csv)}
    comparisons = {
        round(float(row["VCE_V"]), 12): row for row in read_csv(args.comparison_csv)
    }
    output_rows: list[dict[str, object]] = []
    point_reports: list[dict[str, object]] = []

    for index in args.indices:
        token = f"vce_{index:03d}"
        edge_rows = read_csv(args.edge_root / token / "sg_edges.csv")
        source_rows = read_csv(args.source_root / f"{token}.csv")
        terminal = terminals[index]
        vce = float(terminal["vce_V"])
        comparison = comparisons[round(vce, 12)]
        collector = float(terminal["collector_A_per_um"])
        sentaurus_collector = float(comparison["sentaurus_Ic_A_per_um"])

        sections: dict[str, dict[str, float | int]] = {}
        for name, cut_um in DEFAULT_CUTS_UM:
            section = section_current(edge_rows, axis="y", cut_um=cut_um)
            section["cut_um"] = cut_um
            section["terminal_closure_relative_error"] = relative_error(
                -float(section["total_A_per_um"]), collector
            )
            section["sdevice_terminal_absolute_log10_error"] = abs(
                math.log10(
                    abs(float(section["total_A_per_um"]) / sentaurus_collector)
                )
            )
            sections[name] = section
            output_rows.append(
                {
                    "VCE_V": vce,
                    "section": name,
                    "cut_y_um": cut_um,
                    **section,
                    "vela_collector_A_per_um": collector,
                    "sdevice_collector_A_per_um": sentaurus_collector,
                }
            )

        upper = sections["base_side"]
        lower = sections["collector_side"]
        recombination = sum(
            float(row["total_integrated_A_per_um"])
            for row in source_rows
            if float(upper["cut_um"]) < float(row["y_um"]) <= float(lower["cut_um"])
        )
        electron_change = float(lower["electron_A_per_um"]) - float(
            upper["electron_A_per_um"]
        )
        hole_change = float(lower["hole_A_per_um"]) - float(
            upper["hole_A_per_um"]
        )
        total_spread = abs(
            float(lower["total_A_per_um"]) - float(upper["total_A_per_um"])
        )
        scale = max(abs(collector), 1.0e-30)
        point_reports.append(
            {
                "index": index,
                "VCE_V": vce,
                "sections": sections,
                "slab": {
                    "bounds_y_um": [upper["cut_um"], lower["cut_um"]],
                    "recombination_integral_A_per_um": recombination,
                    "electron_current_change_A_per_um": electron_change,
                    "hole_current_change_A_per_um": hole_change,
                    "electron_continuity_relative_error": relative_error(
                        electron_change, recombination, 1.0e-30
                    ),
                    "hole_continuity_relative_error": relative_error(
                        -hole_change, recombination, 1.0e-30
                    ),
                    "total_current_spread_A_per_um": total_spread,
                    "total_current_spread_relative_to_terminal": total_spread / scale,
                },
            }
        )

    terminal_closures = [
        float(section["terminal_closure_relative_error"])
        for point in point_reports
        for section in point["sections"].values()
    ]
    sdevice_errors = [
        float(section["sdevice_terminal_absolute_log10_error"])
        for point in point_reports
        for section in point["sections"].values()
    ]
    electron_closures = [
        float(point["slab"]["electron_continuity_relative_error"])
        for point in point_reports
    ]
    hole_closures = [
        float(point["slab"]["hole_continuity_relative_error"])
        for point in point_reports
    ]
    spreads = [
        float(point["slab"]["total_current_spread_relative_to_terminal"])
        for point in point_reports
    ]
    thresholds = json.loads(args.thresholds.read_text(encoding="utf-8"))[
        "wp3_wp5_vela_comparison"
    ]["conservative_section_flux_gate"]["thresholds"]
    checks = {
        "section_to_vela_terminal": max(terminal_closures) <= thresholds[
            "maximum_section_to_vela_terminal_relative_error"
        ],
        "section_to_sdevice_terminal": max(sdevice_errors) <= thresholds[
            "maximum_section_to_sdevice_terminal_absolute_log10_error"
        ],
        "electron_continuity": max(electron_closures) <= thresholds[
            "maximum_carrier_continuity_relative_error"
        ],
        "hole_continuity": max(hole_closures) <= thresholds[
            "maximum_carrier_continuity_relative_error"
        ],
        "section_total_current_spread": max(spreads) <= thresholds[
            "maximum_section_total_current_spread_relative_to_terminal"
        ],
    }
    report = {
        "schema_version": 1,
        "orientation": "positive section current points from top (y=0) toward collector (increasing y)",
        "method": "sum production SG particle line fluxes across the node partition; no nodal-vector reconstruction",
        "thresholds": thresholds,
        "observed_maxima": {
            "section_to_vela_terminal_relative_error": max(terminal_closures),
            "section_to_sdevice_terminal_absolute_log10_error": max(sdevice_errors),
            "electron_continuity_relative_error": max(electron_closures),
            "hole_continuity_relative_error": max(hole_closures),
            "section_total_current_spread_relative_to_terminal": max(spreads),
        },
        "checks": checks,
        "pass": all(checks.values()),
        "points": point_reports,
    }

    args.output_root.mkdir(parents=True, exist_ok=True)
    with (args.output_root / "section_currents.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)
    (args.output_root / "conservative_section_acceptance.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    summary_rows = [
        {
            "VCE_V": point["VCE_V"],
            "electron_closure_error": point["slab"]["electron_continuity_relative_error"],
            "hole_closure_error": point["slab"]["hole_continuity_relative_error"],
            "total_current_spread": point["slab"]["total_current_spread_relative_to_terminal"],
        }
        for point in point_reports
    ]
    with (args.output_root / "conservative_section_points.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(json.dumps({"pass": report["pass"], **report["observed_maxima"]}))
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
