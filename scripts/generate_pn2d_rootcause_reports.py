#!/usr/bin/env python3
"""Generate pn2d IV/BV root-cause report artifacts.

This script is a thin analysis/probing layer. It does not modify solver APIs
or default pn2d deck behavior.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import subprocess
import time
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]


def read_csv(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace").replace("\x00", "")
    return list(csv.DictReader(io.StringIO(text)))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def f(row: dict[str, str], key: str, default: float = 0.0) -> float:
    raw = row.get(key, "")
    if raw is None or raw == "":
        return default
    return float(raw)


def strict_newton(rows: list[dict[str, str]]) -> bool:
    for row in rows:
        if row.get("solver_method") != "gummel_newton":
            return False
        if row.get("handoff_stage") != "newton":
            return False
        if int(row.get("newton_iterations", "0")) <= 0:
            return False
    return True


def interp(rows: list[dict[str, str]], x_key: str, y_key: str, x: float) -> float:
    pts: list[tuple[float, float]] = []
    for row in rows:
        try:
            xv = float(row.get(x_key, ""))
            yv = float(row.get(y_key, ""))
        except (TypeError, ValueError):
            continue
        pts.append((xv, yv))
    pts.sort()
    if not pts:
        raise ValueError(f"No points in curve for {y_key}")
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for i in range(1, len(pts)):
        x0, y0 = pts[i - 1]
        x1, y1 = pts[i]
        if x0 <= x <= x1:
            if x1 == x0:
                return y0
            t = (x - x0) / (x1 - x0)
            return y0 + t * (y1 - y0)
    return pts[-1][1]


def parse_threshold_from_tag(tag: str) -> float | None:
    if "th0p" not in tag:
        return None
    try:
        part = tag.split("th0p", 1)[1]
        digits = ""
        for ch in part:
            if ch.isdigit():
                digits += ch
            else:
                break
        if not digits:
            return None
        return float(f"0.{digits}")
    except Exception:
        return None


def orders_of_magnitude(ref: float, cand: float) -> float:
    eps = 1.0e-300
    ar = max(abs(ref), eps)
    ac = max(abs(cand), eps)
    return abs(math.log10(ac / ar))


def run_case(runner: Path, config: Path, cwd: Path, timeout_s: int = 45) -> str:
    try:
        subprocess.run(
            [str(runner), "--config", str(config)],
            cwd=cwd,
            check=True,
            timeout=timeout_s,
        )
        return "ok"
    except subprocess.TimeoutExpired:
        return "timeout"
    except subprocess.CalledProcessError:
        return "failed"


def build_reports(root: Path, runner: Path, reuse_bv_curves: bool = False) -> None:
    reports_dir = root / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    contact_summary = read_csv(root / "pn2d_contact_relax_summary.csv")
    baseline_row = next((r for r in contact_summary if r.get("tag") == "baseline"), None)
    if baseline_row is None:
        raise RuntimeError("Missing baseline row in pn2d_contact_relax_summary.csv")

    ref_iv_csv = root / "reference" / "vela" / "pn2d_iv.csv"
    ref_bv_csv = root / "reference" / "vela" / "pn2d_bv.csv"
    ref_iv_rows = read_csv(ref_iv_csv)
    ref_bv_rows = read_csv(ref_bv_csv)

    baseline_summary = {
        "source": str(root),
        "checks": {
            "iv_orders_0p2_to_0p3": float(baseline_row["iv_window_orders_0p2_to_0p3"]),
            "iv_pass_orders_le_0p50": float(baseline_row["iv_window_orders_0p2_to_0p3"]) <= 0.50,
            "bv_orders_at_0p05": float(baseline_row["bv_orders_at_0p05"]),
            "bv_pass_orders_lt_0p15": float(baseline_row["bv_orders_at_0p05"]) < 0.15,
            "strict_newton_iv": strict_newton(ref_iv_rows),
            "strict_newton_bv": strict_newton(ref_bv_rows),
            "strict_newton_all": strict_newton(ref_iv_rows) and strict_newton(ref_bv_rows),
            "terminal_sum_abs_A_per_um_at_0p3": float(baseline_row["terminal_sum_abs_A_per_um_at_0p3"]),
            "terminal_sum_pass_le_1e_18": float(baseline_row["terminal_sum_abs_A_per_um_at_0p3"]) <= 1.0e-18,
        },
    }
    (reports_dir / "fresh_baseline_summary.json").write_text(
        json.dumps(baseline_summary, indent=2) + "\n", encoding="utf-8"
    )

    iv_matrix_rows: list[dict[str, Any]] = []
    for row in contact_summary:
        tag = row["tag"]
        collapsed = tag.startswith("n_only_")
        iv_matrix_rows.append(
            {
                "tag": tag,
                "description": row["description"],
                "iv_orders_0p2_to_0p3": row["iv_window_orders_0p2_to_0p3"],
                "bv_orders_at_0p05": row["bv_orders_at_0p05"],
                "candidate_ratio_0p29_over_0p30": row["candidate_ratio_0p29_over_0p30"],
                "ratio_abs_delta_to_reference": row["ratio_abs_delta_to_reference"],
                "terminal_sum_abs_A_per_um_at_0p3": row["terminal_sum_abs_A_per_um_at_0p3"],
                "strict_newton_handoff_all": row["strict_newton_handoff_all"],
                "n_only_axis_collapsed": str(collapsed).lower(),
                "n_only_collapse_reason": (
                    "Anode-driven IV keeps Cathode local Vbias at 0 V, so n_contact_only"
                    " threshold/strength do not provide a continuous effective sweep axis"
                    if collapsed
                    else ""
                ),
            }
        )
    write_csv(
        reports_dir / "iv_candidate_matrix.csv",
        iv_matrix_rows,
        [
            "tag",
            "description",
            "iv_orders_0p2_to_0p3",
            "bv_orders_at_0p05",
            "candidate_ratio_0p29_over_0p30",
            "ratio_abs_delta_to_reference",
            "terminal_sum_abs_A_per_um_at_0p3",
            "strict_newton_handoff_all",
            "n_only_axis_collapsed",
            "n_only_collapse_reason",
        ],
    )

    iv_sem_rows: list[dict[str, Any]] = []
    for row in contact_summary:
        tag = row["tag"]
        threshold = parse_threshold_from_tag(tag)
        cdir = root / "candidates" / tag
        for side in ["cathode", "anode"]:
            path = cdir / f"pn2d_iv_fine_{side}_{tag}.csv"
            if not path.exists():
                continue
            for curve_row in read_csv(path):
                bias = f(curve_row, "bias_V")
                bias_contact = curve_row.get("bias_contact", "")
                current_contact = curve_row.get("current_contact", "")
                contact = current_contact
                local_vbias = bias if contact == bias_contact else 0.0
                contact_polarity = "n_contact" if contact.lower() == "cathode" else "p_contact"
                n_only_active = (
                    tag.startswith("n_only_")
                    and contact_polarity == "n_contact"
                    and threshold is not None
                    and local_vbias >= threshold
                )
                iv_sem_rows.append(
                    {
                        "candidate_tag": tag,
                        "bias_V": bias,
                        "bias_contact": bias_contact,
                        "contact_under_test": contact,
                        "contact_polarity": contact_polarity,
                        "local_Vbias_V": local_vbias,
                        "terminal_voltage_difference_V": bias,
                        "n_only_threshold_V": "" if threshold is None else threshold,
                        "n_only_condition_active": str(bool(n_only_active)).lower(),
                        "degeneracy_flag": str(
                            bool(tag.startswith("n_only_") and contact_polarity == "n_contact" and abs(local_vbias) < 1.0e-15)
                        ).lower(),
                        "degeneracy_reason": (
                            "Anode-biased sweep grounds Cathode locally, collapsing n_contact_only threshold/strength axis"
                            if tag.startswith("n_only_") and contact_polarity == "n_contact" and abs(local_vbias) < 1.0e-15
                            else ""
                        ),
                        "current_total_A_per_um": f(curve_row, "current_total_A_per_um"),
                        "current_electron_A_per_um": f(curve_row, "current_electron_A_per_um"),
                        "current_hole_A_per_um": f(curve_row, "current_hole_A_per_um"),
                        "current_electron_drift_A_per_um": f(curve_row, "current_electron_drift_A_per_um"),
                        "current_electron_diffusion_A_per_um": f(curve_row, "current_electron_diffusion_A_per_um"),
                        "current_hole_drift_A_per_um": f(curve_row, "current_hole_drift_A_per_um"),
                        "current_hole_diffusion_A_per_um": f(curve_row, "current_hole_diffusion_A_per_um"),
                    }
                )
    write_csv(
        reports_dir / "iv_contact_boundary_semantics.csv",
        iv_sem_rows,
        [
            "candidate_tag",
            "bias_V",
            "bias_contact",
            "contact_under_test",
            "contact_polarity",
            "local_Vbias_V",
            "terminal_voltage_difference_V",
            "n_only_threshold_V",
            "n_only_condition_active",
            "degeneracy_flag",
            "degeneracy_reason",
            "current_total_A_per_um",
            "current_electron_A_per_um",
            "current_hole_A_per_um",
            "current_electron_drift_A_per_um",
            "current_electron_diffusion_A_per_um",
            "current_hole_drift_A_per_um",
            "current_hole_diffusion_A_per_um",
        ],
    )

    ref_vela_dir = root / "reference" / "vela"
    base_cfg_path = ref_vela_dir / "simulation_bv.json"
    base_cfg = json.loads(base_cfg_path.read_text(encoding="utf-8"))
    ref_bv_curve = read_csv(root / "reference" / "reference_curves" / "pn2d_bv_reference.csv")

    cases = [
        {"name": "m_default_bgn_none_srh_none_ii_none", "mobility": "default", "bgn": "none", "recomb": ["srh"], "ii": "none"},
        {"name": "m_ct_bgn_none_srh_none_ii_none", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["srh"], "ii": "none"},
        {"name": "m_ct_bgn_slotboom_srh_none_ii_none", "mobility": "caughey_thomas", "bgn": "slotboom", "recomb": ["srh"], "ii": "none"},
        {"name": "m_ct_bgn_none_recomb_none_ii_none", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["none"], "ii": "none"},
        {"name": "m_ct_bgn_none_recomb_srh_ii_none", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["srh"], "ii": "none"},
        {"name": "m_ct_bgn_none_recomb_srh_auger_ii_none", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["srh", "auger"], "ii": "none"},
        {"name": "m_ct_bgn_none_recomb_auger_ii_none", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["auger"], "ii": "none"},
        {"name": "m_ct_bgn_none_recomb_srh_ii_selberherr", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["srh"], "ii": "selberherr"},
        {"name": "m_ct_bgn_none_recomb_srh_auger_ii_selberherr", "mobility": "caughey_thomas", "bgn": "none", "recomb": ["srh", "auger"], "ii": "selberherr"},
    ]

    probe_biases = [0.05, 0.10, 0.20, 0.50]
    bv_matrix_rows: list[dict[str, Any]] = []
    bv_decomp_rows: list[dict[str, Any]] = []
    case_cfg_dir = reports_dir / "bv_case_configs"
    case_cfg_dir.mkdir(parents=True, exist_ok=True)
    case_csv_dir = reports_dir / "bv_case_curves"
    case_csv_dir.mkdir(parents=True, exist_ok=True)
    run_id = str(int(time.time()))

    for case in cases:
        cfg = json.loads(json.dumps(base_cfg))
        for key in [
            "mesh_file",
            "node_doping_file",
            "node_mobility_file",
            "node_lifetime_file",
            "interface_charge_file",
        ]:
            raw = cfg.get(key)
            if not isinstance(raw, str):
                continue
            p = Path(raw)
            if p.is_absolute():
                continue
            candidate = (ref_vela_dir / p).resolve()
            if candidate.exists():
                cfg[key] = str(candidate)
        cfg["output_csv"] = ""
        cfg.setdefault("solver", {})["method"] = "gummel_newton"
        cfg["solver"]["recombination"] = case["recomb"]
        cfg["solver"]["bandgap_narrowing"] = case["bgn"]
        if case["mobility"] == "default":
            cfg["solver"].pop("mobility", None)
        cfg["solver"].setdefault("impact_ionization", {})
        cfg["solver"]["impact_ionization"]["model"] = case["ii"]
        cfg.setdefault("sweep", {})["start"] = 0.0
        cfg["sweep"]["stop"] = 0.05
        cfg["sweep"]["step"] = 0.05
        cfg["sweep"]["contact"] = "Cathode"
        cfg["sweep"]["current_contact"] = "Cathode"

        cfg_path = case_cfg_dir / f"simulation_bv_{case['name']}.json"
        for b in probe_biases:
            cfg_b = json.loads(json.dumps(cfg))
            cfg_b["sweep"]["start"] = 0.0
            cfg_b["sweep"]["stop"] = b
            cfg_b["sweep"]["step"] = max(b, 0.05)
            bias_tag = str(b).replace('.', 'p')
            out_csv = (case_csv_dir / f"pn2d_bv_{case['name']}_{run_id}_{bias_tag}.csv").resolve()
            run_status = "missing"
            if reuse_bv_curves:
                matches = sorted(case_csv_dir.glob(f"pn2d_bv_{case['name']}_*_{bias_tag}.csv"))
                if matches:
                    out_csv = matches[-1]
                    run_status = "reused"
            else:
                cfg_b["output_csv"] = str(out_csv)
                cfg_path.write_text(json.dumps(cfg_b, indent=2) + "\n", encoding="utf-8")
                run_status = run_case(runner=runner, config=cfg_path, cwd=ref_vela_dir)

            have_curve = out_csv.exists() and out_csv.stat().st_size > 0
            if have_curve:
                curve = read_csv(out_csv)
                cand_total = interp(curve, "bias_V", "current_total_A_per_um", b)
                ref_total = interp(ref_bv_curve, "bias_V", "current_total", b)
                cand_e = interp(curve, "bias_V", "current_electron_A_per_um", b)
                cand_h = interp(curve, "bias_V", "current_hole_A_per_um", b)
                cand_ed = interp(curve, "bias_V", "current_electron_drift_A_per_um", b)
                cand_ef = interp(curve, "bias_V", "current_electron_diffusion_A_per_um", b)
                cand_hd = interp(curve, "bias_V", "current_hole_drift_A_per_um", b)
                cand_hf = interp(curve, "bias_V", "current_hole_diffusion_A_per_um", b)
                orders = orders_of_magnitude(ref_total, cand_total)
                cancel = abs(cand_ed) / max(abs(cand_e), 1.0e-300)
                residual = cand_total - (cand_e + cand_h)
            else:
                ref_total = interp(ref_bv_curve, "bias_V", "current_total", b)
                cand_total = math.nan
                cand_e = math.nan
                cand_h = math.nan
                cand_ed = math.nan
                cand_ef = math.nan
                cand_hd = math.nan
                cand_hf = math.nan
                orders = math.nan
                cancel = math.nan
                residual = math.nan

            bv_matrix_rows.append(
                {
                    "case": case["name"],
                    "bias_V": b,
                    "mobility": case["mobility"],
                    "bandgap_narrowing": case["bgn"],
                    "recombination": "+".join(case["recomb"]),
                    "impact_ionization": case["ii"],
                    "run_status": run_status,
                    "reference_total": ref_total,
                    "candidate_total_A_per_um": cand_total,
                    "candidate_electron_A_per_um": cand_e,
                    "candidate_hole_A_per_um": cand_h,
                    "orders_of_magnitude": orders,
                }
            )
            bv_decomp_rows.append(
                {
                    "case": case["name"],
                    "bias_V": b,
                    "run_status": run_status,
                    "candidate_total_A_per_um": cand_total,
                    "electron_total_A_per_um": cand_e,
                    "electron_drift_A_per_um": cand_ed,
                    "electron_diffusion_A_per_um": cand_ef,
                    "hole_total_A_per_um": cand_h,
                    "hole_drift_A_per_um": cand_hd,
                    "hole_diffusion_A_per_um": cand_hf,
                    "electron_cancel_ratio_abs_drift_over_total": cancel,
                    "terminal_sum_residual_A_per_um": residual,
                }
            )

    write_csv(
        reports_dir / "bv_physics_matrix.csv",
        bv_matrix_rows,
        [
            "case",
            "bias_V",
            "mobility",
            "bandgap_narrowing",
            "recombination",
            "impact_ionization",
            "run_status",
            "reference_total",
            "candidate_total_A_per_um",
            "candidate_electron_A_per_um",
            "candidate_hole_A_per_um",
            "orders_of_magnitude",
        ],
    )
    write_csv(
        reports_dir / "bv_drift_diff_decomposition.csv",
        bv_decomp_rows,
        [
            "case",
            "bias_V",
            "run_status",
            "candidate_total_A_per_um",
            "electron_total_A_per_um",
            "electron_drift_A_per_um",
            "electron_diffusion_A_per_um",
            "hole_total_A_per_um",
            "hole_drift_A_per_um",
            "hole_diffusion_A_per_um",
            "electron_cancel_ratio_abs_drift_over_total",
            "terminal_sum_residual_A_per_um",
        ],
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=REPO / "build" / "pn2d_curve_rootcause_20260527",
        help="Root directory holding fresh reference and candidate scan outputs.",
    )
    parser.add_argument(
        "--runner",
        type=Path,
        default=REPO / "build" / "vela_example_runner.exe",
        help="Path to vela_example_runner executable.",
    )
    parser.add_argument(
        "--reuse-bv-curves",
        action="store_true",
        help="Do not run new BV cases; reuse existing per-bias BV case curves if available.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_reports(
        root=args.root.resolve(),
        runner=args.runner.resolve(),
        reuse_bv_curves=bool(args.reuse_bv_curves),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())