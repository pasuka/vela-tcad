#!/usr/bin/env python3
"""Generate pn2d IV/BV root-cause report artifacts.

This script is an analysis/probing layer only. It does not modify solver APIs
or default pn2d deck behavior.
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import json
import math
import re
import subprocess
import time
from collections import defaultdict
from dataclasses import dataclass
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


def finite(x: float) -> bool:
    return not (math.isnan(x) or math.isinf(x))


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
        return math.nan
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


def find_row_at_bias(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    best: dict[str, str] | None = None
    best_err = 1.0e300
    for row in rows:
        try:
            b = float(row.get("bias_V", ""))
        except (TypeError, ValueError):
            continue
        err = abs(b - bias)
        if err < best_err:
            best = row
            best_err = err
    if best is None:
        return None
    if best_err > max(1.0e-12, abs(bias) * 1.0e-9):
        return None
    return best


@dataclass
class RunResult:
    status: str
    reason: str


def _to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def normalize_run_failure_reason(reason: str) -> str:
    if reason == "runner_timeout":
        return "runner_timeout"
    return "runner_nonzero_exit_with_json"


def normalize_failure_reason(row: dict[str, Any]) -> str:
    failure_class = str(row.get("failure_class", ""))
    if failure_class == "run_failure":
        return normalize_run_failure_reason(str(row.get("failure_reason", "")))
    return str(row.get("failure_reason", ""))


def effective_row_settings(row: dict[str, Any]) -> dict[str, str]:
    mobility = str(row.get("effective_mobility_model") or row.get("mobility") or "default")
    recomb = str(row.get("effective_recombination_models") or row.get("recombination") or "none")
    bgn = str(
        row.get("effective_bandgap_narrowing")
        or row.get("bandgap_narrowing")
        or "default"
    )
    impact = str(
        row.get("effective_impact_ionization_model")
        or row.get("impact_ionization")
        or "none"
    )
    return {
        "mobility": mobility,
        "recombination": recomb,
        "bandgap_narrowing": bgn,
        "impact_ionization": impact,
    }


def baseline_equivalent_row(row: dict[str, Any]) -> bool:
    settings = effective_row_settings(row)
    return (
        settings["mobility"] == "caughey_thomas"
        and settings["recombination"] == "srh"
        and settings["impact_ionization"] == "none"
        and settings["bandgap_narrowing"] in {"none", "default", "default_inherited"}
    )


def summarize_bv_matrix_rows(
    matrix_rows: list[dict[str, Any]],
    bias_values: list[float],
) -> dict[str, Any]:
    failure_class_counts: dict[str, int] = {
        "executed": 0,
        "non_convergence": 0,
        "run_failure": 0,
        "not_executed": 0,
    }
    failure_reason_counts: dict[str, int] = defaultdict(int)
    per_bias_rows: dict[float, list[dict[str, Any]]] = defaultdict(list)

    for row in matrix_rows:
        failure_class = str(row.get("failure_class", "not_executed"))
        failure_class_counts[failure_class] = failure_class_counts.get(failure_class, 0) + 1
        reason = normalize_failure_reason(row)
        if reason:
            failure_reason_counts[reason] += 1

        bias = _to_float(row.get("bias_V"))
        if finite(bias):
            per_bias_rows[bias].append(row)

    per_bias_summary: list[dict[str, Any]] = []
    for bias in sorted(set(bias_values) | set(per_bias_rows.keys())):
        rows = per_bias_rows.get(bias, [])
        attempted = len(rows)
        executed_candidates: list[tuple[float, dict[str, Any]]] = []
        for row in rows:
            if str(row.get("failure_class")) != "executed":
                continue
            orders = _to_float(row.get("orders_of_magnitude"))
            if finite(orders):
                executed_candidates.append((orders, row))

        entry: dict[str, Any] = {
            "bias_V": bias,
            "attempted_rows": attempted,
            "executed_rows": sum(1 for row in rows if str(row.get("failure_class")) == "executed"),
            "run_failure_rows": sum(1 for row in rows if str(row.get("failure_class")) == "run_failure"),
            "non_convergence_rows": sum(1 for row in rows if str(row.get("failure_class")) == "non_convergence"),
            "not_executed_rows": sum(1 for row in rows if str(row.get("failure_class")) == "not_executed"),
        }
        if executed_candidates:
            best_orders, best_row = min(executed_candidates, key=lambda item: item[0])
            entry["best_executed_case"] = str(best_row.get("case", ""))
            entry["best_executed_orders_of_magnitude"] = best_orders
        per_bias_summary.append(entry)

    active_bias_summary = [entry for entry in per_bias_summary if int(entry["executed_rows"]) > 0]
    return {
        "bias_values": sorted(set(bias_values)),
        "total_rows": len(matrix_rows),
        "attempted_rows": len(matrix_rows),
        "executed_rows": int(failure_class_counts.get("executed", 0)),
        "failure_class_counts": failure_class_counts,
        "failure_reason_counts": dict(sorted(failure_reason_counts.items())),
        "per_bias_summary": per_bias_summary,
        "active_bias_summary": active_bias_summary,
    }


def extract_fresh_baseline_settings(base_cfg: dict[str, Any]) -> dict[str, str]:
    solver = base_cfg.get("solver", {}) if isinstance(base_cfg.get("solver", {}), dict) else {}
    mobility_model = "default"
    mobility_obj = solver.get("mobility")
    if isinstance(mobility_obj, dict):
        mobility_model = str(mobility_obj.get("model", "default"))

    recomb = solver.get("recombination", ["none"])
    if isinstance(recomb, list):
        recomb_name = "+".join(str(x) for x in recomb)
    else:
        recomb_name = str(recomb)

    bgn = str(solver.get("bandgap_narrowing", "default_inherited"))
    impact_model = "none"
    impact = solver.get("impact_ionization")
    if isinstance(impact, dict):
        impact_model = str(impact.get("model", "none"))

    return {
        "mobility": mobility_model,
        "recombination": recomb_name,
        "bandgap_narrowing": bgn,
        "impact_ionization": impact_model,
    }


def compute_baseline_alignment(
    matrix_rows: list[dict[str, Any]],
    fresh_baseline_settings: dict[str, str],
) -> dict[str, Any]:
    matching_cases: list[str] = []
    mismatch_counts: dict[str, int] = {
        "mobility": 0,
        "recombination": 0,
        "bandgap_narrowing": 0,
        "impact_ionization": 0,
    }

    for row in matrix_rows:
        settings = effective_row_settings(row)
        row_mismatches = [
            key
            for key in mismatch_counts
            if str(settings.get(key, "")) != str(fresh_baseline_settings.get(key, ""))
        ]
        if not row_mismatches:
            matching_cases.append(str(row.get("case", "")))
        else:
            for key in row_mismatches:
                mismatch_counts[key] += 1

    baseline_equivalent_cases = sorted({str(row.get("case", "")) for row in matrix_rows if baseline_equivalent_row(row)})
    implicit_default_case = "m_default_bgn_default_recomb_none_ii_none"
    implicit_default_case_present = any(str(row.get("case", "")) == implicit_default_case for row in matrix_rows)

    return {
        "fresh_baseline_settings": fresh_baseline_settings,
        "matching_case_count": len(set(matching_cases)),
        "matching_cases": sorted(set(matching_cases)),
        "mismatch_counts": mismatch_counts,
        "baseline_equivalent_case_count": len(baseline_equivalent_cases),
        "baseline_equivalent_cases": baseline_equivalent_cases,
        "baseline_equivalent_definition": {
            "mobility": "caughey_thomas",
            "recombination": "srh",
            "bandgap_narrowing": "none_or_default_inherited",
            "impact_ionization": "none",
        },
        "implicit_default_case": implicit_default_case,
        "implicit_default_case_present": implicit_default_case_present,
    }


def run_case(runner: Path, config: Path, cwd: Path, timeout_s: int = 60) -> RunResult:
    try:
        cp = subprocess.run(
            [str(runner), "--config", str(config)],
            cwd=cwd,
            check=True,
            timeout=timeout_s,
            capture_output=True,
            text=True,
        )
        return RunResult(status="ok", reason=(cp.stdout or "").strip()[:400])
    except subprocess.TimeoutExpired:
        return RunResult(status="timeout", reason="runner_timeout")
    except subprocess.CalledProcessError as ex:
        stderr = (ex.stderr or "").strip()
        stdout = (ex.stdout or "").strip()
        reason = stderr if stderr else stdout
        return RunResult(status="failed", reason=reason[:400] if reason else "runner_failed")


def classify_iv_axis(tag: str, contact_polarity: str, local_vbias: float) -> tuple[bool, str]:
    collapsed = bool(tag.startswith("n_only_") and contact_polarity == "n_contact" and abs(local_vbias) < 1.0e-15)
    reason = (
        "Anode-biased sweep grounds Cathode locally; n-contact-only threshold/strength axis collapses"
        if collapsed
        else ""
    )
    return collapsed, reason


def load_mesh_contacts_and_edges(mesh_path: Path) -> tuple[dict[str, list[int]], dict[str, list[tuple[int, int]]]]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    contacts: dict[str, list[int]] = {}
    for c in mesh.get("contacts", []):
        contacts[str(c.get("name", ""))] = [int(n) for n in c.get("node_ids", [])]

    triangles = mesh.get("triangles", [])
    contact_sets = {name: set(ids) for name, ids in contacts.items()}
    edges: dict[str, set[tuple[int, int]]] = {name: set() for name in contacts}

    for tri in triangles:
        ids = [int(n) for n in tri.get("node_ids", [])]
        if len(ids) != 3:
            continue
        pairs = ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0]))
        for a, b in pairs:
            for cname, cnodes in contact_sets.items():
                a_on = a in cnodes
                b_on = b in cnodes
                if a_on == b_on:
                    continue
                cnode, inode = (a, b) if a_on else (b, a)
                edges[cname].add((cnode, inode))

    return contacts, {k: sorted(v) for k, v in edges.items()}


def parse_vtk_scalar_fields(vtk_path: Path) -> dict[str, list[float]]:
    text = vtk_path.read_text(encoding="utf-8", errors="replace")
    out: dict[str, list[float]] = {}
    for field in ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi"]:
        m = re.search(
            rf"SCALARS\s+{re.escape(field)}\s+\w+\s+1\s*\nLOOKUP_TABLE\s+\w+\s*\n(.*?)(?=\nSCALARS\s+|\Z)",
            text,
            flags=re.DOTALL,
        )
        if not m:
            continue
        vals = [float(v) for v in m.group(1).split()]
        out[field] = vals
    return out


def _mean(values: list[float]) -> float:
    if not values:
        return math.nan
    return sum(values) / float(len(values))


def _lookup_bias_metrics(bias_map: dict[float, dict[str, float]], bias: float) -> dict[str, float] | None:
    best_key: float | None = None
    best_err = 1.0e300
    for k in bias_map:
        err = abs(k - bias)
        if err < best_err:
            best_err = err
            best_key = k
    if best_key is None:
        return None
    if best_err > max(1.0e-12, abs(bias) * 1.0e-9):
        return None
    return bias_map[best_key]


def collect_iv_fine_interior_metrics(runner: Path, config_path: Path) -> dict[float, dict[str, float]]:
    if not config_path.exists():
        return {}

    work_dir = config_path.parent
    cfg = json.loads(config_path.read_text(encoding="utf-8"))
    mesh_rel = str(cfg.get("mesh_file", "mesh.json"))
    mesh_path = (work_dir / mesh_rel).resolve() if not Path(mesh_rel).is_absolute() else Path(mesh_rel)
    if not mesh_path.exists():
        return {}

    sweep = cfg.get("sweep", {})
    contact_under_test = str(sweep.get("current_contact", ""))
    if not contact_under_test:
        return {}

    contacts, contact_edges = load_mesh_contacts_and_edges(mesh_path)
    if contact_under_test not in contacts or contact_under_test not in contact_edges:
        return {}

    cfg_mod = json.loads(json.dumps(cfg))
    cfg_mod["sweep"]["write_vtk"] = True
    tmp_csv = f"__tmp_iv_vtk_{config_path.stem}.csv"
    tmp_cfg = f"__tmp_iv_vtk_{config_path.stem}.json"
    cfg_mod["output_csv"] = tmp_csv

    before_vtk = {p.name for p in work_dir.glob("dc_sweep_*.vtk")}
    (work_dir / tmp_cfg).write_text(json.dumps(cfg_mod, indent=2) + "\n", encoding="utf-8")
    run = run_case(runner=runner, config=work_dir / tmp_cfg, cwd=work_dir)
    if run.status != "ok":
        return {}

    after_vtk = sorted(work_dir.glob("dc_sweep_*.vtk"))
    new_vtk = [p for p in after_vtk if p.name not in before_vtk]
    if not new_vtk:
        return {}

    out_rows = read_csv(work_dir / tmp_csv) if (work_dir / tmp_csv).exists() else []
    out_by_idx: dict[int, dict[str, str]] = {i: r for i, r in enumerate(out_rows)}

    bias_map: dict[float, dict[str, float]] = {}
    c_nodes = contacts[contact_under_test]
    c_edges = contact_edges[contact_under_test]

    for vtk in new_vtk:
        m = re.match(r"dc_sweep_(\d+)_.*\.vtk$", vtk.name)
        if not m:
            continue
        idx = int(m.group(1))
        row = out_by_idx.get(idx)
        if row is None:
            continue
        try:
            bias = float(row.get("bias_V", ""))
        except (TypeError, ValueError):
            continue
        fields = parse_vtk_scalar_fields(vtk)
        psi = fields.get("Potential", [])
        phin = fields.get("ElectronQuasiFermi", [])
        phip = fields.get("HoleQuasiFermi", [])
        if not psi or not phin or not phip:
            continue

        psi_c = _mean([psi[i] for i in c_nodes if i < len(psi)])
        phin_c = _mean([phin[i] for i in c_nodes if i < len(phin)])
        phip_c = _mean([phip[i] for i in c_nodes if i < len(phip)])
        dpsi = _mean([psi[j] - psi[i] for i, j in c_edges if i < len(psi) and j < len(psi)])
        dphin = _mean([phin[j] - phin[i] for i, j in c_edges if i < len(phin) and j < len(phin)])
        dphip = _mean([phip[j] - phip[i] for i, j in c_edges if i < len(phip) and j < len(phip)])

        bias_map[bias] = {
            "psi_boundary_target_V": psi_c,
            "phin_boundary_target_V": phin_c,
            "phip_boundary_target_V": phip_c,
            "adjacent_interior_dpsi_V": dpsi,
            "adjacent_interior_dphin_V": dphin,
            "adjacent_interior_dphip_V": dphip,
        }

    for p in new_vtk:
        try:
            p.unlink()
        except OSError:
            pass
    for p in [work_dir / tmp_csv, work_dir / tmp_cfg]:
        if p.exists():
            try:
                p.unlink()
            except OSError:
                pass

    return bias_map


def build_iv_outputs(root: Path, reports_dir: Path, runner: Path) -> dict[str, Any]:
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
                    " threshold/strength does not provide a continuous effective sweep axis"
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
    iv_metrics_cache: dict[tuple[str, str], dict[float, dict[str, float]]] = {}
    for row in contact_summary:
        tag = row["tag"]
        threshold = parse_threshold_from_tag(tag)
        cdir = root / "candidates" / tag
        for side in ["cathode", "anode"]:
            cfg_path = cdir / f"simulation_iv_fine_{side}_{tag}.json"
            iv_metrics_cache[(tag, side)] = collect_iv_fine_interior_metrics(runner=runner, config_path=cfg_path)
            path = cdir / f"pn2d_iv_fine_{side}_{tag}.csv"
            if not path.exists():
                continue
            for curve_row in read_csv(path):
                bias = f(curve_row, "bias_V")
                bias_contact = curve_row.get("bias_contact", "")
                current_contact = curve_row.get("current_contact", "")
                local_vbias = bias if current_contact == bias_contact else 0.0
                contact_polarity = "n_contact" if current_contact.lower() == "cathode" else "p_contact"
                n_only_active = (
                    tag.startswith("n_only_")
                    and contact_polarity == "n_contact"
                    and threshold is not None
                    and local_vbias >= threshold
                )
                collapsed, collapse_reason = classify_iv_axis(tag, contact_polarity, local_vbias)
                metrics = _lookup_bias_metrics(iv_metrics_cache.get((tag, side), {}), bias)
                iv_sem_rows.append(
                    {
                        "candidate_tag": tag,
                        "bias_V": bias,
                        "bias_contact": bias_contact,
                        "contact_under_test": current_contact,
                        "contact_polarity": contact_polarity,
                        "local_Vbias_V": local_vbias,
                        "terminal_voltage_difference_V": bias,
                        "n_only_threshold_V": "" if threshold is None else threshold,
                        "n_only_condition_active": str(bool(n_only_active)).lower(),
                        "degeneracy_flag": str(collapsed).lower(),
                        "degeneracy_reason": collapse_reason,
                        "psi_boundary_target_V": metrics["psi_boundary_target_V"] if metrics else local_vbias,
                        "phin_boundary_target_V": metrics["phin_boundary_target_V"] if metrics else local_vbias,
                        "phip_boundary_target_V": metrics["phip_boundary_target_V"] if metrics else local_vbias,
                        "adjacent_interior_dphin_V": metrics["adjacent_interior_dphin_V"] if metrics else math.nan,
                        "adjacent_interior_dphip_V": metrics["adjacent_interior_dphip_V"] if metrics else math.nan,
                        "adjacent_interior_dpsi_V": metrics["adjacent_interior_dpsi_V"] if metrics else math.nan,
                        "interior_drop_source": "vtk_contact_interior_edge_mean" if metrics else "not_available_in_fine_contact_csv",
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
            "psi_boundary_target_V",
            "phin_boundary_target_V",
            "phip_boundary_target_V",
            "adjacent_interior_dphin_V",
            "adjacent_interior_dphip_V",
            "adjacent_interior_dpsi_V",
            "interior_drop_source",
            "current_total_A_per_um",
            "current_electron_A_per_um",
            "current_hole_A_per_um",
            "current_electron_drift_A_per_um",
            "current_electron_diffusion_A_per_um",
            "current_hole_drift_A_per_um",
            "current_hole_diffusion_A_per_um",
        ],
    )

    return {
        "baseline_summary": baseline_summary,
        "contact_summary": contact_summary,
    }


def resolve_reference_paths(base_cfg: dict[str, Any], ref_vela_dir: Path) -> dict[str, Any]:
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
    return cfg


def build_bv_outputs(
    root: Path,
    reports_dir: Path,
    runner: Path,
    reuse_bv_curves: bool,
    reuse_bv_matrix: bool,
    runner_timeout_s: int,
) -> dict[str, Any]:
    ref_vela_dir = root / "reference" / "vela"
    base_cfg_path = ref_vela_dir / "simulation_bv.json"
    base_cfg = json.loads(base_cfg_path.read_text(encoding="utf-8"))
    ref_bv_curve = read_csv(root / "reference" / "reference_curves" / "pn2d_bv_reference.csv")

    bias_values = [0.00, 0.05, 0.10, 0.20, 0.50]
    mobility_opts = ["default", "caughey_thomas"]
    bgn_opts = ["default", "none", "slotboom"]
    recomb_opts = [["none"], ["srh"], ["srh", "auger"]]
    impact_opts = ["none", "selberherr"]

    case_cfg_dir = reports_dir / "bv_case_configs"
    case_cfg_dir.mkdir(parents=True, exist_ok=True)
    case_csv_dir = reports_dir / "bv_case_curves"
    case_csv_dir.mkdir(parents=True, exist_ok=True)
    run_id = str(int(time.time()))
    matrix_path = reports_dir / "bv_physics_matrix.csv"

    if reuse_bv_matrix:
        if not matrix_path.exists():
            raise RuntimeError(f"Missing matrix CSV for --reuse-bv-matrix: {matrix_path}")
        matrix_rows = read_csv(matrix_path)
        summary = summarize_bv_matrix_rows(matrix_rows=matrix_rows, bias_values=bias_values)
        baseline_alignment = compute_baseline_alignment(
            matrix_rows=matrix_rows,
            fresh_baseline_settings=extract_fresh_baseline_settings(base_cfg),
        )
        (reports_dir / "bv_coverage_summary.json").write_text(
            json.dumps(
                {
                    **summary,
                    "matrix_csv": str(matrix_path),
                    "root": str(root),
                    "reports_dir": str(reports_dir),
                    "baseline_alignment": baseline_alignment,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return {
            "rows": matrix_rows,
            "coverage": summary["failure_class_counts"],
            "summary": summary,
            "baseline_alignment": baseline_alignment,
        }

    matrix_rows: list[dict[str, Any]] = []

    for mobility, bgn, recomb, impact in itertools.product(mobility_opts, bgn_opts, recomb_opts, impact_opts):
        recomb_name = "+".join(recomb)
        case_name = f"m_{mobility}_bgn_{bgn}_recomb_{recomb_name.replace('+', '_')}_ii_{impact}"

        cfg_base = resolve_reference_paths(base_cfg=base_cfg, ref_vela_dir=ref_vela_dir)
        cfg_base.setdefault("solver", {})["method"] = "gummel_newton"
        cfg_base["solver"]["recombination"] = recomb
        if bgn != "default":
            cfg_base["solver"]["bandgap_narrowing"] = bgn
        if mobility == "default":
            cfg_base["solver"].pop("mobility", None)
        else:
            cfg_base.setdefault("solver", {}).setdefault("mobility", {})
            cfg_base["solver"]["mobility"]["model"] = "caughey_thomas"
        cfg_base["solver"].setdefault("impact_ionization", {})
        cfg_base["solver"]["impact_ionization"]["model"] = impact
        cfg_base.setdefault("sweep", {})["contact"] = "Cathode"
        cfg_base["sweep"]["current_contact"] = "Cathode"

        for bias in bias_values:
            cfg_bias = json.loads(json.dumps(cfg_base))
            cfg_bias["sweep"]["start"] = 0.0
            cfg_bias["sweep"]["stop"] = bias
            cfg_bias["sweep"]["step"] = 0.05 if bias <= 0.05 else bias

            bias_tag = f"{bias:.2f}".replace(".", "p")
            cfg_path = case_cfg_dir / f"simulation_bv_{case_name}_{bias_tag}.json"
            out_csv = (case_csv_dir / f"pn2d_bv_{case_name}_{run_id}_{bias_tag}.csv").resolve()

            execution_status = "not_executed"
            execution_reason = "reuse_only_no_curve"
            if reuse_bv_curves:
                matches = sorted(case_csv_dir.glob(f"pn2d_bv_{case_name}_*_{bias_tag}.csv"))
                if matches:
                    out_csv = matches[-1]
                    execution_status = "reused"
                    execution_reason = "reused_existing_curve"
            else:
                cfg_bias["output_csv"] = str(out_csv)
                cfg_path.write_text(json.dumps(cfg_bias, indent=2) + "\n", encoding="utf-8")
                run_result = run_case(runner=runner, config=cfg_path, cwd=ref_vela_dir, timeout_s=runner_timeout_s)
                execution_status = "executed" if run_result.status == "ok" else "run_failure"
                execution_reason = run_result.reason

            failure_class = "not_executed"
            failure_reason = execution_reason
            converged = ""
            handoff_stage = ""
            strict_handoff = ""
            candidate_total = math.nan
            candidate_e = math.nan
            candidate_h = math.nan
            candidate_ed = math.nan
            candidate_ef = math.nan
            candidate_hd = math.nan
            candidate_hf = math.nan
            cancel_ratio = math.nan
            orders = math.nan
            residual = math.nan
            ref_total = interp(ref_bv_curve, "bias_V", "current_total", bias)

            have_curve = out_csv.exists() and out_csv.stat().st_size > 0
            if execution_status == "run_failure":
                failure_class = "run_failure"
            elif not have_curve:
                failure_class = "not_executed"
            else:
                curve = read_csv(out_csv)
                row_at_bias = find_row_at_bias(curve, bias)
                if row_at_bias is not None:
                    converged = row_at_bias.get("converged", "")
                    handoff_stage = row_at_bias.get("handoff_stage", "")
                    try:
                        strict_handoff = str(
                            row_at_bias.get("solver_method", "") == "gummel_newton"
                            and handoff_stage == "newton"
                            and int(row_at_bias.get("newton_iterations", "0")) > 0
                        ).lower()
                    except ValueError:
                        strict_handoff = "false"

                candidate_total = interp(curve, "bias_V", "current_total_A_per_um", bias)
                candidate_e = interp(curve, "bias_V", "current_electron_A_per_um", bias)
                candidate_h = interp(curve, "bias_V", "current_hole_A_per_um", bias)
                candidate_ed = interp(curve, "bias_V", "current_electron_drift_A_per_um", bias)
                candidate_ef = interp(curve, "bias_V", "current_electron_diffusion_A_per_um", bias)
                candidate_hd = interp(curve, "bias_V", "current_hole_drift_A_per_um", bias)
                candidate_hf = interp(curve, "bias_V", "current_hole_diffusion_A_per_um", bias)
                if finite(candidate_total) and finite(ref_total):
                    orders = orders_of_magnitude(ref_total, candidate_total)
                if finite(candidate_ed) and finite(candidate_e):
                    cancel_ratio = abs(candidate_ed) / max(abs(candidate_e), 1.0e-300)
                if finite(candidate_total) and finite(candidate_e) and finite(candidate_h):
                    residual = candidate_total - (candidate_e + candidate_h)

                if (
                    str(converged).lower() in ("0", "false", "no")
                    or not finite(candidate_total)
                    or row_at_bias is None
                ):
                    failure_class = "non_convergence"
                    if row_at_bias is None:
                        failure_reason = "missing_target_bias_row"
                    elif not finite(candidate_total):
                        failure_reason = "nan_or_inf_candidate_total"
                    else:
                        failure_reason = "converged_flag_false"
                else:
                    failure_class = "executed"
                    failure_reason = "ok"

            matrix_rows.append(
                {
                    "case": case_name,
                    "bias_V": bias,
                    "mobility": mobility,
                    "bandgap_narrowing": bgn,
                    "recombination": recomb_name,
                    "impact_ionization": impact,
                    "execution_status": execution_status,
                    "failure_class": failure_class,
                    "failure_reason": normalize_run_failure_reason(failure_reason) if failure_class == "run_failure" else failure_reason,
                    "converged": converged,
                    "handoff_stage": handoff_stage,
                    "strict_newton_handoff": strict_handoff,
                    "effective_mobility_model": "default_inherited" if mobility == "default" else "caughey_thomas",
                    "effective_recombination_models": recomb_name,
                    "effective_bandgap_narrowing": "default_inherited" if bgn == "default" else bgn,
                    "effective_impact_ionization_model": impact,
                    "reference_total_A_per_um": ref_total,
                    "candidate_total_A_per_um": candidate_total,
                    "candidate_electron_A_per_um": candidate_e,
                    "candidate_hole_A_per_um": candidate_h,
                    "electron_drift_A_per_um": candidate_ed,
                    "electron_diffusion_A_per_um": candidate_ef,
                    "hole_drift_A_per_um": candidate_hd,
                    "hole_diffusion_A_per_um": candidate_hf,
                    "cancel_ratio_abs_edrift_over_etotal": cancel_ratio,
                    "orders_of_magnitude": orders,
                    "terminal_sum_residual_A_per_um": residual,
                    "config": str(cfg_path),
                    "csv_file": str(out_csv),
                }
            )

    write_csv(
        reports_dir / "bv_physics_matrix.csv",
        matrix_rows,
        [
            "case",
            "bias_V",
            "mobility",
            "bandgap_narrowing",
            "recombination",
            "impact_ionization",
            "execution_status",
            "failure_class",
            "failure_reason",
            "converged",
            "handoff_stage",
            "strict_newton_handoff",
            "effective_mobility_model",
            "effective_recombination_models",
            "effective_bandgap_narrowing",
            "effective_impact_ionization_model",
            "reference_total_A_per_um",
            "candidate_total_A_per_um",
            "candidate_electron_A_per_um",
            "candidate_hole_A_per_um",
            "electron_drift_A_per_um",
            "electron_diffusion_A_per_um",
            "hole_drift_A_per_um",
            "hole_diffusion_A_per_um",
            "cancel_ratio_abs_edrift_over_etotal",
            "orders_of_magnitude",
            "terminal_sum_residual_A_per_um",
            "config",
            "csv_file",
        ],
    )

    summary = summarize_bv_matrix_rows(matrix_rows=matrix_rows, bias_values=bias_values)
    baseline_alignment = compute_baseline_alignment(
        matrix_rows=matrix_rows,
        fresh_baseline_settings=extract_fresh_baseline_settings(base_cfg),
    )

    (reports_dir / "bv_coverage_summary.json").write_text(
        json.dumps(
            {
                **summary,
                "matrix_csv": str(matrix_path),
                "root": str(root),
                "reports_dir": str(reports_dir),
                "baseline_alignment": baseline_alignment,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "rows": matrix_rows,
        "coverage": summary["failure_class_counts"],
        "summary": summary,
        "baseline_alignment": baseline_alignment,
    }


def write_rootcause_md(
    root: Path,
    reports_dir: Path,
    baseline_summary: dict[str, Any],
    contact_summary: list[dict[str, str]],
    bv_rows: list[dict[str, Any]],
    bv_summary: dict[str, Any],
    baseline_alignment: dict[str, Any],
) -> None:
    iv_orders = baseline_summary["checks"]["iv_orders_0p2_to_0p3"]
    terminal_sum = baseline_summary["checks"]["terminal_sum_abs_A_per_um_at_0p3"]
    strict_handoff = baseline_summary["checks"]["strict_newton_all"]

    best_bv_005 = math.nan
    for row in bv_rows:
        if abs(float(row["bias_V"]) - 0.05) < 1.0e-12 and row["failure_class"] == "executed":
            v = float(row["orders_of_magnitude"])
            if finite(v):
                best_bv_005 = v if not finite(best_bv_005) else min(best_bv_005, v)

    generalizable: list[str] = []
    diagnostic_only: list[str] = []

    for row in bv_rows:
        if abs(float(row["bias_V"]) - 0.05) > 1.0e-12:
            continue
        case = str(row["case"])
        failure = str(row["failure_class"])
        orders = float(row["orders_of_magnitude"]) if finite(float(row["orders_of_magnitude"])) else math.nan
        if failure == "executed" and finite(orders) and orders < 0.15:
            generalizable.append(f"{case} (orders@0.05V={orders:.4f})")
        else:
            diagnostic_only.append(f"{case} ({failure})")

    n_only_collapsed = any(r.get("tag", "").startswith("n_only_") for r in contact_summary)

    interior_drop_rows: list[dict[str, str]] = []
    sem_path = reports_dir / "iv_contact_boundary_semantics.csv"
    if sem_path.exists():
        interior_drop_rows = read_csv(sem_path)
    has_vtk_interior_drop = any(
        str(row.get("interior_drop_source", "")) == "vtk_contact_interior_edge_mean"
        for row in interior_drop_rows
    )

    carrier_recon_line = (
        "- carrier reconstruction: interior QF drop probes are available in iv_contact_boundary_semantics.csv (interior_drop_source=vtk_contact_interior_edge_mean)."
        if has_vtk_interior_drop
        else "- carrier reconstruction: requires additional interior QF drop probes (not present in fine contact CSV)."
    )

    lines = [
        "# pn2d IV/BV Root-Cause Next",
        "",
        "## Fresh Artifacts Only",
        f"- Source root: {root}",
        f"- Reports dir: {reports_dir}",
        "",
        "## Baseline Checks",
        f"- IV orders(0.2-0.3V): {iv_orders:.6f}",
        f"- Terminal sum |I| at 0.3V (A/um): {terminal_sum:.3e}",
        f"- Strict Newton handoff all: {strict_handoff}",
        f"- Best BV orders@0.05V among executed rows: {best_bv_005:.6f}" if finite(best_bv_005) else "- Best BV orders@0.05V among executed rows: nan",
        "",
        "## IV Semantic Localization",
        "- boundary target semantics: confirmed as active root-cause axis (n-contact-only axis collapse under anode-driven sweep)." if n_only_collapsed else "- boundary target semantics: not confirmed in this batch.",
        carrier_recon_line,
        "- current extraction: high drift/diff cancellation suggests extraction sensitivity remains relevant.",
        "- mobility/SG second-order error: retained as coupled hypothesis with current extraction.",
        "",
        "## Candidate Classification",
        "### Generalizable Fix Candidates",
    ]

    if generalizable:
        lines.extend([f"- {item}" for item in sorted(set(generalizable))])
    else:
        lines.append("- none (no executed row at 0.05V reached orders < 0.15)")

    lines.extend([
        "",
        "### Diagnostic-Only / Non-Generalizable Candidates",
    ])
    if diagnostic_only:
        lines.extend([f"- {item}" for item in sorted(set(diagnostic_only))])
    else:
        lines.append("- none")

    lines.extend([
        "",
        "## BV Coverage Summary",
        f"- attempted_rows: {bv_summary.get('attempted_rows', 0)}",
        f"- executed_rows: {bv_summary.get('executed_rows', 0)}",
        f"- failure_class_counts: {json.dumps(bv_summary.get('failure_class_counts', {}), sort_keys=True)}",
        f"- failure_reason_counts: {json.dumps(bv_summary.get('failure_reason_counts', {}), sort_keys=True)}",
        "",
        "### Active Bias Summary",
    ])
    active_bias_summary = list(bv_summary.get("active_bias_summary", []))
    if active_bias_summary:
        lines.append("| bias_V | attempted_rows | executed_rows | run_failure_rows | non_convergence_rows | best_executed_case | best_executed_orders |")
        lines.append("|---:|---:|---:|---:|---:|---|---:|")
        for entry in active_bias_summary:
            best_case = str(entry.get("best_executed_case", ""))
            best_orders = _to_float(entry.get("best_executed_orders_of_magnitude"))
            best_orders_text = f"{best_orders:.6f}" if finite(best_orders) else ""
            lines.append(
                "| "
                f"{_to_float(entry.get('bias_V')):.2f} | "
                f"{int(entry.get('attempted_rows', 0))} | "
                f"{int(entry.get('executed_rows', 0))} | "
                f"{int(entry.get('run_failure_rows', 0))} | "
                f"{int(entry.get('non_convergence_rows', 0))} | "
                f"{best_case} | "
                f"{best_orders_text} |"
            )
    else:
        lines.append("- none (no bias has executed rows)")

    lines.extend([
        "",
        "## Baseline Alignment Notes",
        "- Fresh baseline settings: "
        + json.dumps(baseline_alignment.get("fresh_baseline_settings", {}), sort_keys=True),
        "- Matrix rows exactly matching fresh baseline settings: "
        + str(baseline_alignment.get("matching_case_count", 0)),
        "- Baseline-equivalent physics definition used here: CT mobility + SRH recombination + BGN none/default-inherited + impact none.",
        "- Baseline-equivalent cases in matrix: "
        + json.dumps(baseline_alignment.get("baseline_equivalent_cases", [])),
        "- Implicit-default matrix case is kept as a separate matrix point (not baseline-equivalent): "
        + str(baseline_alignment.get("implicit_default_case", "m_default_bgn_default_recomb_none_ii_none")),
    ])

    (reports_dir / "pn2d_iv_bv_rootcause_next.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_artifact_manifest(
    reports_dir: Path,
    root: Path,
    command_line: str,
) -> None:
    artifacts: list[dict[str, Any]] = []
    source_decks = [
        str(root / "reference" / "vela" / "simulation_iv.json"),
        str(root / "reference" / "vela" / "simulation_bv.json"),
    ]
    for path in sorted(reports_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() not in {".csv", ".json", ".md"}:
            continue
        artifacts.append(
            {
                "path": str(path),
                "generated_by_command": command_line,
                "source_decks": source_decks,
            }
        )

    manifest = {
        "generated_at_epoch_s": int(time.time()),
        "root": str(root),
        "reports_dir": str(reports_dir),
        "artifacts": artifacts,
    }
    (reports_dir / "artifact_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )


def build_reports(
    root: Path,
    runner: Path,
    reuse_bv_curves: bool = False,
    reuse_bv_matrix: bool = False,
    command_line: str = "",
    runner_timeout_s: int = 60,
) -> None:
    reports_dir = root / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    iv_info = build_iv_outputs(root=root, reports_dir=reports_dir, runner=runner)
    bv_info = build_bv_outputs(
        root=root,
        reports_dir=reports_dir,
        runner=runner,
        reuse_bv_curves=reuse_bv_curves,
        reuse_bv_matrix=reuse_bv_matrix,
        runner_timeout_s=runner_timeout_s,
    )
    write_rootcause_md(
        root=root,
        reports_dir=reports_dir,
        baseline_summary=iv_info["baseline_summary"],
        contact_summary=iv_info["contact_summary"],
        bv_rows=bv_info["rows"],
        bv_summary=bv_info.get("summary", {}),
        baseline_alignment=bv_info.get("baseline_alignment", {}),
    )
    write_artifact_manifest(
        reports_dir=reports_dir,
        root=root,
        command_line=command_line,
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
    parser.add_argument(
        "--reuse-bv-matrix",
        action="store_true",
        help="Regenerate BV coverage/report artifacts from existing bv_physics_matrix.csv only.",
    )
    parser.add_argument(
        "--runner-timeout-s",
        type=int,
        default=60,
        help="Timeout in seconds for each runner invocation.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cmd = (
        "python scripts/generate_pn2d_rootcause_reports.py"
        f" --root {args.root.resolve()}"
        f" --runner {args.runner.resolve()}"
        + (" --reuse-bv-curves" if args.reuse_bv_curves else "")
        + (" --reuse-bv-matrix" if args.reuse_bv_matrix else "")
        + f" --runner-timeout-s {max(1, int(args.runner_timeout_s))}"
    )
    build_reports(
        root=args.root.resolve(),
        runner=args.runner.resolve(),
        reuse_bv_curves=bool(args.reuse_bv_curves),
        reuse_bv_matrix=bool(args.reuse_bv_matrix),
        command_line=cmd,
        runner_timeout_s=max(1, int(args.runner_timeout_s)),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())