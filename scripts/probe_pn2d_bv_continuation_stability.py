#!/usr/bin/env python3
"""Focused pn2d BV continuation stability probe for impact=none candidates.

This script is diagnostic/reporting only. It does not change solver code or
default deck behavior.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]

TARGET_CASES = [
    "m_default_bgn_default_recomb_srh_auger_ii_none",
    "m_default_bgn_none_recomb_srh_auger_ii_none",
    "m_caughey_thomas_bgn_default_recomb_srh_ii_none",
    "m_caughey_thomas_bgn_none_recomb_srh_ii_none",
    "m_default_bgn_slotboom_recomb_srh_auger_ii_none",
]

VARIANTS: list[tuple[float, int]] = [
    (0.05, 180),
    (0.025, 180),
    (0.01, 300),
]

ORDER_BIASES = [0.05, 0.10, 0.20, 0.50]


def read_csv(path: Path) -> list[dict[str, str]]:
    text = path.read_text(encoding="utf-8", errors="replace").replace("\x00", "")
    return list(csv.DictReader(io.StringIO(text)))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def finite(value: float) -> bool:
    return not (math.isnan(value) or math.isinf(value))


def find_row_at_bias(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    best: dict[str, str] | None = None
    best_err = 1.0e300
    for row in rows:
        b = to_float(row.get("bias_V", ""))
        if not finite(b):
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


def strict_newton(rows: list[dict[str, str]]) -> bool:
    saw_any = False
    for row in rows:
        saw_any = True
        if row.get("solver_method") != "gummel_newton":
            return False
        if row.get("handoff_stage") != "newton":
            return False
        if int(row.get("newton_iterations", "0")) <= 0:
            return False
    return saw_any


def orders_of_magnitude(reference: float, candidate: float) -> float:
    eps = 1.0e-300
    return abs(math.log10(max(abs(candidate), eps) / max(abs(reference), eps)))


def interp(rows: list[dict[str, str]], x_key: str, y_key: str, x: float) -> float:
    pts: list[tuple[float, float]] = []
    for row in rows:
        xv = to_float(row.get(x_key, ""))
        yv = to_float(row.get(y_key, ""))
        if not (finite(xv) and finite(yv)):
            continue
        pts.append((xv, yv))
    pts.sort(key=lambda p: p[0])
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


def normalize_run_failure_reason(reason: str) -> str:
    if reason == "runner_timeout":
        return "runner_timeout"
    return "runner_nonzero_exit_with_json"


@dataclass
class RunResult:
    status: str
    reason: str


def run_case(runner: Path, config: Path, cwd: Path, timeout_s: int) -> RunResult:
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


def parse_case(case_name: str) -> tuple[str, str, list[str], str]:
    m = re.match(r"^m_(.+?)_bgn_(.+?)_recomb_(.+?)_ii_(.+)$", case_name)
    if not m:
        raise ValueError(f"Unrecognized case format: {case_name}")
    mobility = m.group(1)
    bgn = m.group(2)
    recomb = m.group(3).split("_")
    impact = m.group(4)
    return mobility, bgn, recomb, impact


def build_case_config(base_cfg: dict[str, Any], ref_vela_dir: Path, case_name: str, stop: float, step: float, out_csv: Path) -> dict[str, Any]:
    mobility, bgn, recomb, impact = parse_case(case_name)
    cfg = resolve_reference_paths(base_cfg=base_cfg, ref_vela_dir=ref_vela_dir)
    cfg.setdefault("solver", {})["method"] = "gummel_newton"
    cfg["solver"]["recombination"] = recomb

    if bgn == "default":
        cfg["solver"].pop("bandgap_narrowing", None)
    else:
        cfg["solver"]["bandgap_narrowing"] = bgn

    if mobility == "default":
        cfg["solver"].pop("mobility", None)
    else:
        cfg["solver"].setdefault("mobility", {})
        cfg["solver"]["mobility"]["model"] = "caughey_thomas"

    cfg["solver"].setdefault("impact_ionization", {})
    cfg["solver"]["impact_ionization"]["model"] = impact

    cfg.setdefault("sweep", {})["mode"] = "bv_reverse"
    cfg["sweep"]["start"] = 0.0
    cfg["sweep"]["stop"] = stop
    cfg["sweep"]["step"] = step
    cfg["sweep"]["contact"] = "Cathode"
    cfg["sweep"]["current_contact"] = "Cathode"
    cfg["output_csv"] = str(out_csv.resolve())
    return cfg


def classify_run(
    run: RunResult,
    curve_rows: list[dict[str, str]],
    stop: float,
    ref_curve: list[dict[str, str]],
) -> tuple[str, str, int, float, str, str, dict[float, str], str]:
    if run.status != "ok":
        return (
            "run_failure",
            normalize_run_failure_reason(run.reason),
            0,
            math.nan,
            "",
            "",
            {b: "" for b in ORDER_BIASES},
            "false",
        )

    rows_sorted = sorted(
        [r for r in curve_rows if finite(to_float(r.get("bias_V", "")))],
        key=lambda r: to_float(r.get("bias_V", "")),
    )
    if not rows_sorted:
        return (
            "non_convergence",
            "missing_curve_rows",
            0,
            math.nan,
            "",
            "",
            {b: "" for b in ORDER_BIASES},
            "false",
        )

    points = len(rows_sorted)
    last = rows_sorted[-1]
    last_bias = to_float(last.get("bias_V", ""))
    last_converged = str(last.get("converged", ""))
    strict_handoff = str(strict_newton(rows_sorted)).lower()

    target = find_row_at_bias(rows_sorted, stop)
    candidate_total = math.nan
    converged = ""
    if target is not None:
        candidate_total = to_float(target.get("current_total_A_per_um", ""))
        converged = str(target.get("converged", ""))

    if target is None:
        status = "non_convergence"
        reason = "missing_target_bias_row"
    elif not finite(candidate_total):
        status = "non_convergence"
        reason = "nan_or_inf_candidate_total"
    elif converged.lower() in {"0", "false", "no"}:
        status = "non_convergence"
        reason = "converged_flag_false"
    else:
        status = "executed"
        reason = "ok"

    orders_out: dict[float, str] = {}
    for bias in ORDER_BIASES:
        row = find_row_at_bias(rows_sorted, bias)
        if row is None:
            orders_out[bias] = ""
            continue
        if str(row.get("converged", "")).lower() in {"0", "false", "no"}:
            orders_out[bias] = ""
            continue
        cand = to_float(row.get("current_total_A_per_um", ""))
        ref = interp(ref_curve, "bias_V", "current_total", bias)
        if not (finite(cand) and finite(ref)):
            orders_out[bias] = ""
            continue
        orders_out[bias] = f"{orders_of_magnitude(ref, cand):.6f}"

    return (
        status,
        reason,
        points,
        last_bias,
        last_converged,
        converged,
        orders_out,
        strict_handoff,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Probe BV continuation stability at 0.50V for selected impact=none cases.")
    parser.add_argument(
        "--root",
        type=Path,
        default=REPO / "build" / "pn2d_curve_rootcause_20260527_round3",
        help="Root directory containing reference/ and reports/ artifacts.",
    )
    parser.add_argument(
        "--runner",
        type=Path,
        default=REPO / "build" / "vela_example_runner.exe",
        help="Path to vela_example_runner executable.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    reports_dir = root / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    ref_vela_dir = root / "reference" / "vela"
    base_cfg_path = ref_vela_dir / "simulation_bv.json"
    ref_curve_path = root / "reference" / "reference_curves" / "pn2d_bv_reference.csv"

    if not args.runner.exists():
        raise FileNotFoundError(f"runner not found: {args.runner}")
    if not base_cfg_path.exists():
        raise FileNotFoundError(f"missing base BV config: {base_cfg_path}")
    if not ref_curve_path.exists():
        raise FileNotFoundError(f"missing BV reference curve: {ref_curve_path}")

    base_cfg = json.loads(base_cfg_path.read_text(encoding="utf-8"))
    ref_curve = read_csv(ref_curve_path)

    cfg_dir = reports_dir / "bv_continuation_configs"
    csv_dir = reports_dir / "bv_continuation_curves"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    csv_dir.mkdir(parents=True, exist_ok=True)

    out_rows: list[dict[str, Any]] = []

    for case_name in TARGET_CASES:
        for step, timeout_s in VARIANTS:
            print(
                f"[probe] case={case_name} stop=0.50 step={step} timeout_s={timeout_s}",
                flush=True,
            )
            step_tag = f"{step:.3f}".rstrip("0").rstrip(".").replace(".", "p")
            cfg_path = cfg_dir / f"simulation_bv_{case_name}_stop0p50_step{step_tag}.json"
            out_csv = csv_dir / f"pn2d_bv_{case_name}_stop0p50_step{step_tag}.csv"

            cfg = build_case_config(
                base_cfg=base_cfg,
                ref_vela_dir=ref_vela_dir,
                case_name=case_name,
                stop=0.50,
                step=step,
                out_csv=out_csv,
            )
            cfg_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")

            run = run_case(
                runner=args.runner.resolve(),
                config=cfg_path.resolve(),
                cwd=ref_vela_dir,
                timeout_s=timeout_s,
            )
            curve_rows = read_csv(out_csv) if out_csv.exists() and out_csv.stat().st_size > 0 else []

            (
                status,
                reason,
                points,
                last_bias,
                last_converged,
                converged_at_stop,
                orders,
                strict_handoff,
            ) = classify_run(
                run=run,
                curve_rows=curve_rows,
                stop=0.50,
                ref_curve=ref_curve,
            )

            out_rows.append(
                {
                    "case": case_name,
                    "stop": f"{0.50:.2f}",
                    "step": f"{step:.3f}".rstrip("0").rstrip("."),
                    "timeout": str(timeout_s),
                    "status": status,
                    "failure_reason_class": reason,
                    "points": points,
                    "last_bias": "" if not finite(last_bias) else f"{last_bias:.12g}",
                    "last_converged": last_converged,
                    "converged_at_stop": converged_at_stop,
                    "orders_at_0p05": orders[0.05],
                    "orders_at_0p10": orders[0.10],
                    "orders_at_0p20": orders[0.20],
                    "orders_at_0p50": orders[0.50],
                    "strict_handoff": strict_handoff,
                    "csv_path": str(out_csv.resolve()),
                    "config_path": str(cfg_path.resolve()),
                }
            )
            print(
                f"[done] case={case_name} step={step} status={status} reason={reason} points={points}",
                flush=True,
            )

    out_path = reports_dir / "bv_continuation_stability.csv"
    write_csv(
        out_path,
        out_rows,
        [
            "case",
            "stop",
            "step",
            "timeout",
            "status",
            "failure_reason_class",
            "points",
            "last_bias",
            "last_converged",
            "converged_at_stop",
            "orders_at_0p05",
            "orders_at_0p10",
            "orders_at_0p20",
            "orders_at_0p50",
            "strict_handoff",
            "csv_path",
            "config_path",
        ],
    )

    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())