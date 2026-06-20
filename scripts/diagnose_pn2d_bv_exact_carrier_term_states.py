#!/usr/bin/env python3
"""Prepare and run exact C++ carrier-term probes for focused BV state variants."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


STATE_FIELDS = [
    "ElectrostaticPotential",
    "eQuasiFermiPotential",
    "hQuasiFermiPotential",
]

VARIANTS = [
    "vela_baseline",
    "vela_psi_sentaurus_qf",
    "vela_qf_shift",
]

NODE_FIELDS = [
    "variant",
    "bias_V",
    "node_id",
    "x",
    "y",
    "support_class",
    "electron_flux",
    "electron_recombination",
    "electron_impact",
    "electron_gauge",
    "electron_boundary",
    "electron_term_sum",
    "electron_residual",
    "electron_adjusted_impact",
    "electron_adjusted_term_sum",
    "electron_adjusted_residual",
    "hole_flux",
    "hole_recombination",
    "hole_impact",
    "hole_gauge",
    "hole_boundary",
    "hole_term_sum",
    "hole_residual",
    "hole_adjusted_impact",
    "hole_adjusted_term_sum",
    "hole_adjusted_residual",
    "impact_electron_source",
    "impact_hole_source",
    "impact_combined_source",
    "electron_residual_over_abs_impact",
    "hole_residual_over_abs_impact",
    "electron_required_impact_multiplier",
    "hole_required_impact_multiplier",
    "electron_impact_multiplier_delta",
    "hole_impact_multiplier_delta",
    "electron_closed_residual",
    "hole_closed_residual",
    "donors_m3",
    "acceptors_m3",
    "net_doping_m3",
    "ni_eff_m3",
    "psi_V",
    "phin_V",
    "phip_V",
    "node_volume_m2",
    "electron_density_m3",
    "hole_density_m3",
    "electron_density_over_ni",
    "hole_density_over_ni",
    "srh_excess_product_m6",
    "srh_excess_over_ni2",
    "srh_taup_n_plus_ni_m3_s",
    "srh_taun_p_plus_ni_m3_s",
    "srh_denominator_m3_s",
    "srh_rate_m3_s",
    "srh_recombination_integral_s_inv",
    "srh_recombination_integral_scaled_s_inv",
    "srh_scaled_over_exact_recombination",
    "exact_recombination_over_abs_impact",
    "inferred_continuity_scale_from_srh",
    "block_psi",
    "block_phin",
    "block_phip",
    "block_combined",
    "adjusted_block_psi",
    "adjusted_block_phin",
    "adjusted_block_phip",
    "adjusted_block_combined",
]

SUMMARY_KEYS = [
    "electron_flux",
    "electron_recombination",
    "electron_impact",
    "electron_residual",
    "electron_adjusted_residual",
    "hole_flux",
    "hole_recombination",
    "hole_impact",
    "hole_residual",
    "hole_adjusted_residual",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--electron-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--hole-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--electron-impact-scale", type=float, default=1.0)
    parser.add_argument("--hole-impact-scale", type=float, default=1.0)
    parser.add_argument("--preserve-contact-qf-on-shift", action="store_true")
    parser.add_argument("--bias-contact", default="Anode")
    parser.add_argument("--ground-contact", default="Cathode")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_scalar_csv(path: Path) -> list[float]:
    values_by_id: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            values_by_id[int(row["node_id"])] = float(row["component0"])
    if not values_by_id:
        return []
    return [values_by_id[index] for index in range(max(values_by_id) + 1)]


def write_scalar_csv(path: Path, values: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "component0"])
        for node_id, value in enumerate(values):
            writer.writerow([node_id, f"{value:.17g}"])


def shifted(values: list[float], delta: float, preserve_nodes: set[int] | None = None) -> list[float]:
    fixed = preserve_nodes or set()
    return [value if index in fixed else value + delta for index, value in enumerate(values)]


def load_vela_state(vtk: Path, node_count: int) -> dict[str, list[float]]:
    state = fluxforms.load_vela_state(vtk, node_count)
    return {
        "ElectrostaticPotential": list(state["psi"]),
        "eQuasiFermiPotential": list(state["phin"]),
        "hQuasiFermiPotential": list(state["phip"]),
    }


def load_sentaurus_state(root: Path, node_count: int) -> dict[str, list[float]]:
    result: dict[str, list[float]] = {}
    for field in STATE_FIELDS:
        values = read_scalar_csv(root / "fields" / f"{field}_region0.csv")
        if len(values) != node_count:
            raise RuntimeError(
                f"{field} has {len(values)} values, expected {node_count}")
        result[field] = values
    return result


def write_state(fields_dir: Path, state: dict[str, list[float]]) -> None:
    for field in STATE_FIELDS:
        write_scalar_csv(fields_dir / f"{field}_region0.csv", state[field])


def variant_states(
    vela: dict[str, list[float]],
    sentaurus: dict[str, list[float]],
    electron_shift: float,
    hole_shift: float,
    preserve_shift_nodes: set[int] | None = None,
) -> dict[str, dict[str, list[float]]]:
    return {
        "vela_baseline": {
            "ElectrostaticPotential": list(vela["ElectrostaticPotential"]),
            "eQuasiFermiPotential": list(vela["eQuasiFermiPotential"]),
            "hQuasiFermiPotential": list(vela["hQuasiFermiPotential"]),
        },
        "vela_psi_sentaurus_qf": {
            "ElectrostaticPotential": list(vela["ElectrostaticPotential"]),
            "eQuasiFermiPotential": list(sentaurus["eQuasiFermiPotential"]),
            "hQuasiFermiPotential": list(sentaurus["hQuasiFermiPotential"]),
        },
        "vela_qf_shift": {
            "ElectrostaticPotential": list(vela["ElectrostaticPotential"]),
            "eQuasiFermiPotential": shifted(
                vela["eQuasiFermiPotential"], electron_shift, preserve_shift_nodes),
            "hQuasiFermiPotential": shifted(
                vela["hQuasiFermiPotential"], hole_shift, preserve_shift_nodes),
        },
    }


def support_rows(path: Path) -> list[dict[str, str]]:
    return [
        row for row in read_csv_rows(path)
        if row.get("support_class") not in (None, "", "inactive")
    ]


def resolve_config_path(base_config: Path, raw: str) -> str:
    path = Path(raw)
    return str((path if path.is_absolute() else base_config.parent / path).resolve())


def build_probe_config(
    base_config: Path,
    output_csv: Path,
    fields_dir: Path,
    bias: float,
    bias_contact: str,
    ground_contact: str,
    electron_impact_scale: float,
    hole_impact_scale: float,
) -> dict[str, Any]:
    cfg = json.loads(base_config.read_text())
    cfg["simulation_type"] = "newton_carrier_term_probe"
    cfg["output_csv"] = str(output_csv.resolve())
    cfg["state_fields_dir"] = str(fields_dir.resolve())
    for key in ["mesh_file", "node_doping_file", "materials_file"]:
        if key in cfg:
            cfg[key] = resolve_config_path(base_config, cfg[key])
    cfg.pop("sweep", None)
    cfg["carrier_term_probe"] = {
        "electron_impact_scale": electron_impact_scale,
        "hole_impact_scale": hole_impact_scale,
    }
    for contact in cfg.get("contacts", []):
        name = contact.get("name")
        if name == bias_contact:
            contact["bias"] = bias
        elif name == ground_contact:
            contact["bias"] = 0.0
    return cfg


def run_runner(runner: Path, config: Path) -> dict[str, Any]:
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=config.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {config}: {result.stderr.strip() or result.stdout.strip()}")
    return json.loads(result.stdout)


def read_probe_rows(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row for row in read_csv_rows(path)}


def finite_float(raw: Any) -> float | None:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values = []
    for row in rows:
        value = finite_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def limited_exp(value: float) -> float:
    return math.exp(min(500.0, max(-500.0, value)))


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None:
        return None
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return None
    if abs(denominator) < 1.0e-300:
        return None
    return numerator / denominator


def srh_decomposition(
    *,
    psi: float,
    phin: float,
    phip: float,
    ni_eff: float,
    vt: float,
    node_volume: float,
    taun: float,
    taup: float,
    continuity_scale: float = 1.0,
    exact_recombination: float | None = None,
    exact_impact: float | None = None,
) -> dict[str, float | None]:
    n = ni_eff * limited_exp((psi - phin) / vt)
    p = ni_eff * limited_exp((phip - psi) / vt)
    excess = ni_eff * ni_eff * math.expm1((phip - phin) / vt)
    taup_term = taup * (n + ni_eff)
    taun_term = taun * (p + ni_eff)
    denominator = taup_term + taun_term
    rate = 0.0 if abs(denominator) < 1.0e-100 else excess / denominator
    integral = rate * node_volume
    scaled_integral = integral / continuity_scale
    exact_over_abs_impact = safe_ratio(exact_recombination, abs(exact_impact or 0.0))
    return {
        "electron_density_m3": n,
        "hole_density_m3": p,
        "electron_density_over_ni": safe_ratio(n, ni_eff),
        "hole_density_over_ni": safe_ratio(p, ni_eff),
        "srh_excess_product_m6": excess,
        "srh_excess_over_ni2": safe_ratio(excess, ni_eff * ni_eff),
        "srh_taup_n_plus_ni_m3_s": taup_term,
        "srh_taun_p_plus_ni_m3_s": taun_term,
        "srh_denominator_m3_s": denominator,
        "srh_rate_m3_s": rate,
        "srh_recombination_integral_s_inv": integral,
        "srh_recombination_integral_scaled_s_inv": scaled_integral,
        "srh_scaled_over_exact_recombination": safe_ratio(scaled_integral, exact_recombination),
        "exact_recombination_over_abs_impact": exact_over_abs_impact,
        "inferred_continuity_scale_from_srh": safe_ratio(integral, exact_recombination),
    }


def impact_closure(
    *,
    flux: float,
    recombination: float,
    impact: float,
    residual: float,
) -> dict[str, float | None]:
    required = None
    if math.isfinite(impact) and abs(impact) >= 1.0e-300:
        required = -(flux + recombination) / impact
    closed = None
    if required is not None:
        closed = flux + recombination + required * impact
    return {
        "required_impact_multiplier": required,
        "impact_multiplier_delta": required - 1.0 if required is not None else None,
        "residual_over_abs_impact": safe_ratio(residual, abs(impact)),
        "closed_residual": closed,
    }


def prefixed(prefix: str, values: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {"count": len(rows)}
    for key in SUMMARY_KEYS:
        values = finite_values(rows, key)
        summary[f"{key}_median"] = statistics.median(values) if values else None
        summary[f"{key}_mean"] = statistics.fmean(values) if values else None
    return summary


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for support_class in sorted({str(row["support_class"]) for row in rows}):
        class_rows = [row for row in rows if str(row["support_class"]) == support_class]
        variants: dict[str, Any] = {}
        for variant in sorted({str(row["variant"]) for row in class_rows}):
            variants[variant] = summarize_rows([
                row for row in class_rows if str(row["variant"]) == variant
            ])
        result[support_class] = {"variants": variants}
    return result


def solver_config(base_cfg: dict[str, Any]) -> dict[str, Any]:
    solver = base_cfg.get("solver", {})
    return solver if isinstance(solver, dict) else {}


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    if args.runner is None and not args.prepare_only:
        raise SystemExit("--runner is required unless --prepare-only is set")

    base_cfg = json.loads(args.base_config.read_text())
    cfg_solver = solver_config(base_cfg)
    temperature_k = float(cfg_solver.get("temperature_K", 300.0))
    vt = fluxforms.sgdiag.K_B_OVER_Q * temperature_k
    taun = float(cfg_solver.get("taun", 1.0e-5))
    taup = float(cfg_solver.get("taup", 3.0e-6))
    mesh_path = Path(resolve_config_path(args.base_config, base_cfg["mesh_file"]))
    nodes, triangles, contact_by_node = fluxforms.sgdiag.read_mesh(mesh_path)
    node_count = len(nodes)
    node_volumes = fluxforms.sgdiag.node_volumes(nodes, triangles)
    preserve_shift_nodes = set(contact_by_node) if args.preserve_contact_qf_on_shift else set()
    vtk = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    vela = load_vela_state(vtk, node_count)
    sentaurus = load_sentaurus_state(args.sentaurus_dir, node_count)
    states = variant_states(
        vela,
        sentaurus,
        args.electron_qf_shift_v,
        args.hole_qf_shift_v,
        preserve_shift_nodes,
    )
    supports = support_rows(args.support_csv)
    support_by_node = {int(row["node_id"]): row for row in supports}
    node_rows: list[dict[str, Any]] = []
    state_summaries: list[dict[str, Any]] = []

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for variant in VARIANTS:
        fields_dir = args.out_dir / "states" / variant / "fields"
        output_csv = args.out_dir / "carrier_terms" / f"{variant}_carrier_terms.csv"
        config_path = args.out_dir / "configs" / f"{variant}_newton_carrier_term_probe.json"
        fields_dir.mkdir(parents=True, exist_ok=True)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        write_state(fields_dir, states[variant])
        cfg = build_probe_config(
            args.base_config,
            output_csv,
            fields_dir,
            args.bias,
            args.bias_contact,
            args.ground_contact,
            args.electron_impact_scale,
            args.hole_impact_scale,
        )
        config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
        item: dict[str, Any] = {
            "variant": variant,
            "fields_dir": str(fields_dir),
            "config": str(config_path),
            "output_csv": str(output_csv),
        }
        if not args.prepare_only:
            assert args.runner is not None
            status = run_runner(args.runner, config_path)
            item["status"] = status
            blocks = status.get("block_residuals", {})
            adjusted_blocks = status.get("adjusted_block_residuals", {})
            probe_rows = read_probe_rows(output_csv)
            for node_id, support in sorted(support_by_node.items()):
                row = dict(probe_rows[node_id])
                psi = states[variant]["ElectrostaticPotential"][node_id]
                phin = states[variant]["eQuasiFermiPotential"][node_id]
                phip = states[variant]["hQuasiFermiPotential"][node_id]
                row.update({
                    "psi_V": psi,
                    "phin_V": phin,
                    "phip_V": phip,
                    "node_volume_m2": node_volumes[node_id],
                })
                row.update(srh_decomposition(
                    psi=psi,
                    phin=phin,
                    phip=phip,
                    ni_eff=float(row["ni_eff_m3"]),
                    vt=vt,
                    node_volume=node_volumes[node_id],
                    taun=taun,
                    taup=taup,
                    exact_recombination=finite_float(row.get("electron_recombination")),
                    exact_impact=finite_float(row.get("electron_impact")),
                ))
                row.update(prefixed("electron", impact_closure(
                    flux=float(row["electron_flux"]),
                    recombination=float(row["electron_recombination"]),
                    impact=float(row["electron_impact"]),
                    residual=float(row["electron_residual"]),
                )))
                row.update(prefixed("hole", impact_closure(
                    flux=float(row["hole_flux"]),
                    recombination=float(row["hole_recombination"]),
                    impact=float(row["hole_impact"]),
                    residual=float(row["hole_residual"]),
                )))
                row.update({
                    "variant": variant,
                    "bias_V": args.bias,
                    "support_class": support.get("support_class", ""),
                    "block_psi": blocks.get("psi"),
                    "block_phin": blocks.get("phin"),
                    "block_phip": blocks.get("phip"),
                    "block_combined": blocks.get("combined"),
                    "adjusted_block_psi": adjusted_blocks.get("psi"),
                    "adjusted_block_phin": adjusted_blocks.get("phin"),
                    "adjusted_block_phip": adjusted_blocks.get("phip"),
                    "adjusted_block_combined": adjusted_blocks.get("combined"),
                })
                node_rows.append(row)
        state_summaries.append(item)

    if node_rows:
        write_csv(args.out_dir / "exact_carrier_term_state_nodes.csv", node_rows, NODE_FIELDS)
    payload = clean_json({
        "bias_V": args.bias,
        "variants": VARIANTS,
        "impact_scales": {
            "electron": args.electron_impact_scale,
            "hole": args.hole_impact_scale,
        },
        "preserve_contact_qf_on_shift": args.preserve_contact_qf_on_shift,
        "preserved_shift_node_count": len(preserve_shift_nodes),
        "selected_node_count": len(support_by_node),
        "node_row_count": len(node_rows),
        "states": state_summaries,
        "support_classes": class_summary(node_rows) if node_rows else {},
    })
    (args.out_dir / "exact_carrier_term_state_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
