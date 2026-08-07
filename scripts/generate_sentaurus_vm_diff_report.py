#!/usr/bin/env python3
"""Generate a machine-readable diff report for Sentaurus VM vs Vela outputs."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from compare_sentaurus_vm_oracle import compare_series

def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


def fail_report(args: argparse.Namespace, item: str, reason: str, error: str) -> int:
    report = {
        "schema": "vela.sentaurus_vm.diff_report.v1",
        "sentaurus_dir": str(args.sentaurus_dir),
        "vela_dir": str(args.vela_dir),
        "items": [],
        "overall_pass": False,
        "failures": [{"item": item, "reason": reason, "error": error}],
        "missing_outputs": [],
        "expectations": {
            "sentaurus_dir": [
                "nodes.csv",
                "elements.csv",
                "contacts.csv",
                "doping.csv",
                "pn2d_sentaurus_inventory.json",
                "reference_curves/pn2d_idvd_reference.csv",
            ],
            "vela_dir": [
                "nodes.csv",
                "elements.csv",
                "contacts.csv",
                "doping.csv",
                "inventory.json",
                "pn2d_idvd.csv",
            ],
        },
    }
    write_report(args.output_json, report)
    print(json.dumps(report, indent=2))
    return 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-dir", required=True, type=Path)
    parser.add_argument("--vela-dir", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--tolerance", required=True, type=float)
    parser.add_argument("--x-tolerance", type=float, default=0.0)
    args = parser.parse_args()
    try:
        if not args.output_json.parent.exists():
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text("", encoding="utf-8")
    except Exception as exc:
        print(f"failed to prepare output report: {exc}", file=sys.stderr)
        return 1

    for name, value in (("tolerance", args.tolerance), ("x_tolerance", args.x_tolerance)):
        if not math.isfinite(value) or value < 0.0:
            return fail_report(
                args,
                "arguments",
                "input_error",
                f"{name} must be finite and >= 0",
            )

    report = {
        "schema": "vela.sentaurus_vm.diff_report.v1",
        "sentaurus_dir": str(args.sentaurus_dir),
        "vela_dir": str(args.vela_dir),
        "items": [],
    }

    report["overall_pass"] = True
    report["failures"] = []
    report["missing_outputs"] = []

    try:
        for name in ["nodes.csv", "elements.csv", "contacts.csv", "doping.csv"]:
            sentaurus_path = args.sentaurus_dir / name
            vela_path = args.vela_dir / name
            if not sentaurus_path.exists():
                report["items"].append({
                    "type": "csv",
                    "name": name,
                    "status": "missing_oracle",
                    "sentaurus_exists": False,
                    "vela_exists": vela_path.exists(),
                })
                report["overall_pass"] = False
                report["failures"].append({"item": name, "reason": "missing_oracle"})
                report["missing_outputs"].append(name)
            elif not vela_path.exists():
                report["items"].append({
                    "type": "csv",
                    "name": name,
                    "status": "missing_candidate",
                    "sentaurus_exists": True,
                    "vela_exists": False,
                })
                report["overall_pass"] = False
                report["failures"].append({"item": name, "reason": "missing_candidate"})
                report["missing_outputs"].append(name)
            else:
                oracle_lines = sentaurus_path.read_text(encoding="utf-8").splitlines()
                vela_lines = vela_path.read_text(encoding="utf-8").splitlines()
                identical = vela_lines == oracle_lines
                report["items"].append({
                    "type": "csv",
                    "name": name,
                    "vela_path": str(vela_path),
                    "oracle_path": str(sentaurus_path),
                    "vela_lines": len(vela_lines),
                    "oracle_lines": len(oracle_lines),
                    "identical": identical,
                })
                if not identical:
                    report["overall_pass"] = False
                    report["failures"].append({"item": name, "reason": "csv_diff"})

        sentaurus_inventory = args.sentaurus_dir / "pn2d_sentaurus_inventory.json"
        if not sentaurus_inventory.exists():
            sentaurus_inventory = args.sentaurus_dir / "nmos2d_sentaurus_inventory.json"
        vela_inventory = args.vela_dir / "inventory.json"
        if not sentaurus_inventory.exists():
            report["items"].append({"type": "json", "name": "inventory", "status": "missing_oracle"})
            report["overall_pass"] = False
            report["failures"].append({"item": "inventory", "reason": "missing_oracle"})
            report["missing_outputs"].append("inventory")
        elif not vela_inventory.exists():
            report["items"].append({"type": "json", "name": "inventory", "status": "missing_candidate"})
            report["overall_pass"] = False
            report["failures"].append({"item": "inventory", "reason": "missing_candidate"})
            report["missing_outputs"].append("inventory")
        else:
            sentaurus_inventory_data = load_json(sentaurus_inventory)
            vela_inventory_data = load_json(vela_inventory)
            inventory_equal = sentaurus_inventory_data == vela_inventory_data
            report["items"].append({
                "type": "json",
                "name": "inventory",
                "sentaurus": sentaurus_inventory_data,
                "vela": vela_inventory_data,
                "identical": inventory_equal,
            })
            if not inventory_equal:
                report["overall_pass"] = False
                report["failures"].append({"item": "inventory", "reason": "json_diff"})

        vela_curve = args.vela_dir / "pn2d_idvd.csv"
        oracle_curve = args.sentaurus_dir / "reference_curves" / "pn2d_idvd_reference.csv"
        if vela_curve.exists() and oracle_curve.exists():
            curve_compare = compare_series(
                vela_curve,
                oracle_curve,
                "voltage_V",
                "current_A_per_um",
                args.tolerance,
                args.x_tolerance,
            )
            report["items"].append({
                "type": "curve",
                "name": "pn2d_idvd",
                "compare": curve_compare,
            })
            if not curve_compare["pass"]:
                report["overall_pass"] = False
                report["failures"].append({"item": "pn2d_idvd", "reason": "curve_diff"})
        elif oracle_curve.exists() or vela_curve.exists():
            report["items"].append({
                "type": "curve",
                "name": "pn2d_idvd",
                "status": "missing_candidate" if oracle_curve.exists() and not vela_curve.exists() else "missing_oracle",
            })
            report["overall_pass"] = False
            report["failures"].append({"item": "pn2d_idvd", "reason": "missing_curve"})
            report["missing_outputs"].append("pn2d_idvd")
        else:
            report["items"].append({
                "type": "curve",
                "name": "pn2d_idvd",
                "status": "missing_both",
            })
            report["overall_pass"] = False
            report["failures"].append({"item": "pn2d_idvd", "reason": "missing_both"})
            report["missing_outputs"].append("pn2d_idvd")
    except Exception as exc:
        return fail_report(args, "report", "input_error", str(exc))

    report["expectations"] = {
        "sentaurus_dir": [
            "nodes.csv",
            "elements.csv",
            "contacts.csv",
            "doping.csv",
            "pn2d_sentaurus_inventory.json",
            "reference_curves/pn2d_idvd_reference.csv",
        ],
        "vela_dir": [
            "nodes.csv",
            "elements.csv",
            "contacts.csv",
            "doping.csv",
            "inventory.json",
            "pn2d_idvd.csv",
        ],
    }

    write_report(args.output_json, report)
    print(json.dumps(report, indent=2))
    return 0 if report["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
