#!/usr/bin/env python3
"""Run exact PN2D Sentaurus avalanche-off/IIC/on process observations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.sentaurus_avalanche_replay import parse_plt
from scripts.pn2d_bv_process_contract import (
    EXACT_BIAS_TOLERANCE_V,
    SCHEMA_ID,
    sha256,
    validate_process_run,
)
from scripts.pn2d_high_bias_process_contract import (
    EXACT_HIGH_BIAS_V,
    SENTAURUS_RELEASE,
)
from scripts.run_pn2d_high_bias_oracle_variant_vm import oracle_deck, oracle_tcl
from scripts.sentaurus_avalanche_controls import (
    validate_biases,
    validate_remote_root,
)


BRANCHES = (
    "avalanche_off",
    "iic_postprocess",
    "avalanche_on",
    "avalanche_on_aval_derivatives",
)
SAFE_FILE_RE = re.compile(r"^[A-Za-z0-9._-]+$")
PROBE_LINE_RE = re.compile(r"^(AVAL_PROBE_[A-Z_]+)\s+(.*)$")
KEY_VALUE_RE = re.compile(r"([A-Za-z0-9_]+)=([^\s]+)")
CURRENT_PLOT_GENERATION_TO_A_PER_UM = 1.6021918e-19 * 1.0e-12


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument(
        "--ssh-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\ssh.exe"),
    )
    parser.add_argument(
        "--scp-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\scp.exe"),
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--template-root",
        type=Path,
        default=Path("build-release/pn2d-minimal6-element-avalanche-replay-20260725"),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--remote-root", required=True)
    parser.add_argument(
        "--biases",
        nargs="+",
        type=float,
        default=list(EXACT_HIGH_BIAS_V),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse complete local or remote branch artifacts after interruption.",
    )
    return parser.parse_args(argv)


def bias_tag(index: int, bias: float) -> str:
    if bias == 0.0:
        return f"snapshot_{index:03d}_zero"
    magnitude = f"{abs(bias):.12g}".replace(".", "p")
    return f"snapshot_{index:03d}_minus{magnitude}"


def validate_process_biases(biases: tuple[float, ...]) -> None:
    if not biases:
        raise ValueError("at least one process bias is required")
    if len(set(biases)) != len(biases):
        raise ValueError("process biases must be unique")
    if any(bias > 0.0 for bias in biases):
        raise ValueError("process biases must be nonpositive")
    if tuple(sorted(biases, reverse=True)) != biases:
        raise ValueError("process biases must be ordered from low to high magnitude")
    if 0.0 in biases and biases[0] != 0.0:
        raise ValueError("equilibrium 0 V must be the first process bias")


def negative_process_biases(biases: tuple[float, ...]) -> tuple[float, ...]:
    negative = tuple(bias for bias in biases if bias < 0.0)
    if not negative:
        raise ValueError("at least one negative process bias is required")
    validate_biases(negative)
    return negative


def exact_solve_block(branch: str, biases: tuple[float, ...]) -> str:
    validate_process_biases(biases)
    lines = [
        "Solve {",
        "  Coupled(Iterations=100) { Poisson }",
        "  Coupled(Iterations=100) { Poisson Electron Hole }",
    ]
    start_index = 0
    if biases[0] == 0.0:
        lines.append(
            f'  Plot(FilePrefix="{bias_tag(0, 0.0)}" NoOverWrite)'
        )
        start_index = 1
    for index, bias in enumerate(biases[start_index:], start=start_index):
        lines.extend(
            [
                "  Quasistationary(",
                (
                    "    InitialStep=1e-4 MinStep=1e-10 MaxStep=0.05"
                    if index == 0
                    else "    InitialStep=1e-3 MinStep=1e-10 MaxStep=0.05"
                ),
                "    Increment=1.2 Decrement=2.0",
                f'    Goal {{ Name="Anode" Voltage={bias:.17g} }}',
                "  ) {",
                "    Coupled { Poisson Electron Hole }",
                "    CurrentPlot(Time=(1))",
                (
                    f'    Plot(FilePrefix="{bias_tag(index, bias)}" '
                    "Time=(1) NoOverWrite)"
                ),
                "  }",
            ]
        )
    lines.append("}")
    return "\n".join(lines)


def _replace_math_open(deck: str, controls: Sequence[str]) -> str:
    if not controls:
        return deck
    replacement = "Math {\n  " + "\n  ".join(controls)
    result, count = re.subn(r"Math\s*\{", replacement, deck, count=1)
    if count != 1:
        raise ValueError("Math block was not found exactly once")
    return result


def _add_currentplot_maxima(deck: str) -> str:
    anchor = """\
  AvalancheGeneration(
    Integrate(Name="AvalancheIntegral" Semiconductor)
  )
"""
    addition = anchor + """\
  ImpactIonization(Maximum(Semiconductor Coordinates))
  ElectricField(Maximum(Semiconductor Coordinates))
"""
    if deck.count(anchor) != 1:
        raise ValueError("CurrentPlot avalanche integral anchor must occur once")
    return deck.replace(anchor, addition, 1)


def make_branch_deck(
    template: str,
    branch: str,
    biases: tuple[float, ...],
) -> str:
    if branch not in BRANCHES:
        raise ValueError(f"unknown branch: {branch}")
    validate_process_biases(biases)
    physics_biases = negative_process_biases(biases)
    builder_variant = (
        "avalanche_disabled" if branch == "avalanche_off" else "explicit_grad_qf"
    )
    deck = oracle_deck(template, builder_variant, physics_biases)
    old_stem = f"runtime_general_tri3_avalanche_probe_{builder_variant}"
    new_stem = f"pn2d_bv_process_{branch}"
    if old_stem not in deck:
        raise ValueError(f"output stem was not found for {branch}")
    deck = deck.replace(old_stem, new_stem)
    controls: list[str] = []
    if branch == "iic_postprocess":
        controls.extend(("ComputeIonizationIntegrals", "AvalPostProcessing"))
    if branch == "avalanche_on_aval_derivatives":
        controls.append("AvalDerivatives")
    if branch in {"avalanche_on", "avalanche_on_aval_derivatives"}:
        controls.extend(
            (
                "CNormPrint",
                "NewtonPlot(Error MinError Residual)",
            )
        )
        file_anchor = f'  Output    = "{new_stem}"\n'
        if deck.count(file_anchor) != 1:
            raise ValueError("File Output anchor must occur once")
        deck = deck.replace(
            file_anchor,
            file_anchor + f'  NewtonPlot = "newton_{branch}_%d_%d_des.tdr"\n',
            1,
        )
    deck = _replace_math_open(deck, controls)
    deck = _add_currentplot_maxima(deck)
    deck, count = re.subn(
        r"Solve\s*\{.*\}\s*$",
        exact_solve_block(branch, biases),
        deck,
        count=1,
        flags=re.S,
    )
    if count != 1:
        raise ValueError("Solve block was not replaced exactly once")
    if "ElectricField)" in re.sub(
        r"ImpactIonization\(Maximum.*?\)",
        "",
        deck,
    ):
        raise ValueError("ElectricField avalanche selector is prohibited")
    return deck.rstrip() + "\n"


def replace_tcl_targets(
    tcl: str,
    old_biases: tuple[float, ...],
    new_biases: tuple[float, ...],
) -> str:
    if old_biases == new_biases:
        return tcl
    old = (
        "foreach candidate {"
        + " ".join(str(bias) for bias in old_biases)
        + "} {"
    )
    new = (
        "foreach candidate {"
        + " ".join(str(bias) for bias in new_biases)
        + "} {"
    )
    if tcl.count(old) != 1:
        raise ValueError("generated Tcl target list was not found exactly once")
    return tcl.replace(old, new, 1)


def make_branch_tcl(template: str, biases: tuple[float, ...]) -> str:
    validate_process_biases(biases)
    physics_biases = negative_process_biases(biases)
    tcl = replace_tcl_targets(
        oracle_tcl(template, physics_biases),
        physics_biases,
        biases,
    )
    return bind_probe_bias_to_anode_contact(tcl)


def bind_probe_bias_to_anode_contact(tcl: str) -> str:
    """Bind process records to the applied Anode QFP, not a domain extremum."""

    old = """\
    set eqfp [$data ReadScalar $::des_data_vertex "eQuasiFermiPotential"]
    set qfp_min 1.0e100
    set qfp_vertex_count [$mesh size_vertex]
    for {set qfp_index 0} {$qfp_index < $qfp_vertex_count} {incr qfp_index} {
        set qfp_value [tcl_cp_get_double $eqfp $qfp_index]
        if {$qfp_value < $qfp_min} {
            set qfp_min $qfp_value
        }
    }

"""
    new = """\
    set eqfp [$data ReadScalar $::des_data_vertex "eQuasiFermiPotential"]
    set qfp_vertex_count [$mesh size_vertex]
    set anode_x 1.0e100
    for {set qfp_index 0} {$qfp_index < $qfp_vertex_count} {incr qfp_index} {
        set qfp_vertex [$mesh vertex $qfp_index]
        set qfp_x [$qfp_vertex coord 0]
        if {$qfp_x < $anode_x} {
            set anode_x $qfp_x
        }
    }
    set anode_qfp_sum 0.0
    set anode_qfp_count 0
    for {set qfp_index 0} {$qfp_index < $qfp_vertex_count} {incr qfp_index} {
        set qfp_vertex [$mesh vertex $qfp_index]
        if {[expr {abs([$qfp_vertex coord 0]-$anode_x)}] < 1.0e-12} {
            set anode_qfp_sum [expr {$anode_qfp_sum
                + [tcl_cp_get_double $eqfp $qfp_index]}]
            incr anode_qfp_count
        }
    }
    if {$anode_qfp_count == 0} {
        return {0.0}
    }
    set anode_qfp [expr {$anode_qfp_sum/double($anode_qfp_count)}]

"""
    if tcl.count(old) != 1:
        raise ValueError("QFP probe-bias observer block was not found exactly once")
    tcl = tcl.replace(old, new, 1)
    old_selector = "abs($qfp_min-$candidate)"
    if tcl.count(old_selector) != 1:
        raise ValueError("QFP probe-bias selector was not found exactly once")
    return tcl.replace(old_selector, "abs($anode_qfp-$candidate)", 1)


def remote_command(remote: str, argv: Sequence[str]) -> list[str]:
    validated = validate_remote_root(remote)
    if not argv or any(
        not isinstance(value, str) or not value or "\n" in value for value in argv
    ):
        raise ValueError("remote argv must contain safe nonempty strings")
    return ["cd", validated, "&&", *argv]


def remote_shell_text(remote: str, argv: Sequence[str]) -> str:
    command = remote_command(remote, argv)
    return (
        f"cd {shlex.quote(command[1])} && "
        + " ".join(shlex.quote(value) for value in command[3:])
    )


def run_checked(argv: Sequence[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(argv),
        check=True,
        capture_output=True,
        text=True,
    )


def probe_sentaurus_release(ssh_bin: Path, ssh_target: str) -> str:
    """Read the release header; O-2018 prints it before rejecting --version."""

    completed = subprocess.run(
        [str(ssh_bin), ssh_target, "sdevice --version 2>&1"],
        check=False,
        capture_output=True,
        text=True,
    )
    match = re.search(
        r"Version\s+([^\s*]+)",
        completed.stdout + completed.stderr,
    )
    if match is None:
        detail = (completed.stdout + completed.stderr).strip()
        raise RuntimeError(
            "Sentaurus release was not found: "
            + (detail[-1000:] if detail else "<empty ssh output>")
        )
    return match.group(1)


def write_ascii(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="ascii", newline="\n")


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        write_ascii(path, "")
        return
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def exact_currentplot_rows(
    plt_path: Path,
    biases: tuple[float, ...],
) -> list[dict[str, float]]:
    names, rows = parse_plt(plt_path)
    voltage_name = next(
        (name for name in names if name == "Anode OuterVoltage"),
        None,
    )
    if voltage_name is None:
        raise ValueError(f"{plt_path}: missing Anode OuterVoltage")
    selected: list[dict[str, float]] = []
    for requested in biases:
        distances = [abs(row[voltage_name] - requested) for row in rows]
        best = min(distances)
        if best > EXACT_BIAS_TOLERANCE_V:
            raise ValueError(
                f"{plt_path}: nearest_bias_substitution requested={requested}, "
                f"nearest_error={best}"
            )
        matches = [
            row
            for row, distance in zip(rows, distances, strict=True)
            if distance <= EXACT_BIAS_TOLERANCE_V
        ]
        selected.append(
            {
                "requested_bias_V": requested,
                "actual_bias_V": matches[-1][voltage_name],
                **matches[-1],
            }
        )
    return selected


def parse_probe_lines(
    run_text: str,
    branch: str,
    actual_by_requested: Mapping[float, float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    fields: list[dict[str, Any]] = []
    aggregates: list[dict[str, Any]] = []
    raw_lines: list[str] = []
    cell_connectivity: dict[tuple[float, int], set[int]] = {}
    node_coordinates: dict[tuple[float, int], list[float]] = {}

    def emit_field(
        *,
        requested: float,
        source_index: int,
        support_kind: str,
        support_key: str,
        centering: str,
        provenance: str,
        carrier: str,
        quantity: str,
        components: list[str],
        unit: str,
        values: list[float],
        coordinates_um: list[float] | None = None,
        connectivity: list[int] | None = None,
    ) -> None:
        record: dict[str, Any] = {
            "branch": branch,
            "requested_bias_V": requested,
            "actual_bias_V": actual_by_requested[requested],
            "support_kind": support_kind,
            "support_key": support_key,
            "centering": centering,
            "provenance": provenance,
            "carrier": carrier,
            "quantity": quantity,
            "components": components,
            "unit": unit,
            "values": values,
            "source": {
                "file": f"{branch}/raw/run.out",
                "dataset": "runtime_tcl",
                "index": source_index,
            },
        }
        if coordinates_um is not None:
            record["coordinates_um"] = coordinates_um
        if connectivity is not None:
            record["connectivity"] = connectivity
        fields.append(record)

    for source_index, line in enumerate(run_text.splitlines()):
        match = PROBE_LINE_RE.match(line)
        if match is None:
            continue
        record_type, body = match.groups()
        values = dict(KEY_VALUE_RE.findall(body))
        if "bias_V" not in values:
            continue
        requested = float(values["bias_V"])
        if requested not in actual_by_requested:
            continue
        raw_lines.append(line)
        if record_type == "AVAL_PROBE_VERTEX":
            key = f"node:{int(values['vertex'])}"
            coordinates = [float(values["x_um"]), float(values["y_um"])]
            node_coordinates[
                (requested, int(values["vertex"]))
            ] = coordinates
            specifications = (
                ("none", "potential", "V", "psi_V"),
                ("electron", "density", "cm^-3", "n_cm3"),
                ("hole", "density", "cm^-3", "p_cm3"),
                ("electron", "quasi_fermi", "V", "eQFP_V"),
                ("hole", "quasi_fermi", "V", "hQFP_V"),
                ("electron", "avalanche_alpha", "cm^-1", "alpha_n_cm_inv"),
                ("hole", "avalanche_alpha", "cm^-1", "alpha_p_cm_inv"),
                (
                    "electron",
                    "avalanche_generation",
                    "cm^-3 s^-1",
                    "generation_n_cm3_s",
                ),
                (
                    "hole",
                    "avalanche_generation",
                    "cm^-3 s^-1",
                    "generation_p_cm3_s",
                ),
                (
                    "total",
                    "avalanche_generation",
                    "cm^-3 s^-1",
                    "generation_total_cm3_s",
                ),
            )
            for carrier, quantity, unit, name in specifications:
                emit_field(
                    requested=requested,
                    source_index=source_index,
                    support_kind="physical_node",
                    support_key=key,
                    centering="vertex",
                    provenance="native",
                    carrier=carrier,
                    quantity=quantity,
                    components=["scalar"],
                    unit=unit,
                    values=[float(values[name])],
                    coordinates_um=coordinates,
                )
        elif record_type == "AVAL_PROBE_ELEMENT":
            key = f"cell:{int(values['element'])}"
            specifications = (
                (
                    "electron",
                    "mobility",
                    "cm^2/(V s)",
                    ["scalar"],
                    ["mu_n_cm2_Vs"],
                ),
                (
                    "hole",
                    "mobility",
                    "cm^2/(V s)",
                    ["scalar"],
                    ["mu_p_cm2_Vs"],
                ),
                (
                    "none",
                    "electric_field",
                    "V/cm",
                    ["x", "y"],
                    ["efield_x_V_cm", "efield_y_V_cm"],
                ),
                (
                    "electron",
                    "quasi_fermi_gradient",
                    "V/cm",
                    ["x", "y"],
                    ["grad_qf_n_x_V_cm", "grad_qf_n_y_V_cm"],
                ),
                (
                    "hole",
                    "quasi_fermi_gradient",
                    "V/cm",
                    ["x", "y"],
                    ["grad_qf_p_x_V_cm", "grad_qf_p_y_V_cm"],
                ),
                (
                    "electron",
                    "current_density",
                    "A/cm^2",
                    ["x", "y"],
                    ["current_n_x_A_cm2", "current_n_y_A_cm2"],
                ),
                (
                    "hole",
                    "current_density",
                    "A/cm^2",
                    ["x", "y"],
                    ["current_p_x_A_cm2", "current_p_y_A_cm2"],
                ),
            )
            for carrier, quantity, unit, components, names in specifications:
                emit_field(
                    requested=requested,
                    source_index=source_index,
                    support_kind="cell",
                    support_key=key,
                    centering="cell",
                    provenance="native",
                    carrier=carrier,
                    quantity=quantity,
                    components=components,
                    unit=unit,
                    values=[float(values[name]) for name in names],
                )
        elif record_type == "AVAL_PROBE_EDGE":
            element = int(values["element"])
            local_edge = int(values["local_edge"])
            start = int(values["start"])
            end = int(values["end"])
            cell_connectivity.setdefault((requested, element), set()).update(
                (start, end)
            )
            key = f"cell:{element}/local_edge:{local_edge}"
            coordinates = [
                float(values["start_x_um"]),
                float(values["start_y_um"]),
                float(values["end_x_um"]),
                float(values["end_y_um"]),
            ]
            for provenance, carrier, name in (
                ("operator_replay", "electron", "sg_jn_A_cm2"),
                ("operator_replay", "hole", "sg_jp_A_cm2"),
                ("reconstructed", "electron", "native_tangent_n_A_cm2"),
                ("reconstructed", "hole", "native_tangent_p_A_cm2"),
            ):
                emit_field(
                    requested=requested,
                    source_index=source_index,
                    support_kind="element_local_edge",
                    support_key=key,
                    centering="element_edge",
                    provenance=provenance,
                    carrier=carrier,
                    quantity="current_density",
                    components=["tangent"],
                    unit="A/cm^2",
                    values=[float(values[name])],
                    coordinates_um=coordinates,
                    connectivity=[start, end],
                )
        elif record_type == "AVAL_PROBE_PROCESS":
            key = f"node:{int(values['vertex'])}"
            specifications = (
                ("electron", "velocity", "cm/s", ["velocity_n_cm_s"]),
                ("hole", "velocity", "cm/s", ["velocity_p_cm_s"]),
                (
                    "total",
                    "current_density",
                    "A/cm^2",
                    ["total_current_x_A_cm2", "total_current_y_A_cm2"],
                ),
                ("electron", "ionization_integral", "1", ["ion_n"]),
                ("hole", "ionization_integral", "1", ["ion_p"]),
                ("total", "ionization_integral", "1", ["ion_mean"]),
                ("none", "doping", "cm^-3", ["doping_cm3"]),
                ("none", "charge_density", "cm^-3", ["space_charge_cm3"]),
                ("total", "srh_recombination", "cm^-3 s^-1", ["srh_cm3_s"]),
            )
            for carrier, quantity, unit, names in specifications:
                emit_field(
                    requested=requested,
                    source_index=source_index,
                    support_kind="physical_node",
                    support_key=key,
                    centering="vertex",
                    provenance="native",
                    carrier=carrier,
                    quantity=quantity,
                    components=(
                        ["x", "y"] if len(names) == 2 else ["scalar"]
                    ),
                    unit=unit,
                    values=[float(values[name]) for name in names],
                )
        elif record_type == "AVAL_PROBE_MEASURE":
            key = (
                f"cell:{int(values['element'])}/"
                f"local_vertex:{int(values['local_vertex'])}"
            )
            for carrier, name in (
                ("electron", "qg_n_A_um"),
                ("hole", "qg_p_A_um"),
                ("total", "qg_total_A_um"),
            ):
                emit_field(
                    requested=requested,
                    source_index=source_index,
                    support_kind="element_local_vertex",
                    support_key=key,
                    centering="element_vertex",
                    provenance="operator_replay",
                    carrier=carrier,
                    quantity="integrated_source",
                    components=["scalar"],
                    unit="A/um",
                    values=[float(values[name])],
                    connectivity=[int(values["vertex"])],
                )
        elif record_type == "AVAL_PROBE_INTEGRAL":
            for carrier, name in (
                ("electron", "qg_n_A_um"),
                ("hole", "qg_p_A_um"),
                ("total", "qg_total_A_um"),
            ):
                aggregates.append(
                    {
                        "branch": branch,
                        "requested_bias_V": requested,
                        "actual_bias_V": actual_by_requested[requested],
                        "carrier": carrier,
                        "quantity": "integrated_source",
                        "unit": "A/um",
                        "value": float(values[name]),
                        "provenance": "operator_replay",
                        "source": {
                            "file": f"{branch}/raw/run.out",
                            "dataset": record_type,
                            "index": source_index,
                        },
                    }
                )
    expected = set(actual_by_requested)
    begin_biases = {
        float(KEY_VALUE_RE.findall(line)[0][1])
        for line in raw_lines
        if line.startswith("AVAL_PROBE_BEGIN ")
    }
    if begin_biases != expected:
        raise ValueError(
            f"{branch}: process records mismatch: expected {expected}, got {begin_biases}"
        )
    for field in fields:
        if field["support_kind"] == "cell":
            key = (
                float(field["requested_bias_V"]),
                int(str(field["support_key"]).split(":", 1)[1]),
            )
            nodes = sorted(cell_connectivity.get(key, ()))
            if not nodes:
                raise ValueError(
                    f"{branch}: missing connectivity for {field['support_key']}"
                )
            field["connectivity"] = nodes
        elif (
            field["support_kind"] == "physical_node"
            and "coordinates_um" not in field
        ):
            key = (
                float(field["requested_bias_V"]),
                int(str(field["support_key"]).split(":", 1)[1]),
            )
            coordinates = node_coordinates.get(key)
            if coordinates is None:
                raise ValueError(
                    f"{branch}: missing coordinates for {field['support_key']}"
                )
            field["coordinates_um"] = coordinates
    return fields, aggregates, raw_lines


def currentplot_aggregates(
    branch: str,
    rows: Sequence[Mapping[str, float]],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    terminal_names = {
        "Anode eCurrent": "electron",
        "Anode hCurrent": "hole",
        "Anode TotalCurrent": "total",
    }
    for source_index, row in enumerate(rows):
        requested = float(row["requested_bias_V"])
        actual = float(row["actual_bias_V"])
        for name, carrier in terminal_names.items():
            if name in row:
                result.append(
                    {
                        "branch": branch,
                        "requested_bias_V": requested,
                        "actual_bias_V": actual,
                        "carrier": carrier,
                        "quantity": "terminal_current",
                        "unit": "A/um",
                        "value": float(row[name]),
                        "provenance": "native",
                        "source": {
                            "file": f"{branch}/raw/pn2d_bv_process_{branch}.plt",
                            "dataset": name,
                            "index": source_index,
                        },
                    }
                )
        generation_names = (
            ("IntegreAvalancheIntegral eAvalancheGeneration", "electron"),
            ("IntegrhAvalancheIntegral hAvalancheGeneration", "hole"),
            ("IntegrAvalancheIntegral AvalancheGeneration", "total"),
        )
        for name, carrier in generation_names:
            if name in row:
                result.append(
                    {
                        "branch": branch,
                        "requested_bias_V": requested,
                        "actual_bias_V": actual,
                        "carrier": carrier,
                        "quantity": "integrated_source",
                        "unit": "A/um",
                        "value": (
                            float(row[name])
                            * CURRENT_PLOT_GENERATION_TO_A_PER_UM
                        ),
                        "provenance": "native",
                        "source": {
                            "file": f"{branch}/raw/pn2d_bv_process_{branch}.plt",
                            "dataset": name,
                            "index": source_index,
                        },
                    }
                )
    return result


def list_remote_files(
    ssh_bin: Path,
    ssh_target: str,
    remote: str,
) -> list[str]:
    completed = run_checked(
        [
            str(ssh_bin),
            ssh_target,
            remote_shell_text(
                remote,
                ["find", ".", "-maxdepth", "1", "-type", "f", "-printf", "%f\\n"],
            ),
        ]
    )
    names = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    if any(SAFE_FILE_RE.fullmatch(name) is None for name in names):
        raise ValueError("remote output contains an unsafe file name")
    return sorted(set(names))


def copy_remote_files(
    ssh_bin: Path,
    scp_bin: Path,
    ssh_target: str,
    remote: str,
    names: Sequence[str],
    destination: Path,
) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for name in names:
        if SAFE_FILE_RE.fullmatch(name) is None:
            raise ValueError(f"unsafe remote artifact name: {name!r}")
    if not names:
        return
    archive_name = "pn2d_process_artifacts.tar.gz"
    run_checked(
        [
            str(ssh_bin),
            ssh_target,
            remote_shell_text(
                remote,
                ["tar", "-czf", archive_name, "--", *names],
            ),
        ]
    )
    archive_path = destination / archive_name
    run_checked(
        [
            str(scp_bin),
            f"{ssh_target}:{remote}/{archive_name}",
            str(archive_path),
        ]
    )
    with tarfile.open(archive_path, "r:gz") as archive:
        members = archive.getmembers()
        archived_names = [member.name for member in members]
        if any(
            not member.isfile()
            or SAFE_FILE_RE.fullmatch(member.name) is None
            or member.name not in names
            for member in members
        ):
            raise ValueError("remote artifact archive contains an unsafe member")
        if sorted(archived_names) != sorted(names):
            raise ValueError("remote artifact archive is incomplete")
        archive.extractall(destination, members=members, filter="data")


def completed_run_text(path: Path) -> bool:
    return path.is_file() and "Sentaurus Device simulation finished" in path.read_text(
        encoding="ascii",
        errors="ignore",
    )


def remote_run_completed(
    ssh_bin: Path,
    ssh_target: str,
    remote: str,
) -> bool:
    completed = subprocess.run(
        [
            str(ssh_bin),
            ssh_target,
            remote_shell_text(
                remote,
                [
                    "grep",
                    "-q",
                    "Sentaurus Device simulation finished",
                    "run.out",
                ],
            ),
        ],
        check=False,
        capture_output=True,
    )
    return completed.returncode == 0


def artifact(root: Path, path: Path) -> dict[str, str]:
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256(path),
    }


def normalized_case(
    root: Path,
    branch: str,
    biases: tuple[float, ...],
    bundle: Path,
    raw: Path,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, str],
]:
    stem = f"pn2d_bv_process_{branch}"
    plt_path = raw / f"{stem}.plt"
    run_path = raw / "run.out"
    rows = exact_currentplot_rows(plt_path, biases)
    actual_by_requested = {
        float(row["requested_bias_V"]): float(row["actual_bias_V"]) for row in rows
    }
    fields, replay_aggregates, process_lines = parse_probe_lines(
        run_path.read_text(encoding="ascii", errors="strict"),
        branch,
        actual_by_requested,
    )
    aggregates = currentplot_aggregates(branch, rows) + replay_aggregates
    normalized = root / branch / "normalized"
    write_csv(normalized / "currentplot_exact.csv", rows)
    write_ascii(
        normalized / "process_records.txt",
        "\n".join(sorted(process_lines)) + "\n",
    )
    write_ascii(
        normalized / "fields.jsonl",
        "".join(
            json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n"
            for row in sorted(
                fields,
                key=lambda row: (
                    row["actual_bias_V"],
                    row["support_kind"],
                    row["support_key"],
                    row["quantity"],
                    row["carrier"],
                ),
            )
        ),
    )
    write_csv(normalized / "aggregate.csv", aggregates)

    bias_records = []
    for index, requested in enumerate(biases):
        prefix = bias_tag(index, requested)
        snapshots = sorted(raw.glob(f"{prefix}*_des.tdr"))
        if len(snapshots) != 1:
            raise ValueError(
                f"{branch} {requested:g} V: expected one snapshot, got {snapshots}"
            )
        bias_records.append(
            {
                "requested_bias_V": requested,
                "actual_bias_V": actual_by_requested[requested],
                "snapshot_tdr": artifact(root, snapshots[0]),
                "currentplot": artifact(root, plt_path),
                "process_record": artifact(root, run_path),
            }
        )
    hashes = {
        path.relative_to(root).as_posix(): sha256(path)
        for path in sorted(normalized.iterdir())
        if path.is_file()
    }
    branch_record = {
        "branch": branch,
        "requested_biases_V": list(biases),
        "bias_records": bias_records,
    }
    return branch_record, fields, aggregates, hashes


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    biases = tuple(args.biases)
    validate_process_biases(biases)
    remote_root = validate_remote_root(args.remote_root)
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_root = args.source_root.resolve()
    template_root = args.template_root.resolve()
    sources = {
        "pn2d_msh.tdr": source_root / "pn2d_msh.tdr",
        "models.par": source_root / "models.par",
    }
    for name, path in sources.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing process input {name}: {path}")
    deck_template = (
        template_root / "runtime_element_avalanche_probe_default.cmd"
    ).read_text(encoding="ascii")
    tcl_template = (
        template_root / "runtime_element_avalanche_probe.tcl"
    ).read_text(encoding="ascii")

    for branch in BRANCHES:
        case = output / branch
        bundle = case / "bundle"
        bundle.mkdir(parents=True, exist_ok=True)
        for name, path in sources.items():
            shutil.copy2(path, bundle / name)
        write_ascii(
            bundle / "runtime_element_avalanche_probe.tcl",
            make_branch_tcl(tcl_template, biases),
        )
        write_ascii(
            bundle / f"pn2d_bv_process_{branch}.cmd",
            make_branch_deck(deck_template, branch, biases),
        )

    if args.dry_run:
        write_ascii(
            output / "dry_run_manifest.json",
            json.dumps(
                {
                    "schema": "vela.pn2d_bv_process_matrix_dry_run.v1",
                    "branches": list(BRANCHES),
                    "biases_V": list(biases),
                    "remote_root": remote_root,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        return 0

    release = probe_sentaurus_release(args.ssh_bin, args.ssh_target)
    if release != SENTAURUS_RELEASE:
        raise RuntimeError(f"expected Sentaurus {SENTAURUS_RELEASE}, got {release}")

    for branch in BRANCHES:
        case = output / branch
        bundle = case / "bundle"
        raw = case / "raw"
        remote = f"{remote_root}/{branch}"
        if args.resume and completed_run_text(raw / "run.out"):
            continue
        run_checked([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote}"])
        for path in sorted(bundle.iterdir()):
            run_checked(
                [
                    str(args.scp_bin),
                    str(path),
                    f"{args.ssh_target}:{remote}/{path.name}",
                ]
            )
        deck_name = f"pn2d_bv_process_{branch}.cmd"
        if args.resume and remote_run_completed(
            args.ssh_bin,
            args.ssh_target,
            remote,
        ):
            returncode = 0
        else:
            completed = subprocess.run(
                [
                    str(args.ssh_bin),
                    args.ssh_target,
                    remote_shell_text(
                        remote,
                        ["sdevice", deck_name],
                    )
                    + " > run.out 2>&1",
                ],
                check=False,
            )
            returncode = completed.returncode
        names = list_remote_files(args.ssh_bin, args.ssh_target, remote)
        retained = [
            name
            for name in names
            if name == "run.out"
            or name.endswith((".plt", ".tdr", ".log", ".cmd", ".par", ".tcl"))
        ]
        copy_remote_files(
            args.ssh_bin,
            args.scp_bin,
            args.ssh_target,
            remote,
            retained,
            raw,
        )
        if returncode:
            raise RuntimeError(
                f"{branch}: sdevice failed with {returncode}; "
                f"artifacts retained in {raw}"
            )

    branch_records: list[dict[str, Any]] = []
    field_records: list[dict[str, Any]] = []
    aggregate_records: list[dict[str, Any]] = []
    normalized_hashes: dict[str, str] = {}
    input_hashes: dict[str, str] = {}
    for branch in BRANCHES:
        case = output / branch
        bundle = case / "bundle"
        branch_record, fields, aggregates, hashes = normalized_case(
            output,
            branch,
            biases,
            bundle,
            case / "raw",
        )
        branch_records.append(branch_record)
        field_records.extend(fields)
        aggregate_records.extend(aggregates)
        normalized_hashes.update(hashes)
        input_hashes.update(
            {
                path.relative_to(output).as_posix(): sha256(path)
                for path in sorted(bundle.iterdir())
                if path.is_file()
            }
        )

    manifest: dict[str, Any] = {
        "schema": SCHEMA_ID,
        "status": "passed",
        "outcome": "sentaurus_process_matrix_available",
        "run_id": output.name,
        "simulator": "sentaurus",
        "release": release,
        "missing_value_policy": "reject",
        "input_hashes": input_hashes,
        "normalized_output_hashes": normalized_hashes,
        "branch_records": branch_records,
        "field_records": field_records,
        "aggregate_records": aggregate_records,
        "newton_attempt_records": [],
    }
    validate_process_run(manifest, base_dir=output)
    write_ascii(
        output / "manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    print(
        json.dumps(
            {
                "outcome": manifest["outcome"],
                "root": str(output),
                "field_records": len(field_records),
                "aggregate_records": len(aggregate_records),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
