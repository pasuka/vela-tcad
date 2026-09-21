"""Shared Sentaurus avalanche control deck and VM validation helpers."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path, PurePosixPath


VARIANTS = {
    "implicit_default": {
        "avalanche": "Avalanche(VanOverstraeten)",
        "aval_dens_grad_qf": False,
    },
    "explicit_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": False,
    },
    "explicit_electric_field": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
        "aval_dens_grad_qf": False,
    },
    "grad_qf_aval_dens_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": True,
    },
}


def validate_biases(biases: tuple[int, ...]) -> None:
    if not biases:
        raise ValueError("at least one bias is required")
    if len(set(biases)) != len(biases):
        raise ValueError("biases must be unique")
    if any(bias >= 0 for bias in biases):
        raise ValueError("all avalanche control biases must be negative")
    if tuple(sorted(biases, reverse=True)) != biases:
        raise ValueError("biases must be ordered from low to high magnitude")


def validate_remote_root(value: str) -> str:
    if re.fullmatch(r"/[A-Za-z0-9._/-]+", value) is None:
        raise ValueError("remote root must be a safe absolute POSIX path")
    path = PurePosixPath(value)
    if (
        not path.is_absolute()
        or str(path) != value
        or any(part in {".", ".."} for part in path.parts)
    ):
        raise ValueError("remote root must be a normalized absolute POSIX path")
    return value


def make_solve_block(biases: tuple[int, ...]) -> str:
    validate_biases(biases)
    lines = [
        "Solve {",
        "  Coupled(Iterations=100) { Poisson }",
        "  Coupled(Iterations=100) { Poisson Electron Hole }",
    ]
    for index, bias in enumerate(biases):
        lines.extend(
            [
                "  Quasistationary(",
                (
                    "    InitialStep=1e-4 MinStep=1e-10 MaxStep=0.05"
                    if index == 0
                    else
                    "    InitialStep=1e-3 MinStep=1e-10 MaxStep=0.05"
                ),
                "    Increment=1.2 Decrement=2.0",
                f'    Goal {{ Name="Anode" Voltage={bias} }}',
                "  ) { Coupled { Poisson Electron Hole } }",
            ]
        )
    lines.append("}")
    return "\n".join(lines)


def make_variant_deck(
    template: str,
    variant: str,
    biases: tuple[int, ...],
) -> str:
    if variant not in VARIANTS:
        raise ValueError(f"unknown avalanche drive variant: {variant}")
    validate_biases(biases)
    default_output = "runtime_element_avalanche_probe_default"
    if default_output not in template:
        raise ValueError("default output stem not found in template")
    result = template.replace(
        default_output,
        f"runtime_element_avalanche_probe_{variant}",
    )
    default_avalanche = "Avalanche(VanOverstraeten)"
    if result.count(default_avalanche) != 1:
        raise ValueError("default Avalanche selector must occur exactly once")
    result = result.replace(
        default_avalanche,
        str(VARIANTS[variant]["avalanche"]),
        1,
    )
    if VARIANTS[variant]["aval_dens_grad_qf"]:
        result, math_count = re.subn(
            r"Math\s*\{",
            "Math {\n  AvalDensGradQF",
            result,
            count=1,
        )
        if math_count != 1:
            raise ValueError("Math block was not found exactly once")
    result, solve_count = re.subn(
        r"Solve\s*\{.*\}\s*$",
        make_solve_block(biases),
        result,
        count=1,
        flags=re.S,
    )
    if solve_count != 1:
        raise ValueError("Solve block was not replaced exactly once")
    return result.rstrip() + "\n"


def make_tcl(template: str, biases: tuple[int, ...]) -> str:
    validate_biases(biases)
    expected = "foreach candidate {-1 -10 -20} {"
    replacement = (
        "foreach candidate {"
        + " ".join(str(bias) for bias in biases)
        + "} {"
    )
    if expected not in template:
        raise ValueError("Tcl target list was not found")
    return template.replace(expected, replacement, 1)


def sentaurus_release(ssh_bin: Path, ssh_target: str) -> str:
    command = (
        'resolved=$(readlink -f "$(command -v sdevice)") && '
        'printf "path=%s\\n" "$resolved" && sdevice --version 2>&1'
    )
    completed = subprocess.run(
        [str(ssh_bin), ssh_target, command],
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    match = re.search(r"Version\s+([^\s*]+)", output)
    if match is None:
        raise RuntimeError("Sentaurus release was not found")
    return match.group(1)
