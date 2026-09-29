#!/usr/bin/env python3
"""Generate ``share/vela/materials/models.par`` from Vela's built-in materials.

``vela device -P`` exports this file, mirroring ``sdevice -P``.  It expresses
the C++ ``MaterialDatabase`` built-ins (Si, SiO2) in Sentaurus ``.par`` syntax
using only the sections the WP3 emitter implements: ``Epsilon``, ``Bandgap``
(``alpha = 0`` so the gap is temperature independent), ``eDOSMass``/``hDOSMass``
Formula 1 with the density of states carried by the additive mass term ``mm``,
and ``ConstantMobility`` with Vela's ``(T/300 K)^-2.2`` exponent.

The built-in silicon set is not thermodynamically closed: ``sqrt(Nc Nv)
exp(-Eg/2kT)`` with the built-in ``Nc``, ``Nv`` and ``Eg`` gives
``6.68e9 cm^-3`` while the built-in ``ni`` is ``1.0e10 cm^-3``.  The emitter
always derives ``ni``, so the round trip reproduces every built-in field
except ``ni`` (see the header comment written into the file).

Usage::

    python scripts/generate_default_models_par.py            # rewrite the file
    python scripts/generate_default_models_par.py --check    # exit 1 if stale
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from vela_sdevice.materials import CONSTANTS, dos_prefactor_cm3, intrinsic_density_cm3  # noqa: E402

MATERIAL_DATABASE = REPO / "src" / "material" / "MaterialDatabase.cpp"
OUTPUT = REPO / "share" / "vela" / "materials" / "models.par"
REFERENCE_T = 300.0
#: (C++ variable, Sentaurus material scope, Vela name)
MATERIALS = (("si", "Silicon", "Si"), ("sio2", "Oxide", "SiO2"))
_ASSIGN = re.compile(r"^\s*(?P<var>si|sio2)\.(?P<field>\w+)\s*=\s*(?P<value>[-+0-9.eE]+)\s*;", re.M)


def builtin_materials(source: Path = MATERIAL_DATABASE) -> dict:
    """Read the ``MaterialDatabase()`` constructor literals (SI units)."""
    text = source.read_text(encoding="utf-8")
    start = text.index("MaterialDatabase::MaterialDatabase()")
    end = text.index("}", text.index('db_["SiO2"]', start))
    values: dict = {}
    for match in _ASSIGN.finditer(text[start:end]):
        values.setdefault(match["var"], {})[match["field"]] = float(match["value"])
    for var in ("si", "sio2"):
        if var not in values:
            raise ValueError(f"MaterialDatabase constructor does not define {var}")
    return values


def effective_mass(dos_m3: float) -> float:
    """Additive Formula 1 mass ``mm`` that reproduces ``dos_m3`` at 300 K."""
    prefactor_m3 = dos_prefactor_cm3(REFERENCE_T) * 1.0e6
    return (dos_m3 / prefactor_m3) ** (2.0 / 3.0)


def _num(value: float) -> str:
    return repr(float(value))


def render(values: dict) -> str:
    si = values["si"]
    sio2 = values["sio2"]
    kT = CONSTANTS["kb_J_per_K"] / CONSTANTS["q_C"] * REFERENCE_T
    nc_cm3 = si["Nc_m3"] / 1.0e6
    nv_cm3 = si["Nv_m3"] / 1.0e6
    derived_ni = intrinsic_density_cm3(nc_cm3, nv_cm3, si["bandgap_eV"], REFERENCE_T)
    me = effective_mass(si["Nc_m3"])
    mh = effective_mass(si["Nv_m3"])
    lines = [
        "* Vela default material parameters (generated; do not edit by hand).",
        "* Source: src/material/MaterialDatabase.cpp built-in Si/SiO2, rendered by",
        "* scripts/generate_default_models_par.py in Sentaurus .par syntax using only the",
        "* sections the `vela device` materials emitter implements.",
        "* Densities of states are carried by the additive Formula 1 mass term mm so that",
        "* Nc(300 K) = 2(2 pi m0 mm kB T/h^2)^1.5 reproduces the built-in values.",
        f"* NOTE: the built-in ni = {si['ni'] / 1.0e6:.6g} cm^-3 is not derivable from the built-in",
        f"* Nc/Nv/Eg (sqrt(Nc Nv) exp(-Eg/2kT) = {derived_ni:.6g} cm^-3 at 300 K, kT = {kT:.6g} eV);",
        "* the emitter always derives ni, so this file reproduces every built-in field except ni.",
        "",
        f'Material = "{MATERIALS[0][1]}" {{',
        "",
        "Epsilon",
        "{ *  Ratio of the permittivities of material and vacuum",
        f"\tepsilon\t= {_num(si['eps_r'])}\t# [1]",
        "}",
        "",
        "Bandgap",
        "{ * Eg = Eg0 + dEg0 + alpha Tpar^2 / (beta + Tpar) - alpha T^2 / (beta + T)",
        "  * alpha = 0 keeps the built-in gap temperature independent.",
        f"\tChi0\t= {_num(si['electron_affinity_eV'])}\t# [eV]",
        f"\tEg0\t= {_num(si['bandgap_eV'])}\t# [eV]",
        "\talpha\t= 0.0\t# [eV K^-1]",
        "\tbeta\t= 636.0\t# [K]",
        f"\tTpar\t= {_num(REFERENCE_T)}\t# [K]",
        "}",
        "",
        "eDOSMass",
        "{ * Formula1: me/m0 = [ (6 * mt)^2 * ml ]^(1/3) + mm, mt = a[Eg(0)/Eg(T)]",
        f"  * a = ml = 0 so me/m0 = mm reproduces the built-in Nc = {nc_cm3:.6g} cm^-3 at 300 K.",
        "\tFormula\t= 1\t# [1]",
        "\ta\t= 0.0\t# [1]",
        "\tml\t= 0.0\t# [1]",
        f"\tmm\t= {_num(me)}\t# [1]",
        "}",
        "",
        "hDOSMass",
        "{ * Formula1: mh/m0 = [(a+bT+cT^2+dT^3+eT^4)/(1+fT+gT^2+hT^3+iT^4)]^(2/3) + mm",
        f"  * a..i = 0 so mh/m0 = mm reproduces the built-in Nv = {nv_cm3:.6g} cm^-3 at 300 K.",
        "\tFormula\t= 1\t# [1]",
        "\ta\t= 0.0\t# [1]",
        "\tb\t= 0.0\t# [K^-1]",
        "\tc\t= 0.0\t# [K^-2]",
        "\td\t= 0.0\t# [K^-3]",
        "\te\t= 0.0\t# [K^-4]",
        "\tf\t= 0.0\t# [K^-1]",
        "\tg\t= 0.0\t# [K^-2]",
        "\th\t= 0.0\t# [K^-3]",
        "\ti\t= 0.0\t# [K^-4]",
        f"\tmm\t= {_num(mh)}\t# [1]",
        "}",
        "",
        "ConstantMobility:",
        "{ * mu_const = mumax (T/T0)^(-Exponent); Vela applies exponent 2.2 to both carriers.",
        f"\tmumax\t= {_num(si['mun'] * 1.0e4)} ,\t{_num(si['mup'] * 1.0e4)}\t# [cm^2/(Vs)]",
        "\tExponent\t= 2.2 ,\t2.2\t# [1]",
        "}",
        "",
        "}",
        "",
        f'Material = "{MATERIALS[1][1]}" {{',
        "",
        "Epsilon",
        "{ *  Ratio of the permittivities of material and vacuum",
        f"\tepsilon\t= {_num(sio2['eps_r'])}\t# [1]",
        "}",
        "",
        "Bandgap",
        "{ * Eg = Eg0 + dEg0 + alpha Tpar^2 / (beta + Tpar) - alpha T^2 / (beta + T)",
        f"\tChi0\t= {_num(sio2['electron_affinity_eV'])}\t# [eV]",
        f"\tEg0\t= {_num(sio2['bandgap_eV'])}\t# [eV]",
        "\talpha\t= 0.0\t# [eV K^-1]",
        "\tbeta\t= 636.0\t# [K]",
        f"\tTpar\t= {_num(REFERENCE_T)}\t# [K]",
        "}",
        "",
        "}",
        "",
    ]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--source", type=Path, default=MATERIAL_DATABASE)
    parser.add_argument("--check", action="store_true", help="exit 1 when the file is stale instead of writing")
    args = parser.parse_args(argv)
    text = render(builtin_materials(args.source))
    if args.check:
        current = args.output.read_text(encoding="utf-8").replace("\r\n", "\n") if args.output.exists() else None
        if current != text:
            print(f"{args.output} is stale; regenerate with scripts/generate_default_models_par.py", file=sys.stderr)
            return 1
        print(f"{args.output} is up to date")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
