#!/usr/bin/env python3
"""Generate ``share/vela/materials/models.par`` (the ``vela device -P`` default library).

The library is written in Sentaurus ``.par`` syntax and covers exactly the
sections the WP3 materials emitter implements.  Its numbers are the published
silicon defaults that Sentaurus Device also ships (and that
``tests/fixtures/sdevice/pn2d_sentaurus2018/models.par`` records), so a deck
whose parameter file omits a section falls back to the same physics a
Sentaurus run would use:

* ``Bandgap`` Varshni law and ``eDOSMass``/``hDOSMass`` Formula 1 masses
  (M. A. Green, J. Appl. Phys. 67, 2944 (1990));
* ``ConstantMobility`` / ``DopingDependence`` Masetti (Masetti, Severi, Solmi,
  IEEE TED 30, 764 (1983)); ``HighFieldDependence`` Canali (Canali et al.,
  IEEE TED 22, 1045 (1975));
* ``Scharfetter`` SRH lifetimes (Fossum & Lee), ``Auger`` (Dziewior & Schmid),
  ``OldSlotboom`` narrowing (Slotboom & de Graaff, Solid-State Electron. 19,
  857 (1976)), ``vanOverstraetendeMan`` ionization (van Overstraeten & de Man,
  Solid-State Electron. 13, 583 (1970)), density-gradient
  ``QuantumPotentialParameters`` (Wettstein, Schenk, Fichtner, IEEE TED 48,
  279 (2001)).

Where Vela cannot represent a Sentaurus default the library states Vela's own
law instead, so ``vela device -P`` never advertises a parameter that would be
rejected: ``ConstantMobility.Exponent`` is 2.2 for both carriers
(``Material::atTemperature``), ``mutunnel`` is omitted, ``Auger.H`` is 0 (no
carrier-density enhancement) and ``HighFieldDependence.betaexp``/``vsatexp``
are 0 (constant beta0/vsat0). The C++ ``MaterialDatabase`` Si/SiO2 built-ins
match the generated no-BGN 300 K snapshot in ``default_materials.json``.
Temperature laws and model activation remain separate from this alignment.

Usage::

    python scripts/generate_default_models_par.py            # rewrite the file
    python scripts/generate_default_models_par.py --check    # exit 1 if stale
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "share" / "vela" / "materials" / "models.par"
MATERIALS_OUTPUT = OUTPUT.with_name("default_materials.json")

#: (section header, comment, [(parameter, value or (electron, hole), unit)]).
#: Values are rendered verbatim so they round-trip bit-for-bit through the
#: parameter IR.
SILICON = [
    ("Epsilon", "Ratio of the permittivities of material and vacuum", [
        ("epsilon", "11.7", "1"),
    ]),
    ("Bandgap", "Eg = Eg0 + dEg0 + alpha Tpar^2/(beta + Tpar) - alpha T^2/(beta + T); "
                "dEg0(<model>) is folded into ni only when that narrowing model is active", [
        ("Chi0", "4.05", "eV"),
        ("Bgn2Chi", "0.5", "1"),
        ("Eg0", "1.16964", "eV"),
        ("dEg0(Bennett)", "0.0000e+00", "eV"),
        ("dEg0(Slotboom)", "-4.7950e-03", "eV"),
        ("dEg0(OldSlotboom)", "-1.5950e-02", "eV"),
        ("dEg0(delAlamo)", "-1.4070e-02", "eV"),
        ("alpha", "4.7300e-04", "eV K^-1"),
        ("beta", "6.3600e+02", "K"),
        ("alpha2", "0.0000e+00", "eV K^-1"),
        ("beta2", "0.0000e+00", "K"),
        ("EgMin", "0.0000e+00", "eV"),
        ("dEgMin", "0.01", "eV"),
        ("Tpar", "0.0000e+00", "K"),
    ]),
    ("OldSlotboom", "deltaEg = Ebgn (ln(N/Nref) + sqrt(ln(N/Nref)^2 + C))", [
        ("Ebgn", "9.0000e-03", "eV"),
        ("Nref", "1.0000e+17", "cm^(-3)"),
        ("C", "0.5", "1"),
    ]),
    ("eDOSMass", "Formula 1: me/m0 = [(6 mt)^2 ml]^(1/3) + mm, mt = a Eg(0)/Eg(T); "
                 "Nc(T) = 2 (2 pi kB T me/h^2)^(3/2)", [
        ("Formula", "1", "1"),
        ("a", "0.1905", "1"),
        ("ml", "0.9163", "1"),
        ("mm", "0.0000e+00", "1"),
    ]),
    ("hDOSMass", "Formula 1: mh/m0 = [(a+bT+cT^2+dT^3+eT^4)/(1+fT+gT^2+hT^3+iT^4)]^(2/3) + mm; "
                 "Nv(T) = 2 (2 pi kB T mh/h^2)^(3/2)", [
        ("Formula", "1", "1"),
        ("a", "0.443587", "1"),
        ("b", "0.003609528", "K^-1"),
        ("c", "0.0001173515", "K^-2"),
        ("d", "1.263218e-06", "K^-3"),
        ("e", "3.025581e-09", "K^-4"),
        ("f", "0.004683382", "K^-1"),
        ("g", "0.0002286895", "K^-2"),
        ("h", "7.469271e-07", "K^-3"),
        ("i", "1.727481e-09", "K^-4"),
        ("mm", "0", "1"),
    ]),
    ("QuantumPotentialParameters", "density-gradient quantum potential (electron component used by Vela)", [
        ("gamma", ("3.6", "5.6"), "1"),
        ("theta", ("0.5", "0.5"), "1"),
        ("xi", ("1", "1"), "1"),
        ("eta", ("1", "1"), "1"),
        ("nu", ("0.0000e+00", "0.0000e+00"), "1"),
    ]),
    ("ConstantMobility:", "mu_const = mumax (T/300 K)^(-Exponent); Vela applies exponent 2.2 to both carriers "
                          "(the Sentaurus electron default 2.5 is not representable)", [
        ("mumax", ("1.4170e+03", "4.7050e+02"), "cm^2/(Vs)"),
        ("Exponent", ("2.2", "2.2"), "1"),
    ]),
    ("DopingDependence:", "Formula 1 (Masetti): mu = mumin1 exp(-Pc/N) + (mu_const - mumin2)/(1+(N/Cr)^alpha) "
                          "- mu1/(1+(Cs/N)^beta)", [
        ("formula", ("1", "1"), "1"),
        ("mumin1", ("52.2", "44.9"), "cm^2/Vs"),
        ("mumin2", ("52.2", "0.0000e+00"), "cm^2/Vs"),
        ("mu1", ("43.4", "29"), "cm^2/Vs"),
        ("Pc", ("0.0000e+00", "9.2300e+16"), "cm^3"),
        ("Cr", ("9.6800e+16", "2.2300e+17"), "cm^3"),
        ("Cs", ("3.4300e+20", "6.1000e+20"), "cm^3"),
        ("alpha", ("0.68", "0.719"), "1"),
        ("beta", ("2", "2"), "1"),
    ]),
    ("HighFieldDependence:", "Canali: mu(F) = mu_low/(1 + (mu_low F/vsat)^beta)^(1/beta); Vela keeps beta0 and vsat0 "
                             "constant in T, hence betaexp = vsatexp = 0 (Sentaurus defaults 0.66, 0.17 and 0.87, 0.52)", [
        ("beta0", ("1.109", "1.213"), "1"),
        ("betaexp", ("0", "0"), "1"),
        ("alpha", ("0.0000e+00", "0.0000e+00"), "1"),
        ("K_dT", ("0.2", "0.2"), "1"),
        ("E0_TrEf", ("4.0000e+03", "4.0000e+03"), "1"),
        ("Ksmooth_TrEf", ("1", "1"), "1"),
        ("Vsat_Formula", ("1", "1"), "1"),
        ("vsat0", ("1.0700e+07", "8.3700e+06"), "1"),
        ("vsatexp", ("0", "0"), "1"),
        ("ku", ("1", "1"), "1"),
        ("kv", ("1", "1"), "1"),
    ]),
    ("Scharfetter", "tau = taumin + (taumax - taumin)/(1 + (N/Nref)^gamma); tau(T) = tau (T/300 K)^Talpha", [
        ("taumin", ("0.0000e+00", "0.0000e+00"), "s"),
        ("taumax", ("1.0000e-05", "3.0000e-06"), "s"),
        ("Nref", ("1.0000e+16", "1.0000e+16"), "cm^(-3)"),
        ("gamma", ("1", "1"), "1"),
        ("Talpha", ("-1.5000e+00", "-1.5000e+00"), "1"),
        ("Tcoeff", ("2.55", "2.55"), "1"),
        ("Etrap", "0.0000e+00", "eV"),
        ("tau0", ("1.0000e-05", "3.0000e-06"), "s"),
    ]),
    ("Auger", "R = (Cn n + Cp p)(np - ni^2), Cn,p = A + B (T/300) + C (T/300)^2; Vela has no "
              "(1 + H exp(-n/N0)) enhancement, hence H = 0 (Sentaurus default H = 3.46667, 8.25688)", [
        ("A", ("6.7000e-32", "7.2000e-32"), "cm^6/s"),
        ("B", ("2.4500e-31", "4.5000e-33"), "cm^6/s"),
        ("C", ("-2.2000e-32", "2.6300e-32"), "cm^6/s"),
        ("H", ("0", "0"), "1"),
        ("N0", ("1.0000e+18", "1.0000e+18"), "cm^(-3)"),
    ]),
    ("vanOverstraetendeMan", "alpha(F) = gamma a exp(-gamma b/F), gamma = tanh(hbarOmega/2kT0)/tanh(hbarOmega/2kT); "
                             "low/high branches switch at E0", [
        ("a(low)", ("7.0300e+05", "1.5820e+06"), "1/cm"),
        ("a(high)", ("7.0300e+05", "6.7100e+05"), "1/cm"),
        ("b(low)", ("1.2310e+06", "2.0360e+06"), "V/cm"),
        ("b(high)", ("1.2310e+06", "1.6930e+06"), "V/cm"),
        ("E0", ("4.0000e+05", "4.0000e+05"), "V/cm"),
        ("hbarOmega", ("0.063", "0.063"), "eV"),
        ("beta(low)", ("0.678925", "0.815009"), "1"),
        ("beta(high)", ("0.678925", "0.677706"), "1"),
        ("lambda", ("6.2000e-07", "4.5000e-07"), "cm"),
    ]),
]

OXIDE = [
    ("Epsilon", "Ratio of the permittivities of material and vacuum", [
        ("epsilon", "3.9", "1"),
    ]),
    ("Bandgap", "temperature-independent gap (alpha = 0); band-edge reference only", [
        ("Chi0", "0.9", "eV"),
        ("Eg0", "9", "eV"),
        ("alpha", "0.0000e+00", "eV K^-1"),
        ("beta", "6.3600e+02", "K"),
        ("Tpar", "3.0000e+02", "K"),
    ]),
]

NITRIDE = [
    ("Epsilon", "Ratio of the permittivities of material and vacuum", [
        ("epsilon", "7.5", "1"),
    ]),
    ("Bandgap", "temperature-independent gap (alpha = 0); band-edge reference only", [
        ("Chi0", "1.9", "eV"),
        ("Eg0", "5", "eV"),
        ("alpha", "0.0000e+00", "eV K^-1"),
        ("beta", "6.3600e+02", "K"),
        ("Tpar", "3.0000e+02", "K"),
    ]),
]

MATERIALS = [("Silicon", SILICON), ("Oxide", OXIDE), ("Nitride", NITRIDE)]

HEADER = """\
* Vela default material parameters for `vela device -P` (generated; do not edit by hand).
* Rendered by scripts/generate_default_models_par.py in Sentaurus .par syntax, restricted to the
* sections the Vela materials emitter implements. Silicon uses the published defaults that
* Sentaurus Device also ships (Green 1990 band gap and DOS masses, Masetti 1983 mobility,
* Canali 1975 high-field saturation, Fossum SRH lifetimes, Dziewior-Schmid Auger,
* Slotboom-de Graaff 1976 narrowing, van Overstraeten-de Man 1970 ionization).
* Where Vela cannot represent the Sentaurus default the value describes Vela's own law:
* ConstantMobility.Exponent = 2.2 for both carriers, no mutunnel, Auger.H = 0,
* HighFieldDependence.betaexp = vsatexp = 0.
* The intrinsic density is always derived: ni = sqrt(Nc Nv) exp(-(Eg(T) + dEg0)/(2 kT)), with
* dEg0(OldSlotboom) folded in only when the OldSlotboom model is active.
"""


def _line(name: str, value, unit: str) -> str:
    text = f"{value[0]} ,\t{value[1]}" if isinstance(value, tuple) else str(value)
    return f"\t{name}\t= {text}\t# [{unit}]"


def render() -> str:
    lines = [HEADER]
    for material, sections in MATERIALS:
        lines.append(f'Material = "{material}" {{')
        lines.append("")
        for section, comment, parameters in sections:
            lines.append(section)
            lines.append(f"{{ * {comment}")
            lines.extend(_line(name, value, unit) for name, value, unit in parameters)
            lines.append("}")
            lines.append("")
        lines.append("}")
        lines.append("")
    return "\n".join(lines)


def render_materials(par_path: Path) -> str:
    """Freeze the emitter's Si/SiO2 scalars at 300 K without active BGN."""
    import json
    from vela_sdevice.materials import emit_materials
    document = emit_materials(par_path, active_models=[], temperature_K=300.0,
                              mesh_materials=["Si", "SiO2"], benchmark="builtin:300K").materials
    return json.dumps(document, indent=2, ensure_ascii=False) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--materials-output", type=Path,
                        help="300 K JSON snapshot (defaults beside --output)")
    parser.add_argument("--check", action="store_true", help="exit 1 when the file is stale instead of writing")
    args = parser.parse_args(argv)
    materials_output = args.materials_output or args.output.with_name("default_materials.json")
    text = render()
    if args.check:
        current = args.output.read_text(encoding="utf-8").replace("\r\n", "\n") if args.output.exists() else None
        if current != text:
            print(f"{args.output} is stale; regenerate with scripts/generate_default_models_par.py", file=sys.stderr)
            return 1
        materials_text = render_materials(args.output)
        current_materials = materials_output.read_text(encoding="utf-8") if materials_output.exists() else None
        if current_materials != materials_text:
            print(f"{materials_output} is stale; regenerate with scripts/generate_default_models_par.py", file=sys.stderr)
            return 1
        print(f"{args.output} is up to date")
        print(f"{materials_output} is up to date")
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    print(f"wrote {args.output}")
    materials_output.parent.mkdir(parents=True, exist_ok=True)
    with materials_output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(render_materials(args.output))
    print(f"wrote {materials_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
