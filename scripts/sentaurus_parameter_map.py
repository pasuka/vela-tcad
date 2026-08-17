#!/usr/bin/env python3
"""Mapping contract between Sentaurus ``.par`` parameters and Vela inputs.

This module is the *contract* layer of the ``.par`` import pipeline.  It sits
between the lossless syntax IR produced by :mod:`sentaurus_parameter_ir` and
any emitter that writes Vela JSON, and it answers exactly one question:

    May this ``.par`` parameter be turned into a Vela input, and with what
    loss of fidelity?

The presence of a Vela JSON field is *not* evidence that a Sentaurus parameter
can be imported.  Two models can share a field name and still solve different
equations.  Every mapping therefore carries an explicit status:

``exact``
    Vela evaluates the same formula from the same parameters.  The imported
    value reproduces Sentaurus at every temperature the formula covers.

``frozen_at_temperature``
    The value is numerically correct only at one temperature, because Vela
    either has no temperature law for it or uses a different one.  The
    emitter must record the temperature it froze at.

``approximated``
    Vela implements a reduced form of the Sentaurus model.  Importing loses
    physics, so it requires an explicit opt-in.

``unsupported_formula``
    The section exists in Vela but this particular Sentaurus ``Formula`` or
    named sub-model is not implemented.

``unsupported_model``
    Vela has no counterpart at all.

Only ``exact`` and ``frozen_at_temperature`` are importable by default.

The matrix deliberately does *not* describe every section in a ``.par`` file.
A ``.par`` file is a candidate library, not an activation list: entries here
are the parameters an emitter is allowed to consider once the corresponding
model has been shown to be active by the SDevice execution IR.  Sections that
are absent from the matrix are reported as ``unmapped`` and, like
``unsupported_*``, block the import unless the caller opts in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

STATUS_EXACT = "exact"
STATUS_FROZEN = "frozen_at_temperature"
STATUS_APPROXIMATED = "approximated"
STATUS_UNSUPPORTED_FORMULA = "unsupported_formula"
STATUS_UNSUPPORTED_MODEL = "unsupported_model"

#: Every status an entry may carry, ordered from most to least faithful.
STATUSES = (
    STATUS_EXACT,
    STATUS_FROZEN,
    STATUS_APPROXIMATED,
    STATUS_UNSUPPORTED_FORMULA,
    STATUS_UNSUPPORTED_MODEL,
)

#: Status used for a parameter that the matrix does not describe at all.
STATUS_UNMAPPED = "unmapped"
STATUS_INACTIVE = "inactive"

MODEL_ALIASES = {
    "DopingDep": "DopingDependence",
    "HighFieldsaturation": "HighFieldSaturation",
    "eHighFieldSaturation": "HighFieldSaturation",
    "eHighFieldsaturation": "HighFieldSaturation",
    "hHighFieldSaturation": "HighFieldSaturation",
    "hHighFieldsaturation": "HighFieldSaturation",
}

SECTION_MODELS = {
    "DopingDependence": "DopingDependence",
    "HighFieldDependence": "HighFieldSaturation",
    "Scharfetter": "SRH",
    "Auger": "Auger",
    "OldSlotboom": "OldSlotboom",
    "vanOverstraetendeMan": "VanOverstraeten",
    "QuantumPotentialParameters": "eQuantumPotential",
}

#: Statuses an emitter may act on without an explicit lossy opt-in.
DEFAULT_ALLOWED_STATUSES = frozenset({STATUS_EXACT, STATUS_FROZEN})

#: Statuses an emitter may act on with ``--allow-lossy``.
LOSSY_ALLOWED_STATUSES = DEFAULT_ALLOWED_STATUSES | {STATUS_APPROXIMATED}


class ParameterMapError(ValueError):
    """Raised when a ``.par`` parameter may not be imported."""


@dataclass(frozen=True)
class MappingEntry:
    """One row of the mapping matrix.

    ``section`` and ``parameter`` are the *raw* Sentaurus spellings as they
    appear in the syntax IR, so the matrix can be joined against IR output
    without any normalisation step that could silently merge distinct names.
    """

    section: str
    parameter: str
    status: str
    #: Dotted path of the Vela input this feeds, or ``None`` when the mapping
    #: is not a direct field write (for example a term of a derived quantity).
    target: Optional[str] = None
    #: Paren variant that qualifies the parameter, e.g. ``dEg0(OldSlotboom)``.
    variant: Optional[str] = None
    #: Sentaurus ``Formula`` value this row is valid for, when the section is
    #: formula-switched.  ``None`` means the row is formula-independent.
    formula: Optional[str] = None
    #: Execution-IR model that must be active before this row may be used.
    #: ``None`` marks a material-level constant that is always meaningful.
    requires_model: Optional[str] = None
    #: Value at which the parameter is inert, i.e. it has no effect on the
    #: solved equations.  A parameter Vela cannot represent is harmless when
    #: Sentaurus is also not using it, so an unsupported row whose observed
    #: value is neutral is demoted to ``exact``.  ``None`` means the parameter
    #: has no known neutral value and any occurrence is meaningful.
    neutral_value: Optional[float] = None
    #: Why the status is what it is.  Every non-``exact`` row must say so.
    note: str = ""

    @property
    def key(self) -> tuple:
        return (self.section, self.parameter, self.variant, self.formula)


def _rows() -> list[MappingEntry]:
    """Build the matrix.

    Grouped by the Vela consumer rather than by ``.par`` file order, because
    the reviewer of this contract cares about "what feeds Vela field X".
    """
    rows: list[MappingEntry] = []

    # -- Material scalars ------------------------------------------------
    # Vela stores a single relative permittivity; Sentaurus Epsilon is the
    # same isotropic quantity for the materials in the corpus.
    rows.append(MappingEntry(
        section="Epsilon", parameter="epsilon",
        status=STATUS_EXACT, target="materials[].eps_r",
    ))

    # Varshni: Eg(T) = Eg0 - alpha T^2 / (T + beta), referred to Tpar.
    # Vela has no Varshni evaluator, so the emitter must evaluate it once and
    # write a single number: correct at the chosen temperature only.
    for name in ("Eg0", "alpha", "beta", "Tpar"):
        rows.append(MappingEntry(
            section="Bandgap", parameter=name,
            status=STATUS_FROZEN, target="materials[].bandgap_eV",
            note=(
                "Vela stores a scalar bandgap_eV and has no Varshni law; the "
                "emitter evaluates Varshni at one temperature. Tpar is the "
                "reference temperature of Eg0 and is 0 K in the corpus, so "
                "Eg0 is not the 300 K value"
            ),
        ))
    rows.append(MappingEntry(
        section="Bandgap", parameter="Chi0",
        status=STATUS_EXACT, target="materials[].electron_affinity_eV",
    ))
    # dEg0 is the bandgap-narrowing zero-doping reference, not a bandgap
    # shift.  Vela folds it into the material ni, and writing it into
    # bandgap_eV as well would double-count it.
    rows.append(MappingEntry(
        section="Bandgap", parameter="dEg0", variant="OldSlotboom",
        status=STATUS_EXACT, target="materials[].ni",
        requires_model="OldSlotboom",
        note=(
            "folded into the material ni, never into bandgap_eV; Vela keeps "
            "the base Varshni gap and the BGN reference separate"
        ),
    ))
    # The other dEg0 references are inert until SDevice selects them; each is
    # gated on its own model so an unused library entry never blocks.
    for variant in ("Slotboom", "Bennett", "delAlamo"):
        rows.append(MappingEntry(
            section="Bandgap", parameter="dEg0", variant=variant,
            status=STATUS_UNSUPPORTED_FORMULA, requires_model=variant,
            note="Vela implements the OldSlotboom reference only",
        ))

    # DOS masses feed Nc/Nv, and Nc/Nv together with the gap feed ni. Both
    # sections are Formula-switched; only Formula 1 is derivable here.
    for name in ("a", "ml", "mm"):
        rows.append(MappingEntry(
            section="eDOSMass", parameter=name,
            status=STATUS_FROZEN, target="materials[].Nc_m3", formula="1",
            note=(
                "Vela stores a scalar Nc_m3; the emitter evaluates the "
                "Formula 1 mass law and Nc(T) at one temperature"
            ),
        ))
    for name in ("a", "b", "c", "d", "e", "f", "g", "h", "i", "mm"):
        rows.append(MappingEntry(
            section="hDOSMass", parameter=name,
            status=STATUS_FROZEN, target="materials[].Nv_m3", formula="1",
            note=(
                "Vela stores a scalar Nv_m3; the emitter evaluates the "
                "Formula 1 mass law and Nv(T) at one temperature"
            ),
        ))
    for section in ("eDOSMass", "hDOSMass"):
        rows.append(MappingEntry(
            section=section, parameter="Formula",
            status=STATUS_EXACT, target=None,
            note="selector consumed by the emitter, not written to Vela",
        ))
        rows.append(MappingEntry(
            section=section, parameter="Nc300" if section == "eDOSMass" else "Nv300",
            status=STATUS_UNSUPPORTED_FORMULA, formula="2",
            note="Formula 2 (direct Nc300/Nv300) is not wired into the emitter",
        ))

    # -- Mobility --------------------------------------------------------
    # ConstantMobility mumax is the lattice mobility at Tpar.
    rows.append(MappingEntry(
        section="ConstantMobility", parameter="mumax",
        status=STATUS_FROZEN, target="materials[].mun / materials[].mup",
        note=(
            "Vela scales both carriers by the same hardcoded (T/T0)^-2.2, "
            "while .par carries separate electron/hole exponents; only the "
            "value at the reference temperature is trustworthy"
        ),
    ))
    rows.append(MappingEntry(
        section="ConstantMobility", parameter="Exponent",
        status=STATUS_APPROXIMATED, target=None,
        note=(
            "Material::atTemperature applies pow(ratio, -2.2) to electrons "
            "and holes alike; the per-carrier exponents (2.5, 2.2) cannot be "
            "represented"
        ),
    ))

    # DopingDependence Formula 1 is Masetti, which Vela implements term for
    # term.  Formula 2 (Arora) is not implemented.
    masetti = {
        "mumin1": "solver.mobility.electron_mumin1_m2_V_s / hole_mumin1_m2_V_s",
        "mumin2": "solver.mobility.electron_mumin2_m2_V_s / hole_mumin2_m2_V_s",
        "mu1": "solver.mobility.electron_mu1_m2_V_s / hole_mu1_m2_V_s",
        "Pc": "solver.mobility.electron_pc_m3 / hole_pc_m3",
        "Cr": "solver.mobility.electron_cr_m3 / hole_cr_m3",
        "Cs": "solver.mobility.electron_cs_m3 / hole_cs_m3",
        "alpha": "solver.mobility.electron_masetti_alpha / hole_masetti_alpha",
        "beta": "solver.mobility.electron_masetti_beta / hole_masetti_beta",
    }
    for name, target in masetti.items():
        rows.append(MappingEntry(
            section="DopingDependence", parameter=name,
            status=STATUS_EXACT, target=target, formula="1",
            requires_model="DopingDependence",
        ))
    rows.append(MappingEntry(
        section="DopingDependence", parameter="Formula",
        status=STATUS_EXACT, target=None,
        requires_model="DopingDependence",
        note="selector consumed by the emitter, not written to Vela",
    ))
    for name in ("Ar_mumin", "Ar_mumax", "Ar_N0", "Ar_alpha",
                 "Ar_expt0", "Ar_expt1", "Ar_expt2", "Ar_expt3"):
        rows.append(MappingEntry(
            section="DopingDependence", parameter=name,
            status=STATUS_UNSUPPORTED_FORMULA, formula="2",
            requires_model="DopingDependence",
            note="Arora (Formula 2) doping dependence is not implemented",
        ))

    # HighFieldDependence: Vela has the Caughey-Thomas shape but no
    # temperature laws for beta or vsat, and no Transferred-Electron terms.
    rows.append(MappingEntry(
        section="HighFieldDependence", parameter="beta0",
        status=STATUS_FROZEN,
        target="solver.mobility.electron_field_beta / hole_field_beta",
        requires_model="HighFieldSaturation",
        note="Vela has no betaexp temperature law; valid at the frozen T",
    ))
    rows.append(MappingEntry(
        section="HighFieldDependence", parameter="vsat0",
        status=STATUS_FROZEN,
        target=("solver.mobility.electron_saturation_velocity_m_s / "
                "hole_saturation_velocity_m_s"),
        requires_model="HighFieldSaturation",
        note="Vela has no vsatexp temperature law; valid at the frozen T",
    ))
    for name in ("betaexp", "vsatexp"):
        rows.append(MappingEntry(
            section="HighFieldDependence", parameter=name,
            status=STATUS_APPROXIMATED, target=None,
            requires_model="HighFieldSaturation",
            note=(
                "temperature exponent has no Vela counterpart; it is absorbed "
                "by freezing beta0/vsat0 and is otherwise dropped"
            ),
        ))
    # Transferred-Electron and driving-force smoothing terms.  Each is inert
    # at the value the corpus ships, so dropping them is lossless there; the
    # neutral value is recorded so a file that changes them fails closed.
    for name, neutral, model in (
        ("alpha", 0.0, "HighFieldSaturation"),
        ("ku", 1.0, "HighFieldSaturation"),
        ("kv", 1.0, "HighFieldSaturation"),
        ("K_dT", None, "TransferredElectron"),
        ("E0_TrEf", None, "TransferredElectron"),
        ("Ksmooth_TrEf", 1.0, "TransferredElectron"),
    ):
        rows.append(MappingEntry(
            section="HighFieldDependence", parameter=name,
            status=STATUS_UNSUPPORTED_MODEL, target=None,
            requires_model=model,
            neutral_value=neutral,
            note=(
                "Transferred-Electron / driving-force smoothing term has no "
                "Vela counterpart"
            ),
        ))
    rows.append(MappingEntry(
        section="HighFieldDependence", parameter="Vsat_Formula",
        status=STATUS_UNSUPPORTED_FORMULA, target=None, formula=None,
        requires_model="HighFieldSaturation",
        neutral_value=1.0,
        note=(
            "saturation-velocity formula selector; only Formula 1 matches the "
            "Caughey-Thomas shape Vela implements"
        ),
    ))

    # -- Recombination ---------------------------------------------------
    scharfetter = {
        "tau0": "solver.srh_doping_dependence.{carrier}.tau_max_s",
        "taumin": "solver.srh_doping_dependence.{carrier}.tau_min_s",
        "taumax": "solver.srh_doping_dependence.{carrier}.tau_max_s",
        "Nref": "solver.srh_doping_dependence.{carrier}.reference_doping_m3",
        "gamma": "solver.srh_doping_dependence.{carrier}.gamma",
    }
    for name, target in scharfetter.items():
        rows.append(MappingEntry(
            section="Scharfetter", parameter=name,
            status=STATUS_EXACT, target=target, requires_model="SRH",
        ))
    rows.append(MappingEntry(
        section="Scharfetter", parameter="Talpha",
        status=STATUS_EXACT,
        target=("solver.srh_doping_dependence.electron_temperature_exponent / "
                "hole_temperature_exponent"),
        requires_model="TempDependence",
        note="power-law lifetime temperature dependence (TempDep)",
    ))
    rows.append(MappingEntry(
        section="Scharfetter", parameter="Tcoeff",
        status=STATUS_UNSUPPORTED_FORMULA, target=None,
        requires_model="ExpTempDependence",
        note="ExpTempDep exponential lifetime law is not implemented",
    ))
    rows.append(MappingEntry(
        section="Scharfetter", parameter="Etrap",
        status=STATUS_UNSUPPORTED_MODEL, target=None, neutral_value=0.0,
        requires_model="SRH",
        note="Vela fixes the SRH trap at mid-gap",
    ))

    for name in ("A", "B", "C"):
        rows.append(MappingEntry(
            section="Auger", parameter=name,
            status=STATUS_APPROXIMATED,
            target="solver.auger_cn_m6_per_s / solver.auger_cp_m6_per_s",
            requires_model="Auger",
            note=(
                "Vela Auger has constant Cn/Cp only; the A/B/C temperature "
                "polynomial collapses to its value at the frozen temperature"
            ),
        ))
    for name in ("H", "N0"):
        rows.append(MappingEntry(
            section="Auger", parameter=name,
            status=STATUS_UNSUPPORTED_MODEL, target=None,
            requires_model="Auger",
            note="carrier-density enhancement term is not implemented",
        ))

    # -- Bandgap narrowing ----------------------------------------------
    rows.append(MappingEntry(
        section="OldSlotboom", parameter="Ebgn",
        status=STATUS_EXACT,
        target="solver.bandgap_narrowing.coefficient_eV",
        requires_model="OldSlotboom",
    ))
    rows.append(MappingEntry(
        section="OldSlotboom", parameter="Nref",
        status=STATUS_EXACT,
        target="solver.bandgap_narrowing.reference_doping_m3",
        requires_model="OldSlotboom",
    ))
    rows.append(MappingEntry(
        section="OldSlotboom", parameter="C",
        status=STATUS_EXACT,
        target="solver.bandgap_narrowing.smoothing",
        requires_model="OldSlotboom",
    ))

    # -- Impact ionization ----------------------------------------------
    overstraeten = [
        ("a", "low", "solver.impact_ionization.electron_a_low_m_inv / hole_a_low_m_inv"),
        ("a", "high", "solver.impact_ionization.electron_a_high_m_inv / hole_a_high_m_inv"),
        ("b", "low", "solver.impact_ionization.electron_b_low_V_m / hole_b_low_V_m"),
        ("b", "high", "solver.impact_ionization.electron_b_high_V_m / hole_b_high_V_m"),
        ("E0", None, "solver.impact_ionization.switch_field_V_m"),
        ("hbarOmega", None, "solver.impact_ionization.phonon_energy_eV"),
    ]
    for name, variant, target in overstraeten:
        rows.append(MappingEntry(
            section="vanOverstraetendeMan", parameter=name,
            status=STATUS_EXACT, target=target, variant=variant,
            requires_model="VanOverstraeten",
        ))
    for name, variant in (("beta", "low"), ("beta", "high"), ("lambda", None)):
        rows.append(MappingEntry(
            section="vanOverstraetendeMan", parameter=name,
            status=STATUS_UNSUPPORTED_FORMULA, target=None, variant=variant,
            requires_model="BandgapDependence",
            note=(
                "BandgapDependence critical field Ecrit = beta*Eg/(q*lambda) "
                "is not implemented; Vela always uses Ecrit = b"
            ),
        ))
    for section, model in (
        ("OkutoCrowell", "OkutoCrowell"),
        ("Lackner", "Lackner"),
        ("UniBo", "UniBo"),
        ("UniBo2", "UniBo2"),
        ("Hatakeyama", "Hatakeyama"),
        ("vanOverstraetendeMan_aniso", "VanOverstraetenAniso"),
    ):
        rows.append(MappingEntry(
            section=section, parameter="*",
            status=STATUS_UNSUPPORTED_MODEL, target=None,
            requires_model=model,
            note=(
                "ImpactIonizationModel implements selberherr and "
                "van_overstraeten only"
            ),
        ))

    # -- Quantum correction ---------------------------------------------
    rows.append(MappingEntry(
        section="QuantumPotentialParameters", parameter="gamma",
        status=STATUS_EXACT,
        target="materials[].electron_quantum_gamma",
        requires_model="eQuantumPotential",
        note="electron component only; the hole component has no Vela path",
    ))
    rows.append(MappingEntry(
        section="QuantumPotentialParameters", parameter="theta",
        status=STATUS_EXACT,
        target="solver.electron_quantum_potential.theta",
        requires_model="eQuantumPotential",
        note="electron component only; Vela has a single scalar theta",
    ))
    for name, neutral in (("xi", 1.0), ("eta", 1.0), ("nu", 0.0)):
        rows.append(MappingEntry(
            section="QuantumPotentialParameters", parameter=name,
            status=STATUS_UNSUPPORTED_MODEL, target=None,
            requires_model="eQuantumPotential",
            neutral_value=neutral,
            note=(
                "Vela's Eq. 231 density-gradient form has no adjustable "
                "quasi-Fermi, electrostatic, or stress weight"
            ),
        ))

    return rows


#: The mapping matrix, as an ordered list of rows.
MATRIX: tuple[MappingEntry, ...] = tuple(_rows())


def _build_index() -> dict[tuple[str, str], list[MappingEntry]]:
    index: dict[tuple[str, str], list[MappingEntry]] = {}
    for entry in MATRIX:
        index.setdefault((entry.section, entry.parameter), []).append(entry)
    return index


_INDEX = _build_index()


def lookup(
    section: str,
    parameter: str,
    variant: Optional[str] = None,
    formula: Optional[str] = None,
) -> Optional[MappingEntry]:
    """Return the matrix row for a parameter, or ``None`` when unmapped.

    Matching is intentionally strict on ``section`` and ``parameter`` so a
    renamed or misspelled parameter can never silently borrow another row's
    status.  ``variant`` and ``formula`` narrow the match: a row that pins a
    variant or formula only matches when the caller supplies the same one,
    while a row that leaves them open matches any caller value.
    """
    if parameter == "formula":
        parameter = "Formula"
    candidates = _INDEX.get((section, parameter))
    if candidates is None:
        candidates = _INDEX.get((section, "*"))
    if not candidates:
        return None

    def score(entry: MappingEntry) -> Optional[int]:
        points = 0
        if entry.variant is not None:
            if entry.variant != variant:
                return None
            points += 2
        elif variant is not None:
            # The caller has a variant the matrix does not describe.  A row
            # that says nothing about variants must not claim it.
            return None
        if entry.formula is not None:
            if entry.formula != formula:
                return None
            points += 1
        return points

    best: Optional[MappingEntry] = None
    best_points = -1
    for entry in candidates:
        points = score(entry)
        if points is not None and points > best_points:
            best, best_points = entry, points
    return best


@dataclass
class Classification:
    """Result of checking one parameter against the matrix."""

    section: str
    parameter: str
    variant: Optional[str]
    formula: Optional[str]
    status: str
    entry: Optional[MappingEntry]
    reason: str = ""

    @property
    def importable_by_default(self) -> bool:
        return self.status in DEFAULT_ALLOWED_STATUSES

    def to_json(self) -> dict:
        return {
            "section": self.section,
            "parameter": self.parameter,
            "variant": self.variant,
            "formula": self.formula,
            "status": self.status,
            "target": self.entry.target if self.entry else None,
            "requires_model": self.entry.requires_model if self.entry else None,
            "reason": self.reason,
        }


_NEUTRAL_TOLERANCE = 1.0e-12


def _is_inert(entry: MappingEntry, values: Optional[Iterable]) -> bool:
    """Return whether an observed value renders the parameter inert.

    Both carriers must sit at the neutral value: a term that is switched off
    for electrons but active for holes still changes the solve.
    """
    if entry.neutral_value is None or values is None:
        return False
    observed = list(values)
    if not observed:
        return False
    for value in observed:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            return False
        if abs(float(value) - entry.neutral_value) > _NEUTRAL_TOLERANCE:
            return False
    return True


def classify(
    section: str,
    parameter: str,
    variant: Optional[str] = None,
    formula: Optional[str] = None,
    active_models: Optional[Iterable[str]] = None,
    values: Optional[Iterable] = None,
) -> Classification:
    """Classify one parameter, honouring the activated-model context.

    ``active_models`` is the set of models the SDevice execution IR proved to
    be switched on.  A ``.par`` file is a candidate library, so a parameter
    whose model is not active is *not* an error: it is simply inert, and is
    reported as such rather than being imported.  Passing ``None`` means "no
    activation context available", which downgrades every model-gated row.

    ``values`` is the parsed value from the syntax IR.  It is used only to
    detect parameters that are switched off in this particular file.
    """
    models = None if active_models is None else {
        MODEL_ALIASES.get(model, model) for model in active_models
    }
    entry = lookup(section, parameter, variant, formula)
    if entry is None:
        section_model = SECTION_MODELS.get(section, section)
        if models is None or section_model not in models:
            return Classification(
                section, parameter, variant, formula, STATUS_INACTIVE, None,
                reason=f"section model {section_model!r} is not active",
            )
        return Classification(
            section, parameter, variant, formula, STATUS_UNMAPPED, None,
            reason="no mapping contract row describes this parameter",
        )

    if entry.requires_model is not None:
        if models is None:
            return Classification(
                section, parameter, variant, formula,
                STATUS_INACTIVE, entry,
                reason=(
                    f"requires model {entry.requires_model!r} but no activated"
                    " model context was supplied"
                ),
            )
        required = MODEL_ALIASES.get(entry.requires_model, entry.requires_model)
        if required not in models:
            return Classification(
                section, parameter, variant, formula,
                STATUS_INACTIVE, entry,
                reason=(
                    f"model {entry.requires_model!r} is not active; the "
                    "parameter is inert and must not enable a new model"
                ),
            )

    if entry.status not in DEFAULT_ALLOWED_STATUSES and _is_inert(entry, values):
        return Classification(
            section, parameter, variant, formula, STATUS_EXACT, entry,
            reason=(
                f"{entry.status} in general, but this file sets it to its "
                f"neutral value {entry.neutral_value:g}, so dropping it "
                "changes nothing"
            ),
        )

    return Classification(
        section, parameter, variant, formula, entry.status, entry,
        reason=entry.note,
    )


@dataclass
class CoverageReport:
    """Per-status tallies plus the rows that block a default import."""

    counts: dict = field(default_factory=dict)
    blocking: list = field(default_factory=list)
    lossy: list = field(default_factory=list)

    def to_json(self) -> dict:
        return {
            "counts": dict(sorted(self.counts.items())),
            "blocking": [item.to_json() for item in self.blocking],
            "lossy": [item.to_json() for item in self.lossy],
        }


def _formula_of(block: dict) -> Optional[str]:
    """Return the ``Formula`` selector of a block as a canonical string.

    Sentaurus writes ``Formula = 1``, which the syntax IR carries as the float
    ``1.0``.  Rendering that back with plain ``str`` would produce ``"1.0"``
    and never match a matrix row, so integral values are narrowed explicitly.
    """
    for parameter in block.get("parameters", []):
        if str(parameter.get("base_name", "")).lower() != "formula":
            continue
        values = parameter.get("values") or []
        if not values:
            return None
        value = values[0]
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value)
    return None


def classify_ir(
    ir: dict,
    active_models: Optional[Iterable[str]] = None,
) -> tuple[list[Classification], CoverageReport]:
    """Classify every parameter of a ``sentaurus_parameter_ir`` document.

    Only the *active* blocks are considered, so a section that a later
    definition shadows does not inflate the coverage numbers.
    """
    results: list[Classification] = []
    report = CoverageReport()
    models = active_models

    for block in ir.get("blocks", []):
        if block.get("shadowed_by") is not None:
            continue
        section = block.get("section", "")
        formula = _formula_of(block)
        for parameter in block.get("parameters", []):
            if parameter.get("shadowed_by") is not None:
                continue
            result = classify(
                section,
                parameter.get("base_name", parameter.get("raw_name", "")),
                parameter.get("variant"),
                formula,
                models,
                parameter.get("values"),
            )
            results.append(result)
            report.counts[result.status] = report.counts.get(result.status, 0) + 1
            if result.status == STATUS_APPROXIMATED:
                report.lossy.append(result)
            elif (result.status not in DEFAULT_ALLOWED_STATUSES
                    and result.status != STATUS_INACTIVE):
                report.blocking.append(result)

    return results, report


def assert_importable(
    results: Iterable[Classification],
    allow_lossy: bool = False,
) -> None:
    """Fail closed unless every classification is importable.

    Only ``inactive`` results are exempt. Unknown parameters in an active
    section are fatal because silently dropping them violates fail-closed
    import semantics.
    """
    allowed = LOSSY_ALLOWED_STATUSES if allow_lossy else DEFAULT_ALLOWED_STATUSES
    problems = [
        item for item in results
        if item.status not in allowed and item.status != STATUS_INACTIVE
    ]
    if not problems:
        return
    lines = [
        f"{len(problems)} .par parameter(s) cannot be imported faithfully:"
    ]
    for item in problems[:20]:
        detail = f" ({item.reason})" if item.reason else ""
        lines.append(
            f"  {item.section}.{item.parameter}: {item.status}{detail}"
        )
    if len(problems) > 20:
        lines.append(f"  ... and {len(problems) - 20} more")
    if not allow_lossy and any(
        item.status == STATUS_APPROXIMATED for item in problems
    ):
        lines.append("re-run with allow_lossy=True to accept approximations")
    raise ParameterMapError("\n".join(lines))


def matrix_as_json() -> list[dict]:
    """Serialise the matrix so it can be diffed and reviewed as data."""
    return [
        {
            "section": entry.section,
            "parameter": entry.parameter,
            "variant": entry.variant,
            "formula": entry.formula,
            "status": entry.status,
            "target": entry.target,
            "requires_model": entry.requires_model,
            "neutral_value": entry.neutral_value,
            "note": entry.note,
        }
        for entry in MATRIX
    ]


def main(argv: Optional[list[str]] = None) -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--dump-matrix", action="store_true",
        help="write the mapping matrix as JSON to stdout",
    )
    parser.add_argument(
        "--classify", metavar="PARAMETER_IR",
        help="classify every parameter of a parameter IR document",
    )
    parser.add_argument(
        "--active-model", action="append", default=None, metavar="NAME",
        help="an activated SDevice model; repeatable",
    )
    parser.add_argument(
        "--allow-lossy", action="store_true",
        help="accept 'approximated' mappings instead of failing",
    )
    args = parser.parse_args(argv)

    if args.dump_matrix:
        print(json.dumps(matrix_as_json(), indent=2))
        return 0

    if not args.classify:
        parser.error("either --dump-matrix or --classify is required")

    from pathlib import Path

    ir = json.loads(Path(args.classify).read_text(encoding="utf-8"))
    results, report = classify_ir(ir, args.active_model)
    print(json.dumps(report.to_json(), indent=2))
    try:
        assert_importable(results, allow_lossy=args.allow_lossy)
    except ParameterMapError as exc:
        print(str(exc))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
