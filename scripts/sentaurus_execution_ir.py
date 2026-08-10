#!/usr/bin/env python3
"""SDevice ``.cmd`` summary -> Vela execution IR (``vela.sentaurus_execution_ir.v1``).

The execution IR is the stable contract between the Sentaurus frontend and any
Vela orchestration. It records:

``electrodes``      terminal names and their initial bias
``initial_solve``   the ordered pre-sweep solve stages
``stages``          the ordered ``Quasistationary`` sweeps
``physics``         models classified against an explicit support whitelist
``analysis``        the analysis kind derived from the IR, never supplied
                    externally
``unsupported``     fail-closed report for models Vela cannot honour

Model classification inverts the historical blacklist: a model is accepted only
when it appears in :data:`SUPPORTED_MODELS` or :data:`METADATA_ONLY_MODELS`.
Everything else is reported as unsupported, and :func:`build_execution_ir`
refuses to produce an IR for it unless the caller explicitly opts out.
"""

from __future__ import annotations

from typing import Any, Iterable


EXECUTION_IR_SCHEMA = "vela.sentaurus_execution_ir.v1"

# Models that map onto an implemented Vela physics path.
SUPPORTED_MODELS: dict[str, str] = {
    "Mobility": "mobility model container",
    "DopingDependence": "doping-dependent low-field mobility",
    "DopingDep": "doping-dependent low-field mobility",
    "HighFieldSaturation": "high-field mobility saturation",
    "HighFieldsaturation": "high-field mobility saturation",
    "Eparallel": "high-field saturation driving force (parallel field)",
    "GradQuasiFermi": "high-field saturation driving force (quasi-Fermi gradient)",
    "E2": "high-field saturation driving force variant",
    "ElectricField": "electric-field driving force",
    "Recombination": "recombination model container",
    "SRH": "Shockley-Read-Hall recombination",
    "Auger": "Auger recombination",
    "Avalanche": "impact-ionization generation",
    "VanOverstraeten": "van Overstraeten-de Man ionization coefficients",
    "Band2Band": "band-to-band tunneling generation",
    "EffectiveIntrinsicDensity": "effective intrinsic density container",
    "OldSlotboom": "Slotboom bandgap narrowing",
    "Fermi": "Fermi-Dirac carrier statistics",
}

# Models recorded for traceability that provably do not change the Vela solve.
METADATA_ONLY_MODELS: dict[str, str] = {
    "Trap": "trap kinetics imported as metadata only; no trap occupancy is solved",
    "Traps": "trap kinetics imported as metadata only; no trap occupancy is solved",
}

# Models known to change the solved equations in ways Vela cannot reproduce.
KNOWN_UNSUPPORTED_MODELS: dict[str, str] = {
    "Thermodynamic": "Vela does not solve the lattice temperature equation",
    "Hydrodynamic": "Vela does not solve carrier energy balance equations",
    "eTemperature": "carrier temperature transport is not implemented",
    "hTemperature": "carrier temperature transport is not implemented",
    "IALMob": "IALMob surface-orientation mobility is not implemented",
    "Enormal": "normal-field surface mobility degradation is not calibrated "
               "against the Sentaurus Enormal model",
}

ANALYSIS_KINDS = ("equilibrium", "iv", "bv", "cv")


class ExecutionIrError(ValueError):
    """Raised when an SDevice summary cannot be turned into an execution IR."""


def classify_models(models: Iterable[str]) -> dict[str, list[dict[str, str]]]:
    """Split model tokens into supported / metadata-only / unsupported."""
    supported: list[dict[str, str]] = []
    metadata_only: list[dict[str, str]] = []
    unsupported: list[dict[str, str]] = []
    for model in sorted(dict.fromkeys(models)):
        if model in SUPPORTED_MODELS:
            supported.append({"model": model, "mapping": SUPPORTED_MODELS[model]})
        elif model in METADATA_ONLY_MODELS:
            metadata_only.append({"model": model, "reason": METADATA_ONLY_MODELS[model]})
        elif model in KNOWN_UNSUPPORTED_MODELS:
            unsupported.append({"model": model, "reason": KNOWN_UNSUPPORTED_MODELS[model]})
        else:
            unsupported.append({
                "model": model,
                "reason": "model is not in the Vela SDevice support whitelist",
            })
    return {
        "supported": supported,
        "metadata_only": metadata_only,
        "unsupported": unsupported,
    }


def _numeric(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _has_small_signal(cmd_summary: dict[str, Any]) -> bool:
    flags = {str(flag).lower() for flag in cmd_summary.get("math", {}).get("flags", [])}
    parameters = {str(key).lower() for key in cmd_summary.get("math", {}).get("parameters", {})}
    for stage in cmd_summary.get("solve", {}).get("initial_steps", []):
        flags.update(str(item).lower() for item in stage.get("equations", []))
    for sweep in cmd_summary.get("sweeps", []):
        flags.update(str(item).lower() for item in sweep.get("equations", []))
    tokens = flags | parameters
    return bool({"accoupled", "acmethod", "frequency"} & tokens)


def derive_analysis(cmd_summary: dict[str, Any],
                    classification: dict[str, list[dict[str, str]]]) -> dict[str, Any]:
    """Derive the analysis kind from the parsed SDevice content.

    The rules are deliberately explicit so that no caller has to supply an
    external ``kind``. An input that satisfies two mutually exclusive rules is
    rejected rather than resolved by precedence.
    """
    sweeps = cmd_summary.get("sweeps", [])
    supported_models = {item["model"] for item in classification["supported"]}
    has_avalanche = "Avalanche" in supported_models
    small_signal = _has_small_signal(cmd_summary)

    if not sweeps:
        if small_signal:
            raise ExecutionIrError(
                "small-signal analysis requested without any Quasistationary "
                "sweep; cannot derive an analysis kind")
        return {
            "kind": "equilibrium",
            "reason": "no Quasistationary sweep is present",
            "sweep_contacts": [],
        }

    if small_signal and has_avalanche:
        raise ExecutionIrError(
            "ambiguous analysis: the deck requests both small-signal (CV) "
            "solve steps and avalanche generation; split the SDevice command "
            "file into separate analyses")

    contacts = [str(sweep.get("contact")) for sweep in sweeps]
    stops = [_numeric(sweep.get("stop")) for sweep in sweeps]
    if any(stop is None for stop in stops):
        unresolved = [
            sweep.get("stop") for sweep, stop in zip(sweeps, stops) if stop is None
        ]
        raise ExecutionIrError(
            "Quasistationary goal voltage is not numeric "
            f"({unresolved!r}); expand Sentaurus Workbench placeholders with "
            "--template-var before deriving the analysis")

    if small_signal:
        kind = "cv"
        reason = "a small-signal (AC) solve step drives the sweep"
    elif has_avalanche:
        kind = "bv"
        reason = "avalanche generation is enabled during a voltage sweep"
    else:
        kind = "iv"
        reason = "voltage sweep without avalanche or small-signal analysis"

    return {
        "kind": kind,
        "reason": reason,
        "sweep_contacts": contacts,
        "final_contact": contacts[-1],
        "final_voltage": stops[-1],
    }


def _stage_entries(cmd_summary: dict[str, Any]) -> list[dict[str, Any]]:
    stages: list[dict[str, Any]] = []
    for index, step in enumerate(cmd_summary.get("solve", {}).get("initial_steps", [])):
        stages.append({
            "index": index,
            "phase": "initial",
            "type": step.get("type", "Coupled"),
            "equations": list(step.get("equations", [])),
            "parameters": dict(step.get("parameters", {})),
            "depends_on": [index - 1] if index else [],
        })
    offset = len(stages)
    for index, sweep in enumerate(cmd_summary.get("sweeps", [])):
        stage_index = offset + index
        stages.append({
            "index": stage_index,
            "phase": "sweep",
            "type": "Quasistationary",
            "contact": sweep.get("contact"),
            "goal_voltage": _numeric(sweep.get("stop")),
            "equations": list(sweep.get("equations", [])),
            "step_control": dict(sweep.get("step_control", {})),
            "depends_on": [stage_index - 1] if stage_index else [],
        })
    return stages


def build_execution_ir(cmd_summary: dict[str, Any],
                       source: str,
                       models: Iterable[str],
                       allow_unsupported: bool = False) -> dict[str, Any]:
    """Build ``vela.sentaurus_execution_ir.v1`` from a parsed SDevice summary.

    Raises :class:`ExecutionIrError` when the deck uses models outside the
    support whitelist, unless ``allow_unsupported`` is set by a caller that only
    needs the report.
    """
    unresolved = cmd_summary.get("unresolved_placeholders", [])
    if unresolved:
        raise ExecutionIrError(
            f"{source}: unresolved Sentaurus Workbench placeholders "
            f"{', '.join(unresolved)}; supply --template-var for each")

    classification = classify_models(models)
    if classification["unsupported"] and not allow_unsupported:
        details = "; ".join(
            f"{item['model']} ({item['reason']})"
            for item in classification["unsupported"]
        )
        raise ExecutionIrError(
            f"{source}: unsupported SDevice physics: {details}")

    electrodes = [
        {
            "name": str(item.get("name")),
            "voltage": _numeric(item.get("voltage")),
        }
        for item in cmd_summary.get("electrodes", [])
    ]
    if not electrodes:
        raise ExecutionIrError(f"{source}: no Electrode block was found")
    missing_bias = [item["name"] for item in electrodes if item["voltage"] is None]
    if missing_bias:
        raise ExecutionIrError(
            f"{source}: electrodes without a numeric initial voltage: "
            f"{', '.join(missing_bias)}")

    analysis = derive_analysis(cmd_summary, classification)
    stages = _stage_entries(cmd_summary)
    sweep_contacts = {stage.get("contact") for stage in stages if stage["phase"] == "sweep"}
    electrode_names = {item["name"] for item in electrodes}
    unknown_contacts = sorted(
        str(name) for name in sweep_contacts if name not in electrode_names)
    if unknown_contacts:
        raise ExecutionIrError(
            f"{source}: Quasistationary goal references undeclared electrode(s): "
            f"{', '.join(unknown_contacts)}")

    return {
        "schema": EXECUTION_IR_SCHEMA,
        "source": source,
        "files": dict(cmd_summary.get("files", {})),
        "electrodes": electrodes,
        "thermodes": list(cmd_summary.get("thermodes", [])),
        "initial_solve": [stage for stage in stages if stage["phase"] == "initial"],
        "stages": stages,
        "analysis": analysis,
        "physics": {
            "supported": classification["supported"],
            "metadata_only": classification["metadata_only"],
        },
        "unsupported": classification["unsupported"],
        "math": dict(cmd_summary.get("math", {})),
    }
