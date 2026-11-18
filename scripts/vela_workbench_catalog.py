"""Materialize the numerical sd-sim CTest cases for Windows distribution."""
from pathlib import Path
import shutil

import run_bvmethods_ci as bvmethods
import run_pn2d_ci as pn2d
import run_reference_u2022_ci as u2022


ORIGINAL = "sentaurus-original.cmd"
COMPATIBLE_ARGS = ["--physical-constants-profile", "sentaurus_u2022_12",
                   "--threads", "1", "--exit-on-failure"]


def _entry(name, deck, ctest, arguments, solver, **extra):
    return {
        "name": name,
        "arguments": ["--plot-final-only", *arguments, deck],
        "work": f"{Path(deck).stem}_vela",
        "original_script": ORIGINAL,
        "linear_solver": solver,
        "ctest": ctest,
        **extra,
    }


def _stage_pn2d(catalog, reference, entries):
    pn2d.stage_reference(reference)
    campaign = pn2d.campaign
    spec = campaign.read_json(campaign.SPEC)
    names = {"0v": "PN2D Zero Bias", "iv": "PN2D Forward I–V", "bv": "PN2D Reverse BV"}
    cases = [(kind, name, None, None) for kind, name in names.items()]
    experiments = {item["id"]: item for item in spec["experiments"]}
    for variant in ("P0", "P1", "P2", "P3", "D1", "D2", "G1", "G2"):
        for branch in ("forward", "reverse"):
            cases.append((f"{variant}-{branch}", f"PN2D {variant} {branch}", experiments[variant], branch))
    for kind, name, experiment, branch in cases:
        key = f"pn2d_ci_{kind}"
        inputs = catalog / key
        inputs.mkdir()
        if experiment is None:
            source = reference / "base" / kind
            deck = f"pn2d_{kind}_sdevice.cmd"
            shorthand = {"curve": f"pn2d_{kind}.plt", "field": f"pn2d_{kind}_des.vtu"}
            shutil.copyfile(campaign.FIXTURE / "source" / deck, inputs / ORIGINAL)
        else:
            source = reference / "variants" / experiment["id"] / branch
            deck = f"{branch}_sdevice.cmd"
            direction = "正向" if branch == "forward" else "反向"
            shorthand = {"name_zh": f"PN2D {experiment['id']} {direction}"}
            original = pn2d.GRID.sub(r"\g<1>pn2d_msh.tdr\3", (source / deck).read_text(encoding="utf-8"))
            (inputs / ORIGINAL).write_text(original, encoding="utf-8")
        pn2d.stage_inputs(source, inputs, deck, spec, experiment, branch)
        solver = "umfpack" if kind == "bv" or branch == "forward" else "sparselu"
        arguments = pn2d.sdevice_arguments(deck, compensated=kind != "0v")
        entries[key] = _entry(name, deck, key, arguments[:-1], solver, **shorthand)


def _stage_u2022(catalog, entries):
    configurations = (
        ("gate_contacts", "barrier_zero_sdevice", "gate_msh", "Gate: Zero Barrier"),
        ("gate_contacts", "barrier_positive_sdevice", "gate_msh", "Gate: Positive Barrier"),
        ("gate_contacts", "workfunction_sdevice", "gate_msh", "Gate: Work Function"),
        ("schottky_iv", "schottky_n_des", "schottky_n_msh", "Schottky Forward I–V"),
        ("bjt_m1", "bjt_m1_des", "bjt_msh", "BJT M1"),
    )
    for case, stem, grid, label in configurations:
        family = u2022.REFERENCE_ROOT / u2022.CASES[case] / "ci_u2022_12"
        u2022.verify_inputs(family)
        grid_dir = u2022.REFERENCE_ROOT / "sdevice_check_only" if case == "schottky_iv" else family
        ctest = f"reference_u2022_{case}"
        key = f"{ctest}_{stem.removesuffix('_sdevice')}" if case == "gate_contacts" else ctest
        inputs = catalog / key
        u2022.stage_inputs(family, inputs, stem, [grid_dir / f"{grid}.{ext}" for ext in ("grd", "dat")])
        deck = f"{stem}.cmd"
        shutil.copyfile(family / deck, inputs / ORIGINAL)
        entries[key] = _entry(
            label, deck, ctest, COMPATIBLE_ARGS, "pardiso_mkl" if case == "bjt_m1" else "umfpack",
        )


def stage_catalog(catalog, staging):
    entries = {}
    _stage_pn2d(catalog, staging / "pn2d-reference", entries)
    _stage_u2022(catalog, entries)
    key = "bvmethods_ci_full_bv"
    work = staging / "bvmethods"
    bvmethods.stage_inputs(work)
    bvmethods.check_inputs(work / "inputs")
    inputs = catalog / key
    shutil.copytree(work / "inputs", inputs)
    shutil.copyfile(bvmethods.FAMILY / "source" / bvmethods.DECK, inputs / ORIGINAL)
    entries[key] = _entry(
        "BVMethods Full Breakdown", bvmethods.DECK, key, [*COMPATIBLE_ARGS, "--allow-lossy"],
        "pardiso_mkl", stop_criterion="current",
    )
    return entries
