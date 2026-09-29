"""Command-line entry: ``python -m vela_sdevice.materials`` (``PYTHONPATH=scripts``).

Writes ``materials.json`` (C3) and ``materials_report.json`` into ``--out``.
Exit codes follow contract C8: 0 success, 3 missing parameter file, 4 validation
failure (E403/E404/E410), 2 usage error.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from vela_sdevice.materials import (MaterialsError, active_models_from_ir, default_models_par_path,
                                    emit_materials, write_materials)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m vela_sdevice.materials", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--par", type=Path, required=True, help="Sentaurus parameter file (File{Parameter})")
    parser.add_argument("--ir", type=Path, help="execution IR v2 (C1); supplies active models and temperature")
    parser.add_argument("--model", action="append", default=[], help="active matrix model (repeatable; overrides --ir)")
    parser.add_argument("--temperature", type=float, help="lattice temperature in K (overrides --ir)")
    parser.add_argument("--materials", nargs="*", default=[], help="Vela material names present in the mesh")
    parser.add_argument("--mesh-qualification", type=Path, help="C2 report; its `materials` list is used")
    parser.add_argument("--electrode-material", action="append", default=[],
                        help="Electrode{Material=} name whose Bandgap.WorkFunction is needed (repeatable)")
    parser.add_argument("--allow-lossy", action="store_true", help="accept `approximated` parameters")
    parser.add_argument("--band-edge-sensitive", action="store_true",
                        help="absolute band edges are observable (gate/Schottky contacts, heterointerfaces)")
    parser.add_argument("--default-par", type=Path, default=default_models_par_path(),
                        help="Vela default library (share/vela/materials/models.par)")
    parser.add_argument("--benchmark", help="C3 benchmark string (default sdevice:<par stem>)")
    parser.add_argument("--out", type=Path, required=True, help="output directory")
    args = parser.parse_args(argv)

    models = list(args.model)
    temperature = args.temperature
    effective = None
    if args.ir is not None:
        ir = json.loads(args.ir.read_text(encoding="utf-8"))
        effective = ir.get("physics", {}).get("effective", [])
        if not models:
            models = active_models_from_ir(ir)
        if temperature is None:
            temperature = float(ir.get("physics", {}).get("temperature_K", 300.0))
    if temperature is None:
        temperature = 300.0
    materials = list(args.materials)
    if args.mesh_qualification is not None:
        report = json.loads(args.mesh_qualification.read_text(encoding="utf-8"))
        materials.extend(report.get("materials", []))
    if not materials:
        parser.error("no mesh materials given (--materials or --mesh-qualification)")

    try:
        result = emit_materials(args.par, active_models=models, temperature_K=temperature, mesh_materials=materials,
                                electrode_materials=args.electrode_material, allow_lossy=args.allow_lossy,
                                band_edge_sensitive=args.band_edge_sensitive, default_par=args.default_par,
                                benchmark=args.benchmark, effective_models=effective)
    except MaterialsError as problem:
        print(f"{problem.message_id}: {problem.message}", file=sys.stderr)
        return 3 if problem.message_id.startswith("E3") else 4
    materials_path, report_path = write_materials(args.out, result)
    print(json.dumps({"materials": materials_path.as_posix(), "report": report_path.as_posix(),
                      "counts": result.report["counts"], "lossy": len(result.report["lossy"]),
                      "inert": len(result.report["inert"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
