"""Emit WP4 decks from existing JSON contracts; no command parsing or solving."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

from . import DeckEmitError, EmitOptions, emit_stage_plan, write_stage_plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ir", type=Path, required=True)
    parser.add_argument("--mesh-qualification", type=Path, required=True)
    parser.add_argument("--materials-report", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--work-dir", default="device_vela")
    parser.add_argument("--cmd-basename")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--exit-on-failure", action="store_true")
    parser.add_argument("--plot-final-only", action="store_true")
    parser.add_argument("--vtu-binary", action="store_true")
    parser.add_argument("--no-keep-intermediate", action="store_true")
    parser.add_argument("--allow-lossy", action="store_true")
    args = parser.parse_args(argv)
    def read(path: Path) -> dict:
        return json.loads(path.read_text(encoding="utf-8"))
    try:
        ir = read(args.ir)
        result = emit_stage_plan(ir, mesh_qualification=read(args.mesh_qualification),
            materials_report=read(args.materials_report), work_dir=args.work_dir,
            cmd_basename=args.cmd_basename or Path(ir["source"]["path"]).name,
            cmd_sha256=ir["source"]["sha256"], options=EmitOptions(
                exit_on_failure=args.exit_on_failure, threads=args.threads,
                plot_final_only=args.plot_final_only, vtu_binary=args.vtu_binary,
                keep_intermediate=not args.no_keep_intermediate, allow_lossy=args.allow_lossy),
            sidecar=read(args.sidecar) if args.sidecar else None)
        write_stage_plan(args.out, result)
    except DeckEmitError as error:
        print(str(error), file=sys.stderr)
        return int(error.message_id[1])
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"E406: Cannot emit stage plan: {error}", file=sys.stderr)
        return 4
    print(f"Emitted {len(result.decks)} stages to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
