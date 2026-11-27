"""Generate independent SLOT Newton stabilization controls without changing gates."""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path


def configure_stabilization(config: dict, *, regularization: bool,
                            closure_replacement: bool,
                            regularization_scale: float = 1e-8) -> dict:
    if not math.isfinite(regularization_scale) or regularization_scale <= 0:
        raise ValueError("enabled regularization scale must be finite and positive")
    result = copy.deepcopy(config)
    solver = result.setdefault("solver", {})
    solver["carrier_regularization_scale"] = regularization_scale if regularization else 0.0
    closure = solver.setdefault("global_continuity_closure", {})
    closure["linear_system_stabilization"] = (
        "replace_free_row_sum" if closure_replacement else "off")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--carrier-regularization", choices=("on", "off"), required=True)
    parser.add_argument("--closure-replacement", choices=("on", "off"), required=True)
    parser.add_argument("--regularization-scale", type=float, default=1e-8)
    args = parser.parse_args()
    config = configure_stabilization(
        json.loads(args.input.read_text(encoding="utf-8")),
        regularization=args.carrier_regularization == "on",
        closure_replacement=args.closure_replacement == "on",
        regularization_scale=args.regularization_scale)
    # Existing evidence and input files must not be overwritten.
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(config, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
