"""``python -m vela_sdevice.cmd_parser <deck.cmd>`` prints the execution IR v2 and exits with verdict.exit_code."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import parse_cmd_file


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="vela_sdevice.cmd_parser",
                                     description="Parse an SDevice command file into the Vela execution IR v2.")
    parser.add_argument("deck", type=Path)
    parser.add_argument("--output", type=Path, help="write the IR JSON here instead of stdout")
    args = parser.parse_args(argv)
    if not args.deck.is_file():
        print(f"E305 Command file '{args.deck}' not found", file=sys.stderr)
        return 3
    ir = parse_cmd_file(args.deck)
    payload = json.dumps(ir, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        sys.stdout.write(payload)
    return int(ir["verdict"]["exit_code"])


if __name__ == "__main__":
    raise SystemExit(main())
