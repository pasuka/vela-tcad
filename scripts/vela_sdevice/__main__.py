"""Run SDevice command files with Vela via sd-sim or vela device (stdlib only)."""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import sys
from .materials import default_models_par_path
from .orchestrator import DeviceRun, version_info
from .orchestrator.support import RunError


class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.exit(2, f"E202: Option is not supported by {self.prog}: {message}\n")


# These spellings are recognized, but never silently ignored or simulated.
RESERVED_OPTIONS = {
    "-d": "Jacobian/RHS debug export",
    "-n": "Newton log suppression (use --check for input validation)",
    "-L": "Per-material parameter library export",
    "-M": "Mole-fraction-dependent regional parameter export",
    "-S": "SiC parameter export",
    "-r": "Material-library-backed parameter export",
    "-ver": "Installed simulator version selection",
    "-rel": "Installed simulator release selection",
    "-versions": "Installed simulator version listing",
    "--max_threads": "Solver/assembly thread cap (use --threads to set a count)",
    "--tcl": "Tcl command evaluation",
    "--xml": "Additional XML log",
    "--field-names": "Field names and PMI indices",
    "--parameter-names": "Rampable parameter listing",
    "--compiler-version": "Compiler version reporting",
    "-P:All": "Material-selected model parameter export (-P:<material>[:x[:y]])",
    "-L:All": "Material-selected library export (-L:<material>[:x[:y]])",
}


class ReservedOption(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        pending = list(getattr(namespace, self.dest, None) or [])
        pending.append((option_string, self.const))
        setattr(namespace, self.dest, pending)


def reserved_options(parser, argv):
    group = parser.add_argument_group("Reserved options (not implemented; exit code 2)")
    for flag, description in RESERVED_OPTIONS.items():
        group.add_argument(flag, action=ReservedOption, const=description,
                           nargs="?" if flag in ("-ver", "-rel", "--max_threads") else 0,
                           dest="reserved", metavar="VALUE", help=description)
    # Material/interface selectors are a family, not a finite list of flags.
    seen = set(RESERVED_OPTIONS)
    for token in argv:
        if token == "--":
            break
        if token.startswith(("-P:", "-L:")) and token not in seen:
            group.add_argument(token, action="append_const",
                               const=(token, "Material-selected parameter export"),
                               dest="reserved", help=argparse.SUPPRESS)
            seen.add(token)


def resolve_cmd(name):
    path = Path(name)
    for candidate in (path, Path(str(path) + ".cmd"), Path(str(path) + "_des.cmd")):
        if candidate.is_file():
            return candidate.resolve()
    raise RunError("E305", f"Command file '{name}' not found")


def main(argv=None):
    # Windows' default code page also encodes localized OS error details.
    # Keep captured CLI output in the same UTF-8 encoding as the run log.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="backslashreplace")
    argv = list(sys.argv[1:] if argv is None else argv)
    direct = bool(argv and argv[0] == "sd-sim")
    if direct:
        argv.pop(0)
    parser = Parser(prog="sd-sim" if direct else "vela", description=__doc__, allow_abbrev=False)
    if direct:
        device = parser
    else:
        commands = parser.add_subparsers(dest="command", required=True)
        device = commands.add_parser("device", help="Run a preprocessed SDevice command file", allow_abbrev=False)
    device.add_argument("cmd", nargs="?")
    device.add_argument("-v", "--version", action="store_true", dest="version")
    device.add_argument("-P", action="store_true", dest="parameters",
                        help="Export Vela's bundled models.par (deck/material selection not implemented)")
    device.add_argument("-q", action="store_true", dest="quiet")
    reserved_options(device, argv)
    device.add_argument("--check", action="store_true", dest="check")
    device.add_argument("-i", action="store_true", dest="initial_only",
                        help="Vela limitation: run only the first emitted stage; native initial-save semantics not verified")
    device.add_argument("--exit-on-failure", action="store_true")
    device.add_argument("--threads", type=int)
    device.add_argument("--verbose", action="store_true")
    device.add_argument("--allow-lossy", action="store_true")
    keep = device.add_mutually_exclusive_group()
    keep.add_argument("--keep-intermediate", action="store_true", dest="keep_intermediate", default=True)
    keep.add_argument("--no-keep-intermediate", action="store_false", dest="keep_intermediate")
    device.add_argument("--vtu-binary", action="store_true")
    device.add_argument("--plot-final-only", action="store_true")
    # parse_known_args allows reserved options with version/count operands to
    # report their actual unsupported capability before ordinary usage errors.
    options, unknown = parser.parse_known_args(argv)
    pending = list(options.reserved or [])
    if options.parameters and options.cmd:
        pending.append(("-P <commandfile>", "Deck-selected model parameter export"))
    if pending:
        for flag, description in pending:
            print(f"E202: {flag} {description} is not implemented yet; "
                  "reserved for future development. No simulation was started.", file=sys.stderr)
        return 2
    if unknown:
        parser.error("unrecognized arguments: " + " ".join(unknown))
    if options.version:
        version = version_info()
        print(f"Vela {version['version']} ({version['commit']}); SDevice syntax T-2022.03")
        return 0
    if options.parameters:
        try:
            shutil.copyfile(default_models_par_path(), Path.cwd() / "models.par")
        except OSError as error:
            print(f"E601: Cannot write output file 'models.par': {error}", file=sys.stderr)
            return 6
        return 0
    if not options.cmd:
        parser.error("No command file given (reading from stdin is not supported)")
    if options.threads is not None and options.threads < 1:
        parser.error("--threads must be positive")
    try:
        return DeviceRun(resolve_cmd(options.cmd), options).run()
    except RunError as error:
        print(str(error), file=sys.stderr)
        return error.exit_code
    except OSError as error:
        print(f"E601: Cannot read/write frontend files: {error}", file=sys.stderr)
        return 6


if __name__ == "__main__":
    raise SystemExit(main())
