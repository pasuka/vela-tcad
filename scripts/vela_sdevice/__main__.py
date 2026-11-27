"""Run SDevice command files with Vela via sd-sim or vela device (stdlib only)."""
from __future__ import annotations
import argparse
import json
from decimal import InvalidOperation
from pathlib import Path
import shutil
import sys
from .log_levels import FILE_LOG_LEVELS
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
    device.add_argument("--default-par", type=Path, metavar="PAR",
                        help="Use this material defaults library beneath File{Parameter}")
    device.add_argument("--check-arithmetic-input", action="store_true",
                        help="Check DF-ISE or HDF5 TDR numeric input in the extended80 runner without solving")
    device.add_argument("--doping-mode", choices=("auto", "net-only", "reported"), default="auto",
                        help="Resolve Sentaurus doping, explicitly split net doping, or use legacy reported fields")
    device.add_argument("--datexcodes", type=Path, action="append", default=[], metavar="FILE",
                        help="Load Sentaurus species definitions (base before overrides, repeatable)")
    device.add_argument("--check-slot-ialmob", type=Path, metavar="REFERENCE_DECK",
                        help="Check the SLOT IALMob coefficients and geometry against a prepared Vela deck")
    device.add_argument("--slot-dfise-grid", type=Path, metavar="GRD",
                        help="DF-ISE export corresponding to the native SLOT grid")
    device.add_argument("--pn2d-population", action="store_true",
                        help="Use potential and carrier densities with Quasistationary controls for a classical PN2D deck")
    device.add_argument("-i", action="store_true", dest="initial_only",
                        help="Vela limitation: run only the first emitted stage; native initial-save semantics not verified")
    device.add_argument("--exit-on-failure", action="store_true")
    device.add_argument("--threads", type=int)
    device.add_argument("--verbose", action="store_true")
    device.add_argument("--log-level", choices=FILE_LOG_LEVELS, default="debug",
                        help="Minimum severity in the SDevice output log (default: debug)")
    device.add_argument("--allow-lossy", action="store_true")
    device.add_argument("--physical-constants-profile", choices=("codata2018", "sentaurus_u2022_12"),
                        help="must match the solver build; defaults to sentaurus_u2022_12; codata2018 remains an explicit rollback")
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
    if options.check_arithmetic_input:
        if options.initial_only:
            parser.error("--check-arithmetic-input cannot be combined with -i")
        options.check = True
    if bool(options.check_slot_ialmob) != bool(options.slot_dfise_grid):
        parser.error("--check-slot-ialmob and --slot-dfise-grid must be used together")
    if options.check_slot_ialmob and (options.initial_only or options.check_arithmetic_input):
        parser.error("--check-slot-ialmob is an independent input check")
    try:
        cmd = resolve_cmd(options.cmd)
        if options.check_slot_ialmob:
            from .cmd_parser import parse_cmd_file
            from .mesh_import import MeshImportError
            from .orchestrator.slot_ialmob import check_slot_ialmob
            try:
                report = check_slot_ialmob(parse_cmd_file(cmd), cmd, options.slot_dfise_grid.resolve(strict=True),
                                           options.check_slot_ialmob.resolve(strict=True))
            except (KeyError, TypeError, ValueError, InvalidOperation, MeshImportError) as error:
                raise RunError("E407", f"Invalid SLOT IALMob inputs: {error}") from error
            print(json.dumps(report, indent=2))
            return 0
        return DeviceRun(cmd, options).run()
    except RunError as error:
        print(str(error), file=sys.stderr)
        return error.exit_code
    except OSError as error:
        print(f"E601: Cannot read/write frontend files: {error}", file=sys.stderr)
        return 6


if __name__ == "__main__":
    raise SystemExit(main())
