# SDevice runtime packaging

The runtime distribution is built from the `SDeviceRuntime` CMake install
component. `scripts/package_vela.py` consumes offline Python archives and license
texts; it makes no network requests. Python API bindings are not needed.

## Build and package

On Windows, use the UCRT64 Release preset and build the three shipped programs:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target vela vela_example_runner sentaurus_import
python scripts/package_vela.py --build-dir build-release --python-archive build-release/wp11/downloads/python-3.12.10-embed-amd64.zip --python-sha256 4acbed6dd1c744b0376e3b1cf57ce906f9dc9e95e68824584c8099a63025a3c3 --python-license-dir build-release/wp11/licenses/python-runtime --license-root build-release/wp11/licenses --output-dir build-release/wp11/packages
```

This example names the locally retained WP11 inputs. On another machine, supply
their corresponding offline paths. The official Windows input is
[CPython 3.12.10 embeddable amd64](https://www.python.org/downloads/release/python-31210/).
The SHA-256 is checked before installation or extraction starts.

For Ubuntu, build **inside Ubuntu 24.04 x86_64**, using a separate Ninja build
directory and installed C++ dependencies (including Catch2 3 and HDF5).
The packager rejects a different distribution/version instead of labeling its
binaries as an Ubuntu 24.04 package.

```bash
cmake -S . -B build-linux24 -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-linux24 --target vela vela_example_runner sentaurus_import --parallel 4
python3 scripts/package_vela.py --build-dir build-linux24 --python-archive build-release/wp11/downloads/cpython-3.12.12-linux.tar.gz --python-sha256 6249bea270cdd4218a258edb3a1c5602dd9dc0f785b9d29963676276a3838e60 --python-license-dir build-release/wp11/licenses/python-runtime --output-dir build-release/wp11/packages
```

The retained Linux input is
[python-build-standalone release 20260114](https://github.com/astral-sh/python-build-standalone/releases/tag/20260114),
`cpython-3.12.12+20260114-x86_64-unknown-linux-gnu-install_only_stripped.tar.gz`.
Its SHA-256 matches that release's `SHA256SUMS`. The license directory contains
the unmodified `LICENSE*` texts from the corresponding upstream source release,
including `LICENSE.openssl-3.txt`; this supplies notices omitted by the minimal
runtime archives. It is used for both distributions.

When a Windows-created worktree's Git paths cannot be read inside Linux, pass
`-DVELA_SOURCE_REVISION=<verified-full-40-character-commit>` to CMake. Obtain the
commit from Git in the source host; this records it in both the launcher and the
frontend archive. Leave the option empty for ordinary Git checkouts.

## Installed layout and provenance

```text
vela/
  vela[.exe]
  bin/
    vela[.exe]
    vela_example_runner[.exe]
    sentaurus_import[.exe]
    vela_sdevice.pyz
    python/
    *.dll                       # Windows runtime dependency closure
  lib/*.so*                     # Ubuntu runtime dependency closure
  share/vela/
    materials/models.par
    schemas/ templates/ docs/ licenses/
  VERSION
  LICENSE
  THIRD_PARTY_NOTICES
  PACKAGE_MANIFEST.json
```

Both launcher entry points work after relocation. Linux executables use inherited
`$ORIGIN/../lib` RPATH for transitive library resolution. glibc and the ELF loader
remain platform dependencies; HDF5 and the other non-glibc runtime libraries are
included. Windows dependencies are collected with CMake's runtime dependency set
and checked against the actual import tables.

Windows Python's isolated `python312._pth` contains `python312.zip`, `.`, and
`..`, without `import site`. The dot is required for native stdlib extensions
such as `_socket`; all three entries resolve inside the package. The launcher
clears `PYTHONPATH`/`PYTHONHOME` when selecting bundled Python.
It also disables bytecode writes for bundled Python, keeping the package
immutable across concurrent runs. Linux staging and extracted-package smoke
use the native temporary filesystem so Windows-mounted output directories
cannot corrupt the runtime's case-sensitive paths or symlinks.

The packager generates a zip or tar.gz, a SHA-256 sidecar, and an inspectable
manifest beside the archive. The manifest records the source commit, Python
archive hash, and every shipped file's hash/size. It excludes build outputs,
reference archives and original Synopsys artifacts; the frontend's internal
zipapp is audited too.

License generation copies installed MSYS2 license texts or Debian package
copyright records for the actual dependency closure, plus compiled header
dependencies and the supplied Python notices. Missing Windows license mappings
or missing texts stop packaging. Local METIS/MUMPS builds may need
`--license-root`: the WP11 inputs include METIS's Apache notice/license and MUMPS
5.9.1's CeCILL-C, PORD, AMD and LAPACK copyright/license headers from the
[upstream MUMPS source](https://mumps-solver.org/index.php?page=dwnld).
These inputs are retained under the ignored build directory.

## Verification

Supply the produced archive to run the extracted-package test; omitting the
variable explicitly skips its smoke class.

```powershell
$env:VELA_PACKAGE_ARCHIVE = (Resolve-Path build-release/wp11/packages/vela-0.1.0-win64.zip).Path
ctest --preset windows-ucrt64-release -R '^sdevice_(package_smoke|concurrency)$' -V
```

```bash
export VELA_PACKAGE_ARCHIVE="$PWD/build-release/wp11/packages/vela-0.1.0-ubuntu24.04-x86_64.tar.gz"
ctest --test-dir build-linux24 --output-on-failure -R '^sdevice_(package_smoke|concurrency)$' -V
```

The smoke unpacks into a fresh directory containing spaces and Chinese
characters, removes system Python from PATH, runs both launcher entry points,
exports `models.par` with `-P`, and uses that file in the original PN2D IV command
fixture copied onto the synthetic four-node/two-triangle mesh. It verifies the
10 V goal and PLT/VTU/PVD publication, starts the TDR importer, and checks every
manifest hash and bundled license text. Linux additionally verifies that `ldd`
resolves every non-glibc dependency inside the relocated package.

The Linux concurrency matrix runs four distinct commands in one directory and
one command in each of four directories, comparing PLT bytes and all ASCII VTU
numeric arrays against serial runs with zero tolerance. It samples solver child
thread counts with `ps`, tests CLI/Math precedence, and verifies that a duplicate
command returns 7 without changing the first owner's lock. A controlled slow
protocol fixture is used only for the lock timing test; the numerical matrix
uses the real C++ solver. Windows runs the lock test and explicitly skips the
three Linux matrix tests as specified by WP11.

The numerical case uses coordinates in um, uniform donors at 1e15 cm^-3, and
terminal currents in A/um. It is a packaging/concurrency smoke, separate from
the physical PN2D acceptance in [WP10](validation/sdevice_wp10_closeout_2026-10-03.md).
