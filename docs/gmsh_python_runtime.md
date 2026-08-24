# Gmsh Python API runtime

The public offline mesh path pins Gmsh Python API 4.15.2. It does not call
Sentaurus and it never silently falls back to the legacy tensor mesh when Gmsh
is unavailable.

## Probe contract

Run the capability probe before enabling the Gmsh mesh backend:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
D:\msys64\ucrt64\bin\python.exe scripts\probe_gmsh_python.py `
  --output-dir build-release\gmsh-probe `
  --runtime-dir build-release\gmsh-python-runtime
```

On the supported MSYS2 UCRT64 installation, discovery resolves
`lib\gmsh.py`, `bin\libgmsh.dll`, and `bin\gmsh.exe`. The probe copies the
wrapper and shared library into the ignored runtime directory and names the
DLL `gmsh-4.15.dll`, as required by the packaged wrapper. It does not modify
the MSYS2 installation and the staged DLL must not be committed.

Discovery can be made explicit with `--gmsh-python`, `--gmsh-library`, and
`--gmsh-executable`, or the matching `VELA_GMSH_PYTHON`,
`VELA_GMSH_LIBRARY`, and `VELA_GMSH_EXECUTABLE` environment variables.

Success writes `gmsh_probe_manifest.json` with the wrapper, library and CLI
versions, source hashes, loaded library path, fixed topology options, and the
minimal Tri3 smoke-test counts. Failure returns exit code 2, removes any stale
success manifest, and writes strict JSON `fail_report.json`.

## CI installation

CI must install Gmsh 4.15.2 together with its Python wrapper, shared library,
and CLI. The job should then run the same probe command and retain
`gmsh_probe_manifest.json` as an artifact. On MSYS2 UCRT64 the package is:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-gmsh
```

For another platform, provide all three component paths explicitly. A wheel is
acceptable only when it is exactly 4.15.2, imports under the selected Python,
and a matching 4.15.2 `gmsh` executable is present; mixed versions fail closed.
