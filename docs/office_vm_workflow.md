# Remote compute VM workflow

Use `scripts/Invoke-VelaVm.ps1` for SSH file transfer, isolated Vela/Sentaurus
environments and detached tasks. All examples below use placeholders. Real
addresses, ports, usernames, installation paths and validation job records stay
outside the repository.

## Private configuration

Keep connection details in the user's `~/.ssh/config` (on Windows,
`$env:USERPROFILE\.ssh\config`). Replace the example fields locally:

```sshconfig
Host vela-vm
    HostName HOST_ADDRESS
    Port SSH_PORT
    User REMOTE_USER
    IdentityFile ~/.ssh/PRIVATE_KEY_FILENAME
    IdentitiesOnly yes
    BatchMode yes
    StrictHostKeyChecking yes
    ConnectTimeout 15
    ServerAliveInterval 30
    ServerAliveCountMax 3
```

Store machine-specific paths in
`$env:USERPROFILE\.config\vela-vm\config.json`, outside the repository. Start
with [the example](../scripts/remote-vm/config.example.json), replacing its
placeholders. Its fields are `SshHost`, `RemoteBase` (an absolute writable remote
workspace) and `SentaurusRoot` (the directory containing the vendor `bin`).
Use `-ConfigFile` for an alternate private profile. Explicit `-RemoteBase` and
`-SentaurusRoot` override those paths. An explicit SSH alias must match the loaded
profile to avoid applying another machine's paths.

Do not copy the private JSON, SSH files, raw remote logs or connection transcripts
into public documentation. Keep detailed job IDs/source paths in user-local notes
such as `~/.config/vela-vm/validation-*.md`; build logs may also stay in the ignored
`build/` directory. Ignored files must not be force-added to Git.

## Environment behavior

The SSH command transport works with csh or Bash login shells. It clears the
inherited vendor library path before starting Bash, and skips shell startup files.
It does not alter the remote login shell, system profiles or vendor installation.

- **Vela**: `~/.local/bin/micromamba`, root `~/micromamba`, environment
  `vela-gcc16`. `Probe` checks the actual compiler and dependency versions.
  Installing MKL does not itself implement a PARDISO backend in Vela.
- **Sentaurus**: uses privately configured `SentaurusRoot`, inherited license
  variables and optional `~/.config/vela-vm/sentaurus.env` (Bash syntax). No
  micromamba activation. Help output alone does not validate license checkout.
- **Plain**: basic system paths with conflicting library/Python variables cleared.

## Upload and build

Run from the current repository checkout in PowerShell:

```powershell
./scripts/Invoke-VelaVm.ps1 -Action Probe
./scripts/Invoke-VelaVm.ps1 -Action Probe -Environment Sentaurus
./scripts/Invoke-VelaVm.ps1 -Action Sync
```

`Sync` returns `REMOTE_DIR=<private-workspace>/sources/<unique-id>`. It includes
current contents of tracked files and non-ignored untracked files, omits deleted
files, ignored outputs and `.git`, and rejects symlinks/submodules requiring manual
review. Review `git status` first. Every snapshot uses a new directory and cannot
overwrite an active job's source. Snapshots are not Git checkouts.

```powershell
$remote = '<REMOTE_DIR returned by Sync>'
./scripts/Invoke-VelaVm.ps1 -Action Configure -RemoteDir $remote
./scripts/Invoke-VelaVm.ps1 -Action Build -RemoteDir $remote -Jobs 8
./scripts/Invoke-VelaVm.ps1 -Action Test -RemoteDir $remote -Jobs 4
```

Stop after a failed step. Configuration verifies UMFPACK, SPQR, METIS and
HDF5/HighFive state support. Build directory: `build-mamba-gcc16-eigen5-release`.
GCC 16 is required. Check actual Eigen discovery when changing dependencies.
Choose build/test concurrency according to the host allocation; example job
counts do not imply exclusive ownership of a machine's CPUs.

## Detached tasks and file transfer

`Run` is synchronous and propagates failures. `Start` detaches a task with
`nohup`/`setsid` and returns `JOB_ID` and `JOB_DIR`, storing output and exit status.
It survives ordinary SSH disconnection, not VM reboot. Avoid `&` inside
`-Command`: the wrapper must wait for the real simulation process.

```powershell
./scripts/Invoke-VelaVm.ps1 -Action Start -RemoteDir $remote `
  -Command 'bash scripts/remote-vm/build.sh build 8'
./scripts/Invoke-VelaVm.ps1 -Action Status -JobId '<returned-job-id>'

# Use an existing writable parent and a fresh case directory.
$cases = '<absolute-remote-cases-directory>'
./scripts/Invoke-VelaVm.ps1 -Action Upload -LocalPath ./build/my-case `
  -RemotePath $cases
./scripts/Invoke-VelaVm.ps1 -Action Start -Environment Sentaurus `
  -RemoteDir "$cases/my-case" -Command 'sdevice device.cmd'
./scripts/Invoke-VelaVm.ps1 -Action Status -JobId '<returned-job-id>'
./scripts/Invoke-VelaVm.ps1 -Action Fetch -RemotePath "$cases/my-case" `
  -LocalPath ./build/vm-results
```

Use a standalone prepared deck with all mesh/parameter inputs present.
Workbench templates with `@...@` placeholders require preprocessing. Check
license availability and actual scientific outputs before claiming success.
`Status` reports FINISHED plus EXIT_CODE, RUNNING, or UNKNOWN_OR_INTERRUPTED;
the last state is not success.

Remote paths must be absolute ASCII paths without spaces or `..`. Upload/Fetch
use recursive SCP and may overwrite files at the selected destination. Source
snapshots and job outputs remain until explicitly cleaned up. `-DryRun` prints
resolved commands and can expose local/private paths; do not publish its output.

## Future chats and verification

Use the configured SSH alias and private profile; do not print their contents in
public artifacts. Do not substitute another host's paths or license settings.
Reuse installed tools, upload a fresh snapshot as needed, and record source path,
job ID, environment, input deck and results in private notes. Only submit requested
jobs, checking exit codes and scientific output.

Verify changes with PowerShell/Bash syntax checks, the snapshot regression
`tests/regression/test_remote_vm_snapshot.py`, remote environment probes, an
upload/download round-trip, and detached success/failure jobs. Build changes also
need project compilation and relevant CTest coverage. Store identifying evidence
outside tracked files.
