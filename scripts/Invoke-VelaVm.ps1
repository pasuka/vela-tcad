<#
.SYNOPSIS
SSH automation for a privately configured Linux compute host.
.DESCRIPTION
Sync uploads a new snapshot of tracked and non-ignored untracked working files.
It never overwrites an existing snapshot. Run is synchronous; Start detaches a
job and returns its ID. Status reads its log and recorded exit code. Sentaurus
commands must reference prepared standalone decks; Workbench preprocessing is
not implied. No action changes system profiles, installs packages, or kills jobs.
#>
[CmdletBinding()]
param(
    [ValidateSet('Probe','Sync','Upload','Fetch','Configure','Build','Test','Run','Start','Status')]
    [string]$Action = 'Probe',
    [string]$SshHost,
    [string]$RemoteBase,
    [string]$SentaurusRoot,
    [string]$ConfigFile = (Join-Path $env:USERPROFILE '.config/vela-vm/config.json'),
    [ValidateSet('Vela','Sentaurus','Plain')][string]$Environment = 'Vela',
    [string]$RemoteDir,
    [string]$LocalPath,
    [string]$RemotePath,
    [string]$Command,
    [string]$JobId,
    [ValidateRange(1,32)][int]$Jobs = 8,
    [string]$TestRegex = '.',
    [switch]$DryRun
)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
$repoRoot = Split-Path $PSScriptRoot -Parent
# Host identity stays in ~/.ssh/config; machine paths stay outside the repository.
if (Test-Path -LiteralPath $ConfigFile) {
    $settings = Get-Content -LiteralPath $ConfigFile -Raw | ConvertFrom-Json
    if ($SshHost -and $SshHost -ne $settings.SshHost) {
        throw 'SSH alias differs from the private configuration. Supply a matching -ConfigFile.'
    }
    if (!$SshHost) { $SshHost = [string]$settings.SshHost }
    if (!$RemoteBase) { $RemoteBase = [string]$settings.RemoteBase }
    if (!$SentaurusRoot) { $SentaurusRoot = [string]$settings.SentaurusRoot }
}
if (!$SshHost) { $SshHost = 'vela-vm' }

function Quote-Bash([string]$Value) { return "'" + $Value.Replace("'", "'\''") + "'" }
function Invoke-Native([string]$Program, [string[]]$Arguments) {
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Program failed (exit $LASTEXITCODE)" }
}
function Invoke-Remote([string]$Script) {
    if ($DryRun) { Write-Output $Script; return }
    $encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($Script.Replace("`r`n", "`n")))
    # This outer command is valid in csh as well as Bash. Clear the vendor library
    # path before starting Bash, and never source the polluted login profiles.
    $remote = "printf %s $encoded | /usr/bin/base64 -d | /usr/bin/env -u LD_LIBRARY_PATH -u LD_PRELOAD -u BASH_ENV /bin/bash --noprofile --norc"
    Invoke-Native 'ssh' @('-T','-o','BatchMode=yes',$SshHost,$remote)
}
function Check-RemotePath([string]$Path) {
    if ($Path -cnotmatch '\A/[A-Za-z0-9_./-]+\z' -or $Path -match '(^|/)\.\.(/|$)') {
        throw 'Use an absolute remote path with letters, digits, /, _, -, . and no .. components.'
    }
}
function Environment-Script([string]$Choice) {
    $common = @'
set -euo pipefail
export PATH="$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin"
unset LD_LIBRARY_PATH LD_PRELOAD PYTHONPATH PYTHONHOME
'@
    if ($Choice -eq 'Vela') {
        return $common + "`n" + @'
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
eval "$("$HOME/.local/bin/micromamba" shell hook --shell bash)"
micromamba activate vela-gcc16
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
'@
    }
    if ($Choice -eq 'Sentaurus') {
        if (!$SentaurusRoot) { throw 'Set SentaurusRoot in the private config or pass -SentaurusRoot.' }
        Check-RemotePath $SentaurusRoot
        return $common + "`nexport STROOT=$(Quote-Bash $SentaurusRoot)`n" + @'
export PATH="$STROOT/bin:$PATH"
# Optional user-only license/environment settings; never load Vela here.
if [[ -f "$HOME/.config/vela-vm/sentaurus.env" ]]; then
    source "$HOME/.config/vela-vm/sentaurus.env"
fi
'@
    }
    return $common
}

if ($SshHost -notmatch '\A[A-Za-z0-9][A-Za-z0-9_.-]*\z') { throw 'Invalid SSH alias' }
if ($Action -in @('Sync','Start','Status')) {
    if (!$RemoteBase) { throw 'Set RemoteBase in the private config or pass -RemoteBase.' }
    Check-RemotePath $RemoteBase
}
if ($Action -in @('Upload','Fetch')) {
    Check-RemotePath $RemotePath
    if (!$LocalPath) { throw 'Specify -LocalPath' }
    if ($Action -eq 'Upload') {
        $resolved = (Resolve-Path -LiteralPath $LocalPath).Path
        $copyArgs = @('-r', $resolved, "${SshHost}:$RemotePath")
    } else {
        $copyArgs = @('-r', "${SshHost}:$RemotePath", [IO.Path]::GetFullPath($LocalPath))
    }
    if ($DryRun) { $copyArgs; return }
    if ($Action -eq 'Fetch') {
        $parent = Split-Path ([IO.Path]::GetFullPath($LocalPath)) -Parent
        [IO.Directory]::CreateDirectory($parent) | Out-Null
    }
    Invoke-Native 'scp' (@('-o','BatchMode=yes') + $copyArgs)
    return
}
if ($Action -eq 'Sync') {
    $tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
    $target = "$remoteBase/sources/$tag"
    if ($DryRun) { "Upload Git-visible working files from $repoRoot to $target"; return }
    $archive = Join-Path ([IO.Path]::GetTempPath()) "vela-vm-$tag.tar.gz"
    try {
        $python = 'D:\msys64\ucrt64\bin\python.exe'
        if (!(Test-Path -LiteralPath $python)) { throw "Missing $python" }
        Invoke-Native $python @("$PSScriptRoot/remote-vm/snapshot.py",$repoRoot,$archive)
        Invoke-Remote "set -eu`nmkdir -p '$remoteBase/incoming' '$remoteBase/sources'`nmkdir '$target'"
        Invoke-Native 'scp' @('-o','BatchMode=yes',$archive,"${SshHost}:$remoteBase/incoming/$tag.tar.gz")
        Invoke-Remote "set -eu`ntar -xzf '$remoteBase/incoming/$tag.tar.gz' -C '$target'`nrm -- '$remoteBase/incoming/$tag.tar.gz'`nprintf 'REMOTE_DIR=%s\n' '$target'"
    } finally {
        if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive }
    }
    return
}
if ($Action -eq 'Status') {
    if ($JobId -cnotmatch '\A[0-9]{8}-[0-9]{6}-[a-f0-9]{8}\z') { throw 'Specify a returned -JobId' }
    Invoke-Remote (@'
set -eu
job=JOB_PATH
test -d "$job"
if [[ -f "$job/exit-code" ]]; then
    printf 'STATE=FINISHED EXIT_CODE=%s\n' "$(cat "$job/exit-code")"
elif [[ -f "$job/pid" ]] && kill -0 "$(cat "$job/pid")" 2>/dev/null \
    && grep -Fzq "$job/runner.sh" "/proc/$(cat "$job/pid")/cmdline"; then
    echo STATE=RUNNING
else
    echo STATE=UNKNOWN_OR_INTERRUPTED
fi
tail -n 60 "$job/output.log"
'@).Replace('JOB_PATH',(Quote-Bash "$remoteBase/jobs/$JobId"))
    return
}
$body = switch ($Action) {
    'Probe' {
        if ($Environment -eq 'Vela') { @'
hostname
printf 'CONDA_PREFIX=%s\n' "$CONDA_PREFIX"
"$CXX" --version
cmake --version
python -c 'import numpy,h5py,PIL,jsonschema; print("Python dependencies OK")'
micromamba list 'eigen|mkl|highfive|suitesparse'
'@ } elseif ($Environment -eq 'Sentaurus') { 'command -v sdevice sde snmesh; sdevice -h' }
        else { 'hostname; id; pwd' }
    }
    'Configure' { "bash scripts/remote-vm/build.sh configure $Jobs" }
    'Build' { "bash scripts/remote-vm/build.sh build $Jobs" }
    'Test' { "bash scripts/remote-vm/build.sh test $Jobs $(Quote-Bash $TestRegex)" }
    default { if (!$Command) { throw 'Specify -Command' }; $Command }
}
if ($Action -in @('Configure','Build','Test') -and $Environment -ne 'Vela') { throw 'Build actions require Vela' }
if ($Action -ne 'Probe' -and !$RemoteDir) { throw 'Specify -RemoteDir (Sync returns a fresh source directory)' }
$script = (Environment-Script $Environment) + "`n"
if ($RemoteDir) { Check-RemotePath $RemoteDir; $script += "cd -- $(Quote-Bash $RemoteDir)`n" }
$script += "$body`n"
if ($Action -ne 'Start') { Invoke-Remote $script; return }

$tag = (Get-Date -Format 'yyyyMMdd-HHmmss') + '-' + [guid]::NewGuid().ToString('N').Substring(0,8)
$job = "$remoteBase/jobs/$tag"
$payload = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($script.Replace("`r`n", "`n")))
$runner = @'
#!/usr/bin/env bash
job=JOB_PATH
trap 'code=$?; printf "%s\n" "$code" > "$job/exit-code.tmp"; mv "$job/exit-code.tmp" "$job/exit-code"' EXIT
date -u +%FT%TZ > "$job/started-at"
/bin/bash --noprofile --norc "$job/task.sh"
exit $?
'@
$runner = $runner.Replace('JOB_PATH',(Quote-Bash $job))
$runner64 = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($runner.Replace("`r`n", "`n")))
Invoke-Remote @"
set -euo pipefail
mkdir -p '$remoteBase/jobs'
mkdir '$job'
printf %s $payload | base64 -d > '$job/task.sh'
printf %s $runner64 | base64 -d > '$job/runner.sh'
nohup setsid /bin/bash --noprofile --norc '$job/runner.sh' > '$job/output.log' 2>&1 < /dev/null &
printf '%s\n' "`$!" > '$job/pid'
printf 'JOB_ID=%s\nJOB_DIR=%s\n' '$tag' '$job'
"@
