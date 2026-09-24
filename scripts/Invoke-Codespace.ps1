<#
.SYNOPSIS
Control the Vela Codespace from PowerShell/Codex using GitHub CLI.
.DESCRIPTION
List never starts a Codespace. SSH and file-copy actions can resume a stopped
Codespace and consume quota. No action creates an instance or changes billing.
Run executes an explicitly supplied Bash command in the remote repository.
Use -DryRun to inspect an action without contacting GitHub.
#>
[CmdletBinding()]
param(
    [ValidateSet('List', 'Configure', 'Build', 'Test', 'Sync', 'Run', 'Fetch', 'Stop')]
    [string]$Action = 'List',
    [string]$Codespace,
    [string]$RemoteRoot = '/workspaces/vela-tcad',
    [ValidateRange(1, 64)][int]$Jobs = 4,
    [string]$TestRegex = '.',
    [string]$Command,
    [string]$RemotePath,
    [string]$Destination,
    [switch]$DryRun
)
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

function ConvertTo-BashLiteral([string]$Value) {
    # A single quote in Bash is represented by closing the string, an escaped
    # quote, then reopening it. No interpolation of supplied paths or regexes.
    return "'" + $Value.Replace("'", "'\''") + "'"
}

if ($Action -ne 'List' -and [string]::IsNullOrWhiteSpace($Codespace)) {
    throw 'Specify -Codespace using the name returned by -Action List.'
}
if ($Action -eq 'Run' -and [string]::IsNullOrWhiteSpace($Command)) {
    throw 'Run requires -Command (a Bash command).'
}
if ($Action -eq 'Fetch' -and ([string]::IsNullOrWhiteSpace($RemotePath) -or
                            [string]::IsNullOrWhiteSpace($Destination))) {
    throw 'Fetch requires an absolute -RemotePath and a local -Destination.'
}

$ghArgs = @('codespace')
switch ($Action) {
    'List' { $ghArgs += @('list', '--repo', 'pasuka/vela-tcad') }
    'Stop' { $ghArgs += @('stop', '-c', $Codespace) }
    'Fetch' {
        # gh 2.101 quotes literal paths for SCP, but modern Windows OpenSSH
        # treats those quotes as filename characters under SFTP. --expand
        # avoids the added quotes. Allow only shell-inert absolute paths so
        # this remains safe with either SCP or SFTP; never expand user code.
        if ($RemotePath -cnotmatch '\A/[A-Za-z0-9_./-]+\z') {
            throw 'RemotePath must be absolute and contain only ASCII letters, digits, /, _, ., or -.'
        }
        $ghArgs += @('cp', '-c', $Codespace, '--expand', '-r', "remote:$RemotePath", $Destination)
    }
    default {
        $body = switch ($Action) {
            'Configure' { "bash scripts/codespaces/build.sh configure $Jobs" }
            'Build' { "bash scripts/codespaces/build.sh build $Jobs" }
            'Test' { "bash scripts/codespaces/build.sh test $Jobs $(ConvertTo-BashLiteral $TestRegex)" }
            'Sync' {
                @'
if [[ -n "$(git status --porcelain)" ]]; then
    echo 'Remote checkout has local changes; resolve them before Sync.' >&2
    exit 1
fi
git pull --ff-only
git rev-parse HEAD
'@
            }
            'Run' { $Command }
        }
        $script = "set -euo pipefail`ncd -- $(ConvertTo-BashLiteral $RemoteRoot)`n$body`n"
        # Transport a single ASCII argument through PowerShell, gh and SSH.
        $encoded = [Convert]::ToBase64String([Text.Encoding]::UTF8.GetBytes($script))
        $ghArgs += @('ssh', '-c', $Codespace, '--', '-T',
            "printf %s $encoded | base64 --decode | bash")
        if ($DryRun) { Write-Output $script }
    }
}
if ($DryRun) {
    Write-Output ($ghArgs | ConvertTo-Json -Compress)
    return
}

# Reuse the existing installations even in a newly opened PowerShell window.
# Restore PATH on exit; do not modify the user's or system's registry settings.
$savedPath = $env:Path
try {
    foreach ($toolDir in @('D:\msys64\ucrt64\bin', 'C:\Program Files\GitHub CLI')) {
        if (Test-Path -LiteralPath $toolDir) { $env:Path = "$toolDir;$env:Path" }
    }
    $gh = Get-Command gh -CommandType Application -ErrorAction Stop
    & $gh.Source @ghArgs
    if ($LASTEXITCODE -ne 0) { throw "GitHub CLI failed with exit code $LASTEXITCODE." }
}
finally { $env:Path = $savedPath }
