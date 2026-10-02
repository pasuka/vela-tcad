param([Parameter(Mandatory=$true)][string]$Base,[Parameter(Mandatory=$true)][int]$SupervisorPid)
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$env:Path="D:/msys64/ucrt64/bin;D:/msys64/usr/bin;$env:Path"
$root=Join-Path $Base 'postmatrix'
New-Item -ItemType Directory -Force -Path $root | Out-Null
function Save-State($value) {
    [IO.File]::WriteAllText((Join-Path $root 'status.json'),($value | ConvertTo-Json -Depth 8),(New-Object Text.UTF8Encoding($false)))
}
function Run-Checked([string]$label,[string]$program,[string]$arguments) {
    Save-State @{phase=$label;time=[DateTime]::UtcNow.ToString('o')}
    $p=Start-Process -FilePath $program -ArgumentList $arguments -WorkingDirectory "$Base/source" -WindowStyle Hidden -Wait -PassThru -RedirectStandardOutput "$root/$label.stdout.log" -RedirectStandardError "$root/$label.stderr.log"
    [IO.File]::WriteAllText("$root/$label.exit",[string]$p.ExitCode)
    if($p.ExitCode -ne 0) {throw "$label exited $($p.ExitCode); inspect retained logs"}
}
try {
    Save-State @{phase='waiting_for_complete_matrix';supervisor_pid=$SupervisorPid}
    $process=Get-Process -Id $SupervisorPid -ErrorAction SilentlyContinue
    if($process) {$process.WaitForExit()}
    $matrix=Get-Content "$Base/supervisor/summary.json" -Raw | ConvertFrom-Json
    if(-not $matrix.passed -or $matrix.points -ne 816) {throw 'Matrix not fully qualified; regression not started'}
    foreach($item in (Get-Content "$Base/postmatrix_source_hashes.json" -Raw | ConvertFrom-Json)) {
        if((Get-FileHash -LiteralPath $item.Path -Algorithm SHA256).Hash -ne $item.Hash) {throw "Regression input changed: $($item.Path)"}
    }
    Run-Checked 'configure' 'D:/msys64/ucrt64/bin/cmake.exe' '--preset windows-ucrt64-release'
    Run-Checked 'build' 'D:/msys64/ucrt64/bin/cmake.exe' '--build --preset windows-ucrt64-release --parallel 2'
    Run-Checked 'ctest' 'D:/msys64/ucrt64/bin/ctest.exe' "--preset windows-ucrt64-release --output-on-failure --parallel 2 --output-junit $root/ctest.xml"
    [xml]$test=Get-Content "$root/ctest.xml" -Raw
    if([int]$test.testsuite.tests -le 0) {throw 'CTest selected zero tests'}
    Save-State @{phase='regression_passed';tests=[int]$test.testsuite.tests;failures=[int]$test.testsuite.failures;remaining='Review and commit verified changes, then controlled performance work'}
} catch {
    Save-State @{phase='stopped';error=$_.ToString();time=[DateTime]::UtcNow.ToString('o')}
    exit 1
}
