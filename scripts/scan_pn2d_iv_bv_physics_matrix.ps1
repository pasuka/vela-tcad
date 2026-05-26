param(
    [string]$BaseDir = "build/pn2d_tdr_tie_probe",
    [string]$OutputSummary = ""
)

$ErrorActionPreference = "Stop"

if ($OutputSummary -eq "") {
    $OutputSummary = Join-Path $BaseDir "vela/pn2d_iv_bv_physics_matrix_summary.csv"
}

$cases = @(
    @{ name = "default"; kind = "iv"; recombination = @("srh", "auger"); bgn = "slotboom"; mobility = "default" },
    @{ name = "iv_recomb_none"; kind = "iv"; recombination = @("none"); bgn = "slotboom"; mobility = "default" },
    @{ name = "iv_bgn_none"; kind = "iv"; recombination = @("srh", "auger"); bgn = "none"; mobility = "default" },
    @{ name = "iv_srh_tau1e-6"; kind = "iv"; recombination = @("srh"); bgn = "slotboom"; mobility = "default"; taun = 1.0e-6; taup = 1.0e-6 },
    @{ name = "iv_srh_tau1e-8"; kind = "iv"; recombination = @("srh"); bgn = "slotboom"; mobility = "default"; taun = 1.0e-8; taup = 1.0e-8 },
    @{ name = "bv_recomb_none"; kind = "bv"; recombination = @("none"); bgn = "none"; mobility = "promoted_bv" },
    @{ name = "bv_recomb_srh"; kind = "bv"; recombination = @("srh"); bgn = "none"; mobility = "promoted_bv" },
    @{ name = "bv_recomb_srh_auger"; kind = "bv"; recombination = @("srh", "auger"); bgn = "none"; mobility = "promoted_bv" },
    @{ name = "bv_srh_tau1e-6"; kind = "bv"; recombination = @("srh"); bgn = "none"; mobility = "promoted_bv"; taun = 1.0e-6; taup = 1.0e-6 },
    @{ name = "bv_srh_tau1e-8"; kind = "bv"; recombination = @("srh"); bgn = "none"; mobility = "promoted_bv"; taun = 1.0e-8; taup = 1.0e-8 }
)

function Get-InterpolatedValue($rows, [string]$column, [double]$bias, [double]$scale) {
    $points = @(
        $rows |
            Where-Object { $_.bias_V -ne "" -and $_.$column -ne "" } |
            ForEach-Object {
                [pscustomobject]@{
                    bias = [double]$_.bias_V
                    value = [double]$_.$column * $scale
                }
            } |
            Sort-Object bias
    )

    if ($points.Count -eq 0) {
        return [double]::NaN
    }

    foreach ($p in $points) {
        if ([math]::Abs($p.bias - $bias) -le [math]::Max([math]::Abs($bias), 1.0) * 1.0e-12) {
            return $p.value
        }
    }

    for ($i = 0; $i -lt $points.Count - 1; $i++) {
        $p0 = $points[$i]
        $p1 = $points[$i + 1]
        if ($p0.bias -le $bias -and $bias -le $p1.bias -and $p1.bias -ne $p0.bias) {
            $t = ($bias - $p0.bias) / ($p1.bias - $p0.bias)
            return $p0.value + $t * ($p1.value - $p0.value)
        }
    }

    return [double]::NaN
}

function New-CaseConfig([hashtable]$case, [string]$ivBaseConfigPath, [string]$bvBaseConfigPath) {
    $kind = [string]$case.kind
    $name = [string]$case.name
    $baseConfigPath = if ($kind -eq "iv") { $ivBaseConfigPath } else { $bvBaseConfigPath }

    $cfg = Get-Content $baseConfigPath -Raw | ConvertFrom-Json
    $cfg.output_csv = "pn2d_${kind}_${name}.csv"
    $cfg.solver.method = "gummel_newton"
    if ($null -eq $cfg.solver.handoff) {
        $cfg.solver | Add-Member -NotePropertyName "handoff" -NotePropertyValue @{}
    }
    $cfg.solver.handoff.fallback = "none"
    $cfg.solver.handoff.require_gummel_convergence = $true
    $cfg.solver.recombination = @($case.recombination)
    $cfg.solver.bandgap_narrowing = [string]$case.bgn
    if ($case.ContainsKey("taun")) {
        if ($null -eq $cfg.solver.PSObject.Properties["taun"]) {
            $cfg.solver | Add-Member -NotePropertyName "taun" -NotePropertyValue ([double]$case.taun)
        } else {
            $cfg.solver.taun = [double]$case.taun
        }
    }
    if ($case.ContainsKey("taup")) {
        if ($null -eq $cfg.solver.PSObject.Properties["taup"]) {
            $cfg.solver | Add-Member -NotePropertyName "taup" -NotePropertyValue ([double]$case.taup)
        } else {
            $cfg.solver.taup = [double]$case.taup
        }
    }

    if ($kind -eq "bv") {
        if ($null -eq $cfg.solver.impact_ionization) {
            $cfg.solver | Add-Member -NotePropertyName "impact_ionization" -NotePropertyValue @{ model = "none" }
        } else {
            $cfg.solver.impact_ionization.model = "none"
        }
    }

    $configPath = Join-Path $BaseDir "vela/simulation_${kind}_${name}.json"
    ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $configPath
    return [pscustomobject]@{
        path = $configPath
        csv = (Join-Path $BaseDir "vela/pn2d_${kind}_${name}.csv")
    }
}

function Invoke-CaseRun([string]$configPath) {
    $status = "ok"
    $points = ""
    $err = ""

    try {
        $raw = .\build\vela_example_runner.exe --config $configPath | Out-String
        if ($raw) {
            try {
                $parsed = $raw | ConvertFrom-Json
                $points = [string]$parsed.points
            } catch {
                $status = "bad_runner_output"
                $err = "stdout_not_json"
            }
        }
    } catch {
        $status = "runner_failed"
        $err = $_.Exception.Message
    }

    return [pscustomobject]@{
        status = $status
        points = $points
        error = $err
    }
}

function Compare-Case([hashtable]$case, [string]$candidateCsv) {
    $kind = [string]$case.kind
    $name = [string]$case.name
    $referenceCsv = Join-Path $BaseDir "reference_curves/pn2d_${kind}_reference.csv"
    $reportJson = Join-Path $BaseDir "reports/pn2d_${kind}_${name}_comparison.json"
    $reportMd = Join-Path $BaseDir "reports/pn2d_${kind}_${name}_comparison.md"

    if ($kind -eq "iv") {
        $targetBias = 0.29
        $candidateScale = -1.0
        $biasMin = 0.2
        $biasMax = 0.3
    } else {
        $targetBias = 0.05
        $candidateScale = 1.0
        $biasMin = 0.05
        $biasMax = 0.05
    }

    python scripts\compare_reference_curves.py `
        --reference $referenceCsv `
        --candidate $candidateCsv `
        --output-json $reportJson `
        --output-md $reportMd `
        --kind iv `
        --candidate-column current_total_A_per_um `
        --candidate-scale $candidateScale `
        --bias-min $biasMin `
        --bias-max $biasMax | Out-Null

    $report = Get-Content $reportJson -Raw | ConvertFrom-Json
    $refRows = Import-Csv $referenceCsv
    $candRows = Import-Csv $candidateCsv

    $refAtTarget = Get-InterpolatedValue $refRows "current_total" $targetBias 1.0
    $candAtTarget = Get-InterpolatedValue $candRows "current_total_A_per_um" $targetBias $candidateScale

    $ratioAtTarget = [double]::NaN
    if ($refAtTarget -ne 0.0 -and -not [double]::IsNaN($refAtTarget) -and -not [double]::IsNaN($candAtTarget)) {
        $ratioAtTarget = [math]::Abs($candAtTarget / $refAtTarget)
    }

    return [pscustomobject]@{
        points = [int]$report.iv.points_compared
        orders = [double]$report.iv.orders_of_magnitude
        maxRelativeError = [double]$report.iv.max_relative_error
        ratioAtTarget = $ratioAtTarget
        totalAtTarget = $candAtTarget
        reportJson = $reportJson
        reportMd = $reportMd
    }
}

$ivBaseConfig = Join-Path $BaseDir "vela/simulation_iv.json"
$bvBaseConfig = Join-Path $BaseDir "vela/simulation_bv.json"
$rows = @()

foreach ($case in $cases) {
    $cfg = New-CaseConfig $case $ivBaseConfig $bvBaseConfig
    $run = Invoke-CaseRun $cfg.path

    $points = [double]::NaN
    $orders = [double]::NaN
    $maxRelativeError = [double]::NaN
    $ratioAtTarget = [double]::NaN
    $totalAtTarget = [double]::NaN

    if ($run.status -eq "ok" -and (Test-Path $cfg.csv)) {
        $cmp = Compare-Case $case $cfg.csv
        $points = $cmp.points
        $orders = $cmp.orders
        $maxRelativeError = $cmp.maxRelativeError
        $ratioAtTarget = $cmp.ratioAtTarget
        $totalAtTarget = $cmp.totalAtTarget
    }

    $rows += [pscustomobject]@{
        case = [string]$case.name
        kind = [string]$case.kind
        status = [string]$run.status
        points = $points
        orders = $orders
        max_relative_error = $maxRelativeError
        ratio_at_target = $ratioAtTarget
        total_A_per_um = $totalAtTarget
        csv_file = $cfg.csv
        config = $cfg.path
    }

    Write-Host "done $($case.name) status=$($run.status) orders=$orders"
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
