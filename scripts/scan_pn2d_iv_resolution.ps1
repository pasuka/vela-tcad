param(
    [string]$BaseDir = "build/pn2d_tdr_tie_probe",
    [string]$OutputSummary = "",
    [double]$BiasMin = 0.2,
    [double]$BiasMax = 0.3
)

$ErrorActionPreference = "Stop"

if ($OutputSummary -eq "") {
    $OutputSummary = Join-Path $BaseDir "vela/pn2d_iv_resolution_summary.csv"
}

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

function New-ResolutionConfig([double]$step, [string]$baseConfigPath, [string]$tag) {
    $cfg = Get-Content $baseConfigPath -Raw | ConvertFrom-Json
    $cfg.output_csv = "pn2d_iv_resolution_${tag}.csv"

    if ($null -eq $cfg.solver) {
        $cfg | Add-Member -NotePropertyName "solver" -NotePropertyValue @{}
    }
    $cfg.solver.method = "gummel_newton"

    if ($null -eq $cfg.solver.handoff) {
        $cfg.solver | Add-Member -NotePropertyName "handoff" -NotePropertyValue @{}
    }
    $cfg.solver.handoff.fallback = "none"
    $cfg.solver.handoff.require_gummel_convergence = $true

    if ($null -eq $cfg.sweep) {
        throw "simulation_iv.json is missing sweep object"
    }
    $cfg.sweep.step = $step

    $configPath = Join-Path $BaseDir "vela/simulation_iv_resolution_${tag}.json"
    ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $configPath

    return [pscustomobject]@{
        path = $configPath
        csv = Join-Path $BaseDir "vela/pn2d_iv_resolution_${tag}.csv"
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

function Export-AcceptedRows([string]$candidateCsv, [string]$outputCsv) {
    $rows = Import-Csv $candidateCsv
    $accepted = @(
        $rows | Where-Object {
            $_.handoff_stage -eq "newton" -and
            $_.newton_iterations -ne "" -and
            [int]$_.newton_iterations -gt 0 -and
            $_.bias_V -ne "" -and
            $_.current_total_A_per_um -ne ""
        }
    )

    if ($accepted.Count -gt 0) {
        $accepted | Export-Csv -NoTypeInformation -Encoding ASCII $outputCsv
    }

    return $accepted
}

function Compare-Case([string]$candidateCsv, [string]$filteredCsv, [double]$step) {
    $referenceCsv = Join-Path $BaseDir "reference_curves/pn2d_iv_reference.csv"
    $stepTag = ([string]$step).Replace(".", "p")
    $reportJson = Join-Path $BaseDir "reports/pn2d_iv_resolution_${stepTag}_comparison.json"
    $reportMd = Join-Path $BaseDir "reports/pn2d_iv_resolution_${stepTag}_comparison.md"

    python scripts\compare_reference_curves.py `
        --reference $referenceCsv `
        --candidate $filteredCsv `
        --output-json $reportJson `
        --output-md $reportMd `
        --kind iv `
        --candidate-column current_total_A_per_um `
        --candidate-scale -1.0 `
        --bias-min $BiasMin `
        --bias-max $BiasMax | Out-Null

    $report = Get-Content $reportJson -Raw | ConvertFrom-Json
    $refRows = Import-Csv $referenceCsv
    $candRows = Import-Csv $candidateCsv

    $ratio025 = [double]::NaN
    $ratio027 = [double]::NaN
    $ratio029 = [double]::NaN

    foreach ($bias in @(0.25, 0.27, 0.29)) {
        $refAtBias = Get-InterpolatedValue $refRows "current_total" $bias 1.0
        $candAtBias = Get-InterpolatedValue $candRows "current_total_A_per_um" $bias -1.0
        $ratio = [double]::NaN
        if ($refAtBias -ne 0.0 -and -not [double]::IsNaN($refAtBias) -and -not [double]::IsNaN($candAtBias)) {
            $ratio = [math]::Abs($candAtBias / $refAtBias)
        }
        if ($bias -eq 0.25) { $ratio025 = $ratio }
        if ($bias -eq 0.27) { $ratio027 = $ratio }
        if ($bias -eq 0.29) { $ratio029 = $ratio }
    }

    return [pscustomobject]@{
        points = [int]$report.iv.points_compared
        orders = [double]$report.iv.orders_of_magnitude
        maxRelativeError = [double]$report.iv.max_relative_error
        ratio025 = $ratio025
        ratio027 = $ratio027
        ratio029 = $ratio029
    }
}

$ivBaseConfig = Join-Path $BaseDir "vela/simulation_iv.json"
$baseCfg = Get-Content $ivBaseConfig -Raw | ConvertFrom-Json
$promotedStep = [double]$baseCfg.sweep.step

$cases = @(
    @{ tag = "promoted"; step = $promotedStep },
    @{ tag = "step0p02"; step = 0.02 },
    @{ tag = "step0p01"; step = 0.01 }
)

$rows = @()

foreach ($case in $cases) {
    $cfg = New-ResolutionConfig ([double]$case.step) $ivBaseConfig ([string]$case.tag)
    $run = Invoke-CaseRun $cfg.path

    $acceptedCount = 0
    $orders = [double]::NaN
    $points = [double]::NaN
    $maxRelativeError = [double]::NaN
    $ratio025 = [double]::NaN
    $ratio027 = [double]::NaN
    $ratio029 = [double]::NaN

    if ($run.status -eq "ok" -and (Test-Path $cfg.csv)) {
        $filteredCsv = Join-Path $BaseDir "vela/pn2d_iv_resolution_${($case.tag)}_accepted.csv"
        $accepted = Export-AcceptedRows $cfg.csv $filteredCsv
        $acceptedCount = $accepted.Count

        if ($acceptedCount -gt 0 -and (Test-Path $filteredCsv)) {
            $cmp = Compare-Case $cfg.csv $filteredCsv ([double]$case.step)
            $points = $cmp.points
            $orders = $cmp.orders
            $maxRelativeError = $cmp.maxRelativeError
            $ratio025 = $cmp.ratio025
            $ratio027 = $cmp.ratio027
            $ratio029 = $cmp.ratio029
        } else {
            $run.status = "no_accepted_rows"
        }
    }

    $rows += [pscustomobject]@{
        case = [string]$case.tag
        step_V = [double]$case.step
        status = [string]$run.status
        accepted_rows = $acceptedCount
        compared_points = $points
        orders_of_magnitude = $orders
        max_relative_error = $maxRelativeError
        ratio_at_0p25V = $ratio025
        ratio_at_0p27V = $ratio027
        ratio_at_0p29V = $ratio029
        csv_file = $cfg.csv
        config = $cfg.path
    }

    Write-Host "done $($case.tag) step=$($case.step) status=$($run.status) orders=$orders"
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
