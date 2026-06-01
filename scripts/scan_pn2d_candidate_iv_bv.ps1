param(
    [string]$BaseDir = "build/pn2d_current_review",
    [string]$Case = "q_mu0p89_a0p89",
    [string]$OutputSummary = ""
)

$ErrorActionPreference = "Stop"

if ($OutputSummary -eq "") {
    $OutputSummary = Join-Path $BaseDir "vela/pn2d_candidate_iv_bv_summary.csv"
}

$candidates = @{
    q_mu0p89_a0p89 = @{
        muScale = 0.89
        nrefScale = 1.00
        alphaScale = 0.89
    }
}

if (-not $candidates.ContainsKey($Case)) {
    throw "Unknown candidate '$Case'. Known candidates: $($candidates.Keys -join ', ')"
}

$baseElectronMuMinCm2 = 52.2
$baseHoleMuMinCm2 = 44.9
$baseElectronNrefCm3 = 9.68e16
$baseHoleNrefCm3 = 2.23e17
$baseElectronAlpha = 0.68
$baseHoleAlpha = 0.70

function Get-DoubleField($row, [string]$name) {
    return [double]$row.$name
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

function New-CandidateConfig([string]$kind, [string]$baseConfigPath, [string]$csvName, [hashtable]$candidate) {
    $cfg = Get-Content $baseConfigPath -Raw | ConvertFrom-Json
    $cfg.output_csv = $csvName
    $cfg.solver.bandgap_narrowing = "none"
    $cfg.solver.mobility = @{
        model = "caughey_thomas"
        electron_mu_min_m2_V_s = $baseElectronMuMinCm2 * $candidate.muScale
        hole_mu_min_m2_V_s = $baseHoleMuMinCm2 * $candidate.muScale
        electron_nref_m3 = $baseElectronNrefCm3 * $candidate.nrefScale
        hole_nref_m3 = $baseHoleNrefCm3 * $candidate.nrefScale
        electron_alpha = $baseElectronAlpha * $candidate.alphaScale
        hole_alpha = $baseHoleAlpha * $candidate.alphaScale
    }

    $configPath = Join-Path $BaseDir "vela/simulation_${kind}_${Case}.json"
    ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $configPath
    return $configPath
}

function Invoke-CandidateRun([string]$kind, [string]$configPath) {
    $status = "ok"
    $converged = ""
    $points = ""
    $err = ""
    try {
        $raw = .\build\vela_example_runner.exe --config $configPath | Out-String
        if ($raw) {
            try {
                $j = $raw | ConvertFrom-Json
                $converged = [string]$j.converged
                $points = [string]$j.points
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
        converged = $converged
        points = $points
        error = $err
    }
}

function Compare-Candidate([string]$kind, [string]$referenceCsv, [string]$candidateCsv) {
    $reportJson = Join-Path $BaseDir "reports/pn2d_${kind}_${Case}_comparison.json"
    $reportMd = Join-Path $BaseDir "reports/pn2d_${kind}_${Case}_comparison.md"
    if ($kind -eq "iv") {
        python scripts\compare_reference_curves.py `
            --reference $referenceCsv `
            --candidate $candidateCsv `
            --output-json $reportJson `
            --output-md $reportMd `
            --kind iv `
            --candidate-column current_total_A_per_um `
            --candidate-scale -1.0 `
            --bias-min 0.2 `
            --bias-max 0.3 | Out-Null
        $targetBias = 0.29
        $candidateScale = -1.0
    } else {
        python scripts\compare_reference_curves.py `
            --reference $referenceCsv `
            --candidate $candidateCsv `
            --output-json $reportJson `
            --output-md $reportMd `
            --kind iv `
            --candidate-column current_total_A_per_um `
            --candidate-scale 1.0 `
            --bias-min 0.05 `
            --bias-max 0.05 | Out-Null
        $targetBias = 0.05
        $candidateScale = 1.0
    }

    $report = Get-Content $reportJson -Raw | ConvertFrom-Json
    $refRows = Import-Csv $referenceCsv
    $candRows = Import-Csv $candidateCsv
    $refValue = Get-InterpolatedValue $refRows "current_total" $targetBias 1.0
    $candValue = Get-InterpolatedValue $candRows "current_total_A_per_um" $targetBias $candidateScale
    $ratio = [double]::NaN
    if ($refValue -ne 0.0 -and -not [double]::IsNaN($refValue) -and -not [double]::IsNaN($candValue)) {
        $ratio = [math]::Abs($candValue / $refValue)
    }

    return [pscustomobject]@{
        orders = [double]$report.iv.orders_of_magnitude
        ratio = $ratio
        total = $candValue
    }
}

$simulation_iv = Join-Path $BaseDir "vela/simulation_iv.json"
$simulation_bv = Join-Path $BaseDir "vela/simulation_bv.json"
$candidate = $candidates[$Case]
$rows = @()

foreach ($kind in @("iv", "bv")) {
    $baseConfig = if ($kind -eq "iv") { $simulation_iv } else { $simulation_bv }
    $referenceCsv = Join-Path $BaseDir "reference_curves/pn2d_${kind}_reference.csv"
    $csvName = "pn2d_${kind}_${Case}.csv"
    $candidateCsv = Join-Path $BaseDir "vela/$csvName"
    $configPath = New-CandidateConfig $kind $baseConfig $csvName $candidate
    $run = Invoke-CandidateRun $kind $configPath

    $orders = [double]::NaN
    $ratio = [double]::NaN
    $total = [double]::NaN
    if ($run.status -eq "ok" -and (Test-Path $candidateCsv)) {
        $cmp = Compare-Candidate $kind $referenceCsv $candidateCsv
        $orders = $cmp.orders
        $ratio = $cmp.ratio
        $total = $cmp.total
    }

    $rows += [pscustomobject]@{
        case = $Case
        kind = $kind
        status = $run.status
        converged = $run.converged
        points = $run.points
        orders = $orders
        ratio_vs_ref = $ratio
        total_A_per_um = $total
        error = $run.error
        csv_file = $candidateCsv
        config = $configPath
    }

    Write-Host "done $kind $Case status=$($run.status) orders=$orders"
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
