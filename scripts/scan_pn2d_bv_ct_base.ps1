param(
    [string]$BaseConfig = "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json",
    [string]$ReferenceCsv = "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv",
    [string]$OutputSummary = "build/pn2d_recomb_gate/vela/pn2d_bv_ct_scan_summary.csv",
    [ValidateSet("Grid", "ExplicitList")]
    [string]$CaseMode = "Grid",
    [double[]]$MuScales = @(1.0),
    [double[]]$NrefScales = @(1.0),
    [double[]]$AlphaScales = @(1.0),
    [string]$CasesJson = "",
    [string]$TagPrefix = "ct_scan",
    [int]$SecondsPerCase = 0,
    [switch]$UseStartProcess,
    [switch]$IncludeComponentCurrents,
    [string]$BandgapNarrowing = "",
    [string[]]$Recombination = @()
)

$ErrorActionPreference = "Stop"

$baseElectronMuMinCm2 = 52.2
$baseHoleMuMinCm2 = 44.9
$baseElectronNrefCm3 = 9.68e16
$baseHoleNrefCm3 = 2.23e17
$baseElectronAlpha = 0.68
$baseHoleAlpha = 0.70

$base = Get-Content $BaseConfig -Raw | ConvertFrom-Json
$ref = [double]((Import-Csv $ReferenceCsv | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1).current_total)

function Convert-ToCaseList {
    param(
        [string]$Mode,
        [double[]]$Mu,
        [double[]]$Nref,
        [double[]]$Alpha,
        [string]$JsonCases,
        [string]$Prefix
    )

    if ($Mode -eq "ExplicitList") {
        if ([string]::IsNullOrWhiteSpace($JsonCases)) {
            throw "CasesJson is required when CaseMode is ExplicitList"
        }

        $items = $JsonCases | ConvertFrom-Json
        return @($items | ForEach-Object {
            [pscustomobject]@{
                name = [string]$_.name
                muScale = [double]$_.muScale
                nrefScale = [double]$_.nrefScale
                alphaScale = [double]$_.alphaScale
            }
        })
    }

    $cases = @()
    foreach ($muScale in $Mu) {
        foreach ($nrefScale in $Nref) {
            foreach ($alphaScale in $Alpha) {
                $name = "${Prefix}_mu${muScale}_nr${nrefScale}_a${alphaScale}".Replace('.', 'p')
                $cases += [pscustomobject]@{
                    name = $name
                    muScale = $muScale
                    nrefScale = $nrefScale
                    alphaScale = $alphaScale
                }
            }
        }
    }

    return $cases
}

function Set-MobilityConfig {
    param(
        $Config,
        [double]$MuScale,
        [double]$NrefScale,
        [double]$AlphaScale
    )

    $Config.solver.mobility = @{
        model = "caughey_thomas"
        electron_mu_min_m2_V_s = $baseElectronMuMinCm2 * $MuScale
        hole_mu_min_m2_V_s = $baseHoleMuMinCm2 * $MuScale
        electron_nref_m3 = $baseElectronNrefCm3 * $NrefScale
        hole_nref_m3 = $baseHoleNrefCm3 * $NrefScale
        electron_alpha = $baseElectronAlpha * $AlphaScale
        hole_alpha = $baseHoleAlpha * $AlphaScale
    }
}

function Invoke-ExampleRunner {
    param(
        [string]$ConfigPath,
        [string]$Tag,
        [int]$TimeoutSeconds,
        [switch]$StartProcessMode
    )

    if ($StartProcessMode) {
        $stdoutPath = Join-Path (Split-Path -Parent $ConfigPath) "scan_tmp_${Tag}.out"
        $stderrPath = Join-Path (Split-Path -Parent $ConfigPath) "scan_tmp_${Tag}.err"
        Remove-Item $stdoutPath -ErrorAction SilentlyContinue
        Remove-Item $stderrPath -ErrorAction SilentlyContinue

        $proc = Start-Process -FilePath ".\build\vela_example_runner.exe" `
            -ArgumentList @("--config", $ConfigPath) `
            -NoNewWindow -PassThru `
            -RedirectStandardOutput $stdoutPath `
            -RedirectStandardError $stderrPath

        $finished = $true
        try {
            Wait-Process -Id $proc.Id -Timeout $TimeoutSeconds
        } catch {
            $finished = $false
        }

        if (-not $finished) {
            Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
            return [pscustomobject]@{ status = "timeout"; converged = ""; points = ""; error = "timeout_after_${TimeoutSeconds}s" }
        }

        $raw = Get-Content $stdoutPath -Raw -ErrorAction SilentlyContinue
        if ($raw) {
            try {
                $j = $raw | ConvertFrom-Json
                return [pscustomobject]@{ status = "ok"; converged = [string]$j.converged; points = [string]$j.points; error = "" }
            } catch {
                return [pscustomobject]@{ status = "bad_runner_output"; converged = ""; points = ""; error = "stdout_not_json" }
            }
        }

        return [pscustomobject]@{ status = "ok"; converged = ""; points = ""; error = "" }
    }

    try {
        $raw = .\build\vela_example_runner.exe --config $ConfigPath | Out-String
        if ($raw) {
            try {
                $j = $raw | ConvertFrom-Json
                return [pscustomobject]@{ status = "ok"; converged = [string]$j.converged; points = [string]$j.points; error = "" }
            } catch {
                return [pscustomobject]@{ status = "bad_runner_output"; converged = ""; points = ""; error = "stdout_not_json" }
            }
        }
        return [pscustomobject]@{ status = "ok"; converged = ""; points = ""; error = "" }
    } catch {
        return [pscustomobject]@{ status = "runner_failed"; converged = ""; points = ""; error = $_.Exception.Message }
    }
}

$outputDir = Split-Path -Parent $OutputSummary
if (-not [string]::IsNullOrWhiteSpace($outputDir)) {
    New-Item -ItemType Directory -Force -Path $outputDir | Out-Null
}

$cases = Convert-ToCaseList -Mode $CaseMode -Mu $MuScales -Nref $NrefScales -Alpha $AlphaScales -JsonCases $CasesJson -Prefix $TagPrefix
$rows = @()

foreach ($case in $cases) {
    $tag = $case.name
    $csvName = "pn2d_bv_${tag}.csv"
    $cfgPath = Join-Path $outputDir "simulation_bv_${tag}.json"

    $cfg = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
    $cfg.output_csv = $csvName
    Set-MobilityConfig -Config $cfg -MuScale $case.muScale -NrefScale $case.nrefScale -AlphaScale $case.alphaScale

    if (-not [string]::IsNullOrWhiteSpace($BandgapNarrowing)) {
        $cfg.solver.bandgap_narrowing = $BandgapNarrowing
    }
    if ($Recombination.Count -gt 0) {
        $cfg.solver.recombination = $Recombination
    }

    ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $cfgPath

    $run = Invoke-ExampleRunner -ConfigPath $cfgPath -Tag $tag -TimeoutSeconds $SecondsPerCase -StartProcessMode:$UseStartProcess
    $status = $run.status
    $converged = $run.converged
    $points = $run.points
    $err = $run.error

    $total = [double]::NaN
    $ratio = [double]::NaN
    $orders = [double]::NaN
    $eTotal = [double]::NaN
    $eDrift = [double]::NaN
    $eDiff = [double]::NaN
    $hTotal = [double]::NaN
    $hDiff = [double]::NaN

    if ($status -eq "ok") {
        $csvPath = Join-Path $outputDir $csvName
        if (Test-Path $csvPath) {
            $row = Import-Csv $csvPath | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1
            if ($null -eq $row) {
                $status = "missing_0p05_row"
            } else {
                $total = [double]$row.current_total_A_per_um
                if ($ref -ne 0.0) {
                    $ratio = [math]::Abs($total / $ref)
                    if ($ratio -gt 0.0) {
                        $orders = [math]::Abs([math]::Log10($ratio))
                    }
                }

                if ($IncludeComponentCurrents) {
                    $eTotal = [double]$row.current_electron_A_per_um
                    $eDrift = [double]$row.current_electron_drift_A_per_um
                    $eDiff = [double]$row.current_electron_diffusion_A_per_um
                    $hTotal = [double]$row.current_hole_A_per_um
                    $hDiff = [double]$row.current_hole_diffusion_A_per_um
                }
            }
        } else {
            $status = "missing_csv"
        }
    }

    $rowObject = [ordered]@{
        case = $tag
        mu_scale = $case.muScale
        nref_scale = $case.nrefScale
        alpha_scale = $case.alphaScale
        status = $status
        converged = $converged
        points = $points
        total_A_per_um = $total
        ratio_vs_ref = $ratio
        orders = $orders
        error = $err
        config = $cfgPath
        csv_file = Join-Path $outputDir $csvName
    }

    if ($IncludeComponentCurrents) {
        $rowObject.e_total_A_per_um = $eTotal
        $rowObject.e_drift_A_per_um = $eDrift
        $rowObject.e_diff_A_per_um = $eDiff
        $rowObject.h_total_A_per_um = $hTotal
        $rowObject.h_diff_A_per_um = $hDiff
    }

    $rows += [pscustomobject]$rowObject
    Write-Host "done $tag status=$status orders=$orders"
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$ranked = $rows | Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } | Sort-Object orders
Write-Host "summary=$OutputSummary"
$ranked | Select-Object -First 5 case,mu_scale,nref_scale,alpha_scale,orders,ratio_vs_ref,total_A_per_um | Format-Table -AutoSize