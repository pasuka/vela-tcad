param(
    [string]$BaseConfig = "build/pn2d_recomb_gate/vela/simulation_bv.json",
    [string]$ReferenceCsv = "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv",
    [string]$OutputSummary = "build/pn2d_recomb_gate/vela/pn2d_bv_ct_quick6_summary.csv"
)

$ErrorActionPreference = "Stop"

$cases = @(
    @{ name = "q_mu0p89_a0p89";       muScale = 0.89; nrefScale = 1.00; alphaScale = 0.89 },
    @{ name = "q_mu0p89_a0p90";       muScale = 0.89; nrefScale = 1.00; alphaScale = 0.90 },
    @{ name = "q_mu0p91_a0p90";       muScale = 0.91; nrefScale = 1.00; alphaScale = 0.90 },
    @{ name = "q_mu0p90_nr1p02_a0p90"; muScale = 0.90; nrefScale = 1.02; alphaScale = 0.90 },
    @{ name = "q_mu0p90_a0p91";       muScale = 0.90; nrefScale = 1.00; alphaScale = 0.91 },
    @{ name = "q_mu0p90_a0p89";       muScale = 0.90; nrefScale = 1.00; alphaScale = 0.89 }
)

$baseElectronMuMinCm2 = 52.2
$baseHoleMuMinCm2 = 44.9
$baseElectronNrefCm3 = 9.68e16
$baseHoleNrefCm3 = 2.23e17
$baseElectronAlpha = 0.68
$baseHoleAlpha = 0.70

$base = Get-Content $BaseConfig -Raw | ConvertFrom-Json
$ref = [double]((Import-Csv $ReferenceCsv | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1).current_total)
$outputDir = Split-Path -Parent $OutputSummary
$rows = @()

foreach ($c in $cases) {
    $tag = $c.name
    $csvName = "pn2d_bv_${tag}.csv"
    $cfgPath = Join-Path $outputDir "simulation_bv_${tag}.json"

    $cfg = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
    $cfg.output_csv = $csvName
    $cfg.solver.bandgap_narrowing = "none"
    $cfg.solver.mobility = @{
        model = "caughey_thomas"
        electron_mu_min_m2_V_s = $baseElectronMuMinCm2 * $c.muScale
        hole_mu_min_m2_V_s = $baseHoleMuMinCm2 * $c.muScale
        electron_nref_m3 = $baseElectronNrefCm3 * $c.nrefScale
        hole_nref_m3 = $baseHoleNrefCm3 * $c.nrefScale
        electron_alpha = $baseElectronAlpha * $c.alphaScale
        hole_alpha = $baseHoleAlpha * $c.alphaScale
    }
    ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $cfgPath

    $status = "ok"
    $converged = ""
    $points = ""
    $err = ""

    try {
        $raw = .\build\vela_example_runner.exe --config $cfgPath | Out-String
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

    $total = [double]::NaN
    $ratio = [double]::NaN
    $orders = [double]::NaN

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
            }
        } else {
            $status = "missing_csv"
        }
    }

    $rows += [pscustomobject]@{
        case = $tag
        mu_scale = $c.muScale
        nref_scale = $c.nrefScale
        alpha_scale = $c.alphaScale
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

    Write-Host "done $tag status=$status orders=$orders"
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } | Sort-Object orders | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
