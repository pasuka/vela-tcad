param(
    [string]$BaseConfig = "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json",
    [string]$ReferenceCsv = "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv",
    [string]$OutputSummary = "build/pn2d_recomb_gate/vela/pn2d_bv_ct_local_scan_summary.csv"
)

$ErrorActionPreference = "Stop"

$muScales = @(0.6, 0.8, 1.0, 1.2)
$nrefScales = @(0.5, 1.0, 2.0)
$alphaScales = @(0.9, 1.0, 1.1)

# Defaults in include/vela/physics/MobilityModel.h, converted to unit_scaling input units.
$baseElectronMuMinCm2 = 52.2
$baseHoleMuMinCm2 = 44.9
$baseElectronNrefCm3 = 9.68e16
$baseHoleNrefCm3 = 2.23e17
$baseElectronAlpha = 0.68
$baseHoleAlpha = 0.70

$base = Get-Content $BaseConfig -Raw | ConvertFrom-Json
$ref = [double]((Import-Csv $ReferenceCsv | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1).current_total)

$rows = @()

foreach ($muScale in $muScales) {
    foreach ($nrefScale in $nrefScales) {
        foreach ($alphaScale in $alphaScales) {
            $tag = "ct_local_mu${muScale}_nr${nrefScale}_a${alphaScale}".Replace('.', 'p')
            $csvName = "pn2d_bv_${tag}.csv"
            $cfgPath = "build/pn2d_recomb_gate/vela/simulation_bv_${tag}.json"

            $cfg = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
            $cfg.output_csv = $csvName
            $cfg.solver.bandgap_narrowing = "none"
            $cfg.solver.recombination = @("none")
            $cfg.solver.mobility = @{
                model = "caughey_thomas"
                electron_mu_min_m2_V_s = $baseElectronMuMinCm2 * $muScale
                hole_mu_min_m2_V_s = $baseHoleMuMinCm2 * $muScale
                electron_nref_m3 = $baseElectronNrefCm3 * $nrefScale
                hole_nref_m3 = $baseHoleNrefCm3 * $nrefScale
                electron_alpha = $baseElectronAlpha * $alphaScale
                hole_alpha = $baseHoleAlpha * $alphaScale
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
                        # Keep metadata empty if parse fails.
                    }
                }
            } catch {
                $status = "runner_failed"
                $err = $_.Exception.Message
            }

            $total = [double]::NaN
            $ratio = [double]::NaN
            $orders = [double]::NaN
            $eTotal = [double]::NaN
            $eDrift = [double]::NaN
            $eDiff = [double]::NaN
            $hTotal = [double]::NaN
            $hDiff = [double]::NaN

            if ($status -eq "ok") {
                $csvPath = "build/pn2d_recomb_gate/vela/$csvName"
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
                        $eTotal = [double]$row.current_electron_A_per_um
                        $eDrift = [double]$row.current_electron_drift_A_per_um
                        $eDiff = [double]$row.current_electron_diffusion_A_per_um
                        $hTotal = [double]$row.current_hole_A_per_um
                        $hDiff = [double]$row.current_hole_diffusion_A_per_um
                    }
                } else {
                    $status = "missing_csv"
                }
            }

            $rows += [pscustomobject]@{
                mu_scale = $muScale
                nref_scale = $nrefScale
                alpha_scale = $alphaScale
                status = $status
                converged = $converged
                points = $points
                total_A_per_um = $total
                ratio_vs_ref = $ratio
                orders = $orders
                e_total_A_per_um = $eTotal
                e_drift_A_per_um = $eDrift
                e_diff_A_per_um = $eDiff
                h_total_A_per_um = $hTotal
                h_diff_A_per_um = $hDiff
                error = $err
                config = $cfgPath
                csv_file = "build/pn2d_recomb_gate/vela/$csvName"
            }

            Write-Host "done $tag status=$status orders=$orders"
        }
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary

$best = $rows |
    Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } |
    Sort-Object orders

Write-Host "summary=$OutputSummary"
Write-Host "top5:"
$best | Select-Object -First 5 mu_scale, nref_scale, alpha_scale, orders, ratio_vs_ref, total_A_per_um, e_total_A_per_um, h_total_A_per_um | Format-Table -AutoSize
