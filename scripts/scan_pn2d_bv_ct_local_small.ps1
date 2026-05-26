param(
    [string]$BaseConfig = "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json",
    [string]$ReferenceCsv = "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv",
    [string]$OutputSummary = "build/pn2d_recomb_gate/vela/pn2d_bv_ct_local_small_scan_summary.csv",
    [int]$SecondsPerCase = 120
)

$ErrorActionPreference = "Stop"

$muScales = @(0.9, 1.0, 1.1)
$nrefScales = @(0.8, 1.0, 1.2)
$alphaScales = @(0.9, 1.0)

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
            $tag = "ct_small_mu${muScale}_nr${nrefScale}_a${alphaScale}".Replace('.', 'p')
            $csvName = "pn2d_bv_${tag}.csv"
            $cfgPath = "build/pn2d_recomb_gate/vela/simulation_bv_${tag}.json"

            $cfg = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
            $cfg.output_csv = $csvName
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
                $stdoutPath = "build/pn2d_recomb_gate/vela/scan_tmp_${tag}.out"
                $stderrPath = "build/pn2d_recomb_gate/vela/scan_tmp_${tag}.err"
                Remove-Item $stdoutPath -ErrorAction SilentlyContinue
                Remove-Item $stderrPath -ErrorAction SilentlyContinue

                $proc = Start-Process -FilePath ".\\build\\vela_example_runner.exe" `
                    -ArgumentList @("--config", $cfgPath) `
                    -NoNewWindow -PassThru `
                    -RedirectStandardOutput $stdoutPath `
                    -RedirectStandardError $stderrPath

                $finished = $true
                try {
                    Wait-Process -Id $proc.Id -Timeout $SecondsPerCase
                } catch {
                    $finished = $false
                }

                if (-not $finished) {
                    Stop-Process -Id $proc.Id -Force -ErrorAction SilentlyContinue
                    $status = "timeout"
                    $err = "timeout_after_${SecondsPerCase}s"
                } else {
                    $raw = Get-Content $stdoutPath -Raw -ErrorAction SilentlyContinue
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
                }
            } catch {
                $status = "runner_failed"
                $err = $_.Exception.Message
            }

            $total = [double]::NaN
            $ratio = [double]::NaN
            $orders = [double]::NaN
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
                error = $err
                config = $cfgPath
                csv_file = "build/pn2d_recomb_gate/vela/$csvName"
            }
            Write-Host "done $tag status=$status orders=$orders"
        }
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$ranked = $rows | Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } | Sort-Object orders
Write-Host "summary=$OutputSummary"
$ranked | Select-Object -First 5 mu_scale,nref_scale,alpha_scale,orders,ratio_vs_ref,total_A_per_um | Format-Table -AutoSize
