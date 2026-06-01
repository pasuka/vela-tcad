param(
    [string]$BaseConfig = "build/pn2d_recomb_gate/vela/simulation_bv.json",
    [string]$ReferenceCsv = "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv",
    [string]$OutputSummary = "build/pn2d_recomb_gate/vela/pn2d_bv_grid_scan_summary.csv"
)

$ErrorActionPreference = "Stop"

$mobilityModels = @(
    "constant",
    "caughey_thomas",
    "caughey_thomas_field",
    "caughey_thomas_surface",
    "caughey_thomas_field_surface"
)

$bgnModels = @("none", "slotboom")

$base = Get-Content $BaseConfig -Raw | ConvertFrom-Json
$ref = [double]((Import-Csv $ReferenceCsv | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1).current_total)

$rows = @()

foreach ($m in $mobilityModels) {
    foreach ($b in $bgnModels) {
        $tag = "m_${m}__bgn_${b}"
        $csvName = "pn2d_bv_${tag}.csv"
        $cfgPath = "build/pn2d_recomb_gate/vela/simulation_bv_${tag}.json"

        $cfg = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
        $cfg.output_csv = $csvName
        $cfg.solver.mobility.model = $m
        $cfg.solver.bandgap_narrowing = $b
        ($cfg | ConvertTo-Json -Depth 100) | Set-Content -Encoding utf8 $cfgPath

        $status = "ok"
        $converged = ""
        $points = ""
        $err = ""

        try {
            $raw = .\build\vela_example_runner.exe --config $cfgPath | Out-String
            if ($raw) {
                try {
                    $json = $raw | ConvertFrom-Json
                    $converged = [string]$json.converged
                    $points = [string]$json.points
                } catch {
                    # Keep default metadata if parse fails.
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
                if ($null -ne $row) {
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
                } else {
                    $status = "missing_0p05_row"
                }
            } else {
                $status = "missing_csv"
            }
        }

        $rows += [pscustomobject]@{
            mobility        = $m
            bgn             = $b
            status          = $status
            converged       = $converged
            points          = $points
            total_A_per_um  = $total
            ratio_vs_ref    = $ratio
            orders          = $orders
            e_total_A_per_um = $eTotal
            e_drift_A_per_um = $eDrift
            e_diff_A_per_um = $eDiff
            h_total_A_per_um = $hTotal
            h_diff_A_per_um = $hDiff
            error           = $err
            config          = $cfgPath
            csv_file        = "build/pn2d_recomb_gate/vela/$csvName"
        }

        Write-Host "done $tag status=$status orders=$orders"
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary

$top = $rows |
    Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } |
    Sort-Object orders

Write-Host "summary=$OutputSummary"
Write-Host "top3:"
$top | Select-Object -First 3 mobility, bgn, orders, ratio_vs_ref, total_A_per_um, e_total_A_per_um, h_total_A_per_um | Format-Table -AutoSize
