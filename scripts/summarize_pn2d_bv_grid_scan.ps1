param(
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

$ref = [double]((Import-Csv $ReferenceCsv | Where-Object { [double]$_.bias_V -eq 0.05 } | Select-Object -First 1).current_total)
$rows = @()

foreach ($m in $mobilityModels) {
    foreach ($b in $bgnModels) {
        $tag = "m_${m}__bgn_${b}"
        $cfgPath = "build/pn2d_recomb_gate/vela/simulation_bv_${tag}.json"
        $csvPath = "build/pn2d_recomb_gate/vela/pn2d_bv_${tag}.csv"

        $status = "ok"
        $converged = ""
        $points = ""
        $total = [double]::NaN
        $ratio = [double]::NaN
        $orders = [double]::NaN
        $eTotal = [double]::NaN
        $eDrift = [double]::NaN
        $eDiff = [double]::NaN
        $hTotal = [double]::NaN
        $hDiff = [double]::NaN

        if (-not (Test-Path $cfgPath)) {
            $status = "missing_config"
        } elseif (-not (Test-Path $csvPath)) {
            $status = "missing_csv"
        } else {
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

            try {
                $allRows = Import-Csv $csvPath
                if ($allRows.Count -gt 0) {
                    $last = $allRows[-1]
                    $converged = [string]$last.converged
                    $points = [string]$allRows.Count
                }
            } catch {
                # Metadata is optional; keep defaults.
            }
        }

        $rows += [pscustomobject]@{
            mobility = $m
            bgn = $b
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
            error = ""
            config = $cfgPath
            csv_file = $csvPath
        }
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary

$ranked = $rows |
    Where-Object { $_.status -eq "ok" -and -not [double]::IsNaN([double]$_.orders) } |
    Sort-Object orders

Write-Host "summary=$OutputSummary"
Write-Host "top10:"
$ranked | Select-Object mobility, bgn, status, points, total_A_per_um, ratio_vs_ref, orders | Format-Table -AutoSize
