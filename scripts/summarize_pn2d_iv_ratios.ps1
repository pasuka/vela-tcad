param(
    [string]$ReferenceCsv = "build/pn2d_tdr_tie_probe/reference_curves/pn2d_iv_reference.csv",
    [string]$CandidateCsv = "build/pn2d_tdr_tie_probe/vela/pn2d_iv.csv",
    [string]$OutputSummary = "build/pn2d_tdr_tie_probe/vela/pn2d_iv_ratio_summary.csv",
    [double]$BiasMin = 0.2,
    [double]$BiasMax = 0.3,
    [double]$CandidateScale = -1.0
)

$ErrorActionPreference = "Stop"

function Get-InterpolatedRow($rows, [double]$bias) {
    $pts = @($rows | Sort-Object { [double]$_.bias_V })
    for ($i = 0; $i -lt $pts.Count; $i++) {
        $b = [double]$pts[$i].bias_V
        if ([math]::Abs($b - $bias) -le [math]::Max([math]::Abs($bias), 1.0) * 1.0e-12) {
            return $pts[$i]
        }
    }

    for ($i = 0; $i -lt $pts.Count - 1; $i++) {
        $b0 = [double]$pts[$i].bias_V
        $b1 = [double]$pts[$i + 1].bias_V
        if ($b0 -le $bias -and $bias -le $b1 -and $b1 -ne $b0) {
            $t = ($bias - $b0) / ($b1 - $b0)
            $out = [ordered]@{ bias_V = $bias }
            foreach ($name in $pts[$i].PSObject.Properties.Name) {
                if ($name -eq "bias_V") {
                    continue
                }

                $v0 = 0.0
                $v1 = 0.0
                if ([double]::TryParse([string]$pts[$i].$name, [ref]$v0) -and
                    [double]::TryParse([string]$pts[$i + 1].$name, [ref]$v1)) {
                    $out[$name] = $v0 + $t * ($v1 - $v0)
                }
            }
            return [pscustomobject]$out
        }
    }

    return $null
}

$refRows = Import-Csv $ReferenceCsv
$candRows = Import-Csv $CandidateCsv
$rows = @()

foreach ($ref in $refRows) {
    $bias = [double]$ref.bias_V
    if ($bias -lt $BiasMin -or $bias -gt $BiasMax) {
        continue
    }

    $cand = Get-InterpolatedRow $candRows $bias
    if ($null -eq $cand) {
        continue
    }

    $refCurrent = [double]$ref.current_total
    $velaTotal = [double]$cand.current_total_A_per_um * $CandidateScale
    $rows += [pscustomobject]@{
        bias_V = $bias
        reference_A = $refCurrent
        vela_total_A_per_um = $velaTotal
        ratio_vs_ref = if ($refCurrent -ne 0.0) { $velaTotal / $refCurrent } else { [double]::NaN }
        current_electron_A_per_um = [double]$cand.current_electron_A_per_um * $CandidateScale
        current_hole_A_per_um = [double]$cand.current_hole_A_per_um * $CandidateScale
        current_electron_drift_A_per_um = [double]$cand.current_electron_drift_A_per_um * $CandidateScale
        current_electron_diffusion_A_per_um = [double]$cand.current_electron_diffusion_A_per_um * $CandidateScale
        current_hole_diffusion_A_per_um = [double]$cand.current_hole_diffusion_A_per_um * $CandidateScale
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
