$refine2Cases = @'
[
  {"name":"ref2_mu0p88_nr0p95_a0p88","muScale":0.88,"nrefScale":0.95,"alphaScale":0.88},
  {"name":"ref2_mu0p88_nr1p00_a0p88","muScale":0.88,"nrefScale":1.00,"alphaScale":0.88},
  {"name":"ref2_mu0p88_nr1p05_a0p88","muScale":0.88,"nrefScale":1.05,"alphaScale":0.88},
  {"name":"ref2_mu0p90_nr0p95_a0p90","muScale":0.90,"nrefScale":0.95,"alphaScale":0.90},
  {"name":"ref2_mu0p90_nr1p00_a0p90","muScale":0.90,"nrefScale":1.00,"alphaScale":0.90},
  {"name":"ref2_mu0p90_nr1p05_a0p90","muScale":0.90,"nrefScale":1.05,"alphaScale":0.90},
  {"name":"ref2_mu0p92_nr0p95_a0p92","muScale":0.92,"nrefScale":0.95,"alphaScale":0.92},
  {"name":"ref2_mu0p92_nr1p00_a0p92","muScale":0.92,"nrefScale":1.00,"alphaScale":0.92},
  {"name":"ref2_mu0p92_nr1p05_a0p92","muScale":0.92,"nrefScale":1.05,"alphaScale":0.92},
  {"name":"ref2_mu0p90_nr1p00_a0p88","muScale":0.90,"nrefScale":1.00,"alphaScale":0.88},
  {"name":"ref2_mu0p90_nr1p00_a0p92","muScale":0.90,"nrefScale":1.00,"alphaScale":0.92},
  {"name":"ref2_mu0p89_nr1p00_a0p89","muScale":0.89,"nrefScale":1.00,"alphaScale":0.89}
]
'@

& $PSScriptRoot\scan_pn2d_bv_ct_base.ps1 `
    -BaseConfig "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json" `
    -ReferenceCsv "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv" `
    -OutputSummary "build/pn2d_recomb_gate/vela/pn2d_bv_ct_refine2_summary.csv" `
    -CaseMode ExplicitList `
    -CasesJson $refine2Cases `
    -SecondsPerCase 90 `
    -UseStartProcess
