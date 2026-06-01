 $refine6Cases = @'
[
  {"name":"ref6_mu0p88_a0p89","muScale":0.88,"nrefScale":1.00,"alphaScale":0.89},
  {"name":"ref6_mu0p89_a0p89","muScale":0.89,"nrefScale":1.00,"alphaScale":0.89},
  {"name":"ref6_mu0p90_a0p89","muScale":0.90,"nrefScale":1.00,"alphaScale":0.89},
  {"name":"ref6_mu0p88_a0p90","muScale":0.88,"nrefScale":1.00,"alphaScale":0.90},
  {"name":"ref6_mu0p89_a0p90","muScale":0.89,"nrefScale":1.00,"alphaScale":0.90},
  {"name":"ref6_mu0p90_a0p90","muScale":0.90,"nrefScale":1.00,"alphaScale":0.90}
]
'@

& $PSScriptRoot\scan_pn2d_bv_ct_base.ps1 `
    -BaseConfig "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json" `
    -ReferenceCsv "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv" `
    -OutputSummary "build/pn2d_recomb_gate/vela/pn2d_bv_ct_refine6_summary.csv" `
    -CaseMode ExplicitList `
    -CasesJson $refine6Cases `
    -SecondsPerCase 90 `
    -UseStartProcess
