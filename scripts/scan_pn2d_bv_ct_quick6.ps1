$quickCases = @'
[
  {"name":"q_mu0p89_a0p89","muScale":0.89,"nrefScale":1.00,"alphaScale":0.89},
  {"name":"q_mu0p89_a0p90","muScale":0.89,"nrefScale":1.00,"alphaScale":0.90},
  {"name":"q_mu0p91_a0p90","muScale":0.91,"nrefScale":1.00,"alphaScale":0.90},
  {"name":"q_mu0p90_nr1p02_a0p90","muScale":0.90,"nrefScale":1.02,"alphaScale":0.90},
  {"name":"q_mu0p90_a0p91","muScale":0.90,"nrefScale":1.00,"alphaScale":0.91},
  {"name":"q_mu0p90_a0p89","muScale":0.90,"nrefScale":1.00,"alphaScale":0.89}
]
'@

& $PSScriptRoot\scan_pn2d_bv_ct_base.ps1 `
    -BaseConfig "build/pn2d_recomb_gate/vela/simulation_bv.json" `
    -ReferenceCsv "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv" `
    -OutputSummary "build/pn2d_recomb_gate/vela/pn2d_bv_ct_quick6_summary.csv" `
    -CaseMode ExplicitList `
    -CasesJson $quickCases
