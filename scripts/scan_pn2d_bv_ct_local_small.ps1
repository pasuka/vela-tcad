& $PSScriptRoot\scan_pn2d_bv_ct_base.ps1 `
    -BaseConfig "build/pn2d_recomb_gate/vela/simulation_bv_m_caughey_thomas__bgn_none.json" `
    -ReferenceCsv "build/pn2d_recomb_gate/reference_curves/pn2d_bv_reference.csv" `
    -OutputSummary "build/pn2d_recomb_gate/vela/pn2d_bv_ct_local_small_scan_summary.csv" `
    -MuScales @(0.9, 1.0, 1.1) `
    -NrefScales @(0.8, 1.0, 1.2) `
    -AlphaScales @(0.9, 1.0) `
    -TagPrefix "ct_small" `
    -UseStartProcess `
    -SecondsPerCase 120
