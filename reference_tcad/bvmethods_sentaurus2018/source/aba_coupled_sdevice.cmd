* Sentaurus Training BVmethods, Workbench node 4: coupled ABA postprocessing.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "aba_coupled_des.tdr"
  Parameter = "models.par"
  Current   = "aba_coupled_des.plt"
  Output    = "aba_coupled_des.log"
}

Electrode {
  { Name="drain"     Voltage=0.0 }
  { Name="source"    Voltage=0.0 }
  { Name="gate"      Voltage=0.0 Barrier=-0.55 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(DopingDep HighFieldsaturation(GradQuasiFermi) Enormal)
  Recombination(
    SRH(DopingDep)
    Band2Band(E2)
    Avalanche(Eparallel)
  )
  Fermi
}

CurrentPlot {
  AvalancheGeneration(Integrate(Semiconductor))
}

Plot {
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  ElectricField/Vector Potential SpaceCharge
  Doping
  SRH Band2Band
  AvalancheGeneration eAvalancheGeneration hAvalancheGeneration
  eIonIntegral hIonIntegral MeanIonIntegral eAlphaAvalanche hAlphaAvalanche
}

Math {
  Extrapolate
  Iterations=20
  Notdamped=100
  RelErrControl
  ComputeIonizationIntegrals
  BreakAtIonIntegral(3 1.)
  AvalPostProcessing
  Transient=BE
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled { Poisson Electron Hole }
  Quasistationary(
    InitialStep=0.0001 Increment=1.41
    MinStep=1e-7 MaxStep=0.025
    Goal { Name="drain" Voltage=100.0 }
  ) { Coupled { Poisson Electron Hole } }
}
