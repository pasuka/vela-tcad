* Sentaurus Training BVmethods, Workbench node 7: continuation method.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "continuation_des.tdr"
  Parameter = "models.par"
  Current   = "continuation_des.plt"
  Output    = "continuation_des.log"
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
  Recombination(SRH(DopingDep) Band2Band(E2) Avalanche(Eparallel))
  Fermi
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
  AvalDerivatives
  ErrRef(Electron)=1e10
  ErrRef(Hole)=1e10
  BreakCriteria { Current(Contact="drain" AbsVal=1.443e-3) }
  Transient=BE
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled { Poisson Electron Hole }
  Continuation(
    Name="drain"
    InitialVStep=0.1
    Increment=1.41 Decrement=2
    MaxVStep=0.5
    MinVoltage=0 MaxVoltage=100
    MinCurrent=0 MaxCurrent=1.443e-3
    Iadapt=1e-8
  ) { Coupled { Poisson Electron Hole } }
}
