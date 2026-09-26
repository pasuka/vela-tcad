* Sentaurus Training BVmethods, Workbench node 5: external resistor method.
* In 2-D the resistance is interpreted as ohm*um.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "external_resistor_des.tdr"
  Parameter = "models.par"
  Current   = "external_resistor_des.plt"
  Output    = "external_resistor_des.log"
}

Electrode {
  { Name="drain"     Voltage=0.0 Resistor=1e7 }
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
  Quasistationary(
    InitialStep=1e-7 Increment=1.41
    MinStep=1e-10 MaxStep=0.025
    Goal { Name="drain" Voltage=100000.0 }
  ) { Coupled { Poisson Electron Hole } }
}
