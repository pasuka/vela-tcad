* Sentaurus Training BVmethods, Workbench node 8: transient resistor ramp.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "transient_des.tdr"
  Parameter = "models.par"
  Current   = "transient_des.plt"
  Output    = "transient_des.log"
}

Electrode {
  { Name="drain" Voltage=0.0 Voltage=(0.0 at 0.0, 100000.0 at 1.0) Resistor=1e7 }
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

CurrentPlot {
  ImpactIonization(Maximum(Semiconductor Coordinates))
  ElectricField(Maximum(Semiconductor Coordinates))
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
  Transient(
    InitialTime=0 FinalTime=1
    InitialStep=1e-6 Increment=1.41
    MinStep=1e-10 MaxStep=0.01
  ) { Coupled { Poisson Electron Hole } }
}
