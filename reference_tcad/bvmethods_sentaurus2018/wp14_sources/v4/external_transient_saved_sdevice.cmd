* Slow transient diagnostic with the original fixed resistor and physical models.
* Sentaurus Training BVmethods, Workbench node 5: external resistor method.
* In 2-D the resistance is interpreted as ohm*um.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "external_transient_saved_des.tdr"
  Parameter = "models.par"
  Current   = "external_transient_saved_des.plt"
  Output    = "external_transient_saved_des.log"
}

Electrode {
  { Name="drain"     Voltage=0.0 Voltage=(0.0 at 0.0, 100000.0 at 1.0) Resistor=1e7 }
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
  ExtendedPrecision(80)
  Iterations=80
  Notdamped=100
  RelErrControl
  AvalDerivatives
  ErrRef(Electron)=1e10
  ErrRef(Hole)=1e10
  Transient=BE
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled { Poisson Electron Hole }
  Transient(
    InitialTime=0 FinalTime=1
    BreakCriteria { Current(Contact="drain" AbsVal=1.443e-3) }
    InitialStep=1e-7 Increment=1.41
    MinStep=1e-14 MaxStep=0.025

  ) { Coupled { Poisson Electron Hole } }
  Save(FilePrefix="external_transient_seed")
}
