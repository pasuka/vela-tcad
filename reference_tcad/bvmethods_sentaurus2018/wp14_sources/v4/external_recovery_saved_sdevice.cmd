* Same-model fixed-resistor DC continuation after transient branch recovery.
* Sentaurus Training BVmethods, Workbench node 5: external resistor method.
* In 2-D the resistance is interpreted as ohm*um.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "external_recovery_saved_des.tdr"
  Parameter = "models.par"
  Current   = "external_recovery_saved_des.plt"
  Output    = "external_recovery_saved_des.log"
}

Electrode {
  { Name="drain"     Voltage=14513.0674587468 Resistor=1e7 }
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
  ExtendedPrecision(80)
  Digits=8
  Iterations=80
  Notdamped=100
  RelErrControl
  AvalDerivatives
  ErrRef(Electron)=1e10
  ErrRef(Hole)=1e10
  Transient=BE
}

Solve {
  Load(FilePrefix="../external_transient_saved_v2/external_transient_seed")
  Coupled(Iterations=100) { Poisson Electron Hole }
  Quasistationary(
    InitialStep=0.1 Increment=1.2 MinStep=1e-10 MaxStep=0.1
    Goal { Name="drain" Voltage=14300.0 }
  ) { Coupled { Poisson Electron Hole } }
  Save(FilePrefix="external_recovery_saved_below_stop")
  Quasistationary(
    InitialStep=1e-7 Increment=1.41 MinStep=1e-10 MaxStep=0.025
    Goal { Name="drain" Voltage=100000.0 }
    BreakCriteria { Current(Contact="drain" AbsVal=1.443e-3) }
  ) { Coupled { Poisson Electron Hole } }
  Plot(FilePrefix="external_recovery_saved_final")
}
