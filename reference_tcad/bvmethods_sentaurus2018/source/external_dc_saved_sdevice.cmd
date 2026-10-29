* Sentaurus Training BVmethods, Workbench node 5: external resistor method.
* In 2-D the resistance is interpreted as ohm*um.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "external_dc_saved_des.tdr"
  Parameter = "models.par"
  Current   = "external_dc_saved_des.plt"
  Output    = "external_dc_saved_des.log"
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
  Plot(FilePrefix="external_dc_saved_steady")
  Save(FilePrefix="external_dc_saved_steady")
}
