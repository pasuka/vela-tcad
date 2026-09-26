* BVmethods NMOS physical-model ablation oracle.
* Set SRH_VARIANT to Constant or DopingDep and ENORMAL_VARIANT to Off or On
* when materializing the four A/B/C/D decks.  Constant uses an explicit
* 1e-7 s lifetime in full_physics_constant_srh.par; DopingDep uses models.par.
File {
  Grid      = "bvmethods_nmos_msh.tdr"
  Plot      = "full_physics_ablation_des.tdr"
  Parameter = "SRH_PARAMETER_FILE"
  Current   = "full_physics_ablation_des.plt"
  Output    = "full_physics_ablation_des.log"
}

Electrode {
  { Name="drain"     Voltage=0.0 }
  { Name="source"    Voltage=0.0 }
  { Name="gate"      Voltage=0.0 Barrier=-0.55 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(DopingDep HighFieldsaturation(GradQuasiFermi) ENORMAL_MODEL)
  Recombination(
    SRH(SRH_MODEL)
    Band2Band(E2)
    Avalanche(Eparallel)
  )
  Fermi
}

CurrentPlot {
  AvalancheGeneration(Integrate(Semiconductor))
}

Plot {
  eDensity hDensity eQuasiFermi hQuasiFermi
  eMobility hMobility
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  ElectricField/Vector Potential SpaceCharge Doping
  SRH Band2Band AvalancheGeneration
  eAvalancheGeneration hAvalancheGeneration
  eAlphaAvalanche hAlphaAvalanche
}

Math {
  Extrapolate
  Iterations=20
  Notdamped=100
  RelErrControl
  AvalDerivatives
  Transient=BE
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled { Poisson Electron Hole }
  Quasistationary(
    InitialStep=0.0001 Increment=1.41
    MinStep=1e-7 MaxStep=0.025
    Goal { Name="drain" Voltage=6.0 }
  ) { Coupled { Poisson Electron Hole } }
  Set("drain" mode current)
  Quasistationary(
    InitialStep=0.01 Increment=1.2
    MinStep=1e-7 MaxStep=0.1
    Goal { Name="drain" Current=1e-4 }
  ) { Coupled { Poisson Electron Hole } }
}
