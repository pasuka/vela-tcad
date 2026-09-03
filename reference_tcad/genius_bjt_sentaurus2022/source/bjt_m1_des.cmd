File {
  Grid    = "bjt_msh.tdr"
  Plot    = "bjt_m1_des.tdr"
  Current = "bjt_m1"
  Output  = "bjt_m1"
}

Electrode {
  { Name="emitter"   Voltage=0.0 }
  { Name="base"      Voltage=0.0 }
  { Name="collector" Voltage=0.0 }
}

Physics {
  Temperature=300
  Fermi
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(DopingDependence)
  Recombination(SRH(DopingDependence) Auger)
}

Plot {
  Potential ElectricField/Vector SpaceCharge
  eDensity hDensity eCurrent/Vector hCurrent/Vector
  eQuasiFermi hQuasiFermi Doping DonorConcentration AcceptorConcentration
  eMobility hMobility SRHRecombination AugerRecombination
  EffectiveBandGap EffectiveIntrinsicDensity
  ConductionBand ValenceBand BandGap
}

Math {
  Extrapolate
  Derivatives
  RelErrControl
  Digits=5
  ErrRef(Electron)=1e7
  ErrRef(Hole)=1e7
  Iterations=30
  NotDamped=100
}

Solve {
  Poisson
  Coupled(Iterations=100) { Poisson Electron Hole }
  Plot(FilePrefix="bjt_m1_equilibrium")

  NewCurrentPrefix="bjt_m1_base_"
  Quasistationary(
    InitialStep=1e-3 Increment=1.35 Decrement=2.0
    MinStep=1e-7 MaxStep=0.04
    Goal { Name="base" Voltage=0.70 }
  ) { Coupled { Poisson Electron Hole } }
  Plot(FilePrefix="bjt_m1_vbe070")

  NewCurrentPrefix="bjt_m1_collector_"
  Quasistationary(
    InitialStep=0.005 Increment=1.35 Decrement=2.0
    MinStep=1e-7 MaxStep=0.05
    Goal { Name="collector" Voltage=3.0 }
  ) {
    Coupled { Poisson Electron Hole }
    CurrentPlot(Time=(Range=(0 1) Intervals=30))
    Plot(FilePrefix="bjt_m1_vce" Time=(0; 0.333333333333; 0.666666666667; 1) NoOverWrite)
  }
}
