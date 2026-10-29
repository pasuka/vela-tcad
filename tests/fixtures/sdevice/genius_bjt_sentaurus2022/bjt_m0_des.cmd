File {
  Grid    = "bjt_msh.tdr"
  Plot    = "bjt_m0_des.tdr"
  Current = "bjt_m0"
  Output  = "bjt_m0"
}

Electrode {
  { Name="emitter"   Voltage=0.0 }
  { Name="base"      Voltage=0.0 }
  { Name="collector" Voltage=0.0 }
}

Physics {
  Temperature=300
  Recombination(SRH)
}

Plot {
  Potential ElectricField/Vector SpaceCharge
  eDensity hDensity eCurrent/Vector hCurrent/Vector
  eQuasiFermi hQuasiFermi Doping DonorConcentration AcceptorConcentration
  eMobility hMobility SRHRecombination
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
  Plot(FilePrefix="bjt_m0_equilibrium")

  NewCurrentPrefix="bjt_m0_base_"
  Quasistationary(
    InitialStep=1e-3 Increment=1.35 Decrement=2.0
    MinStep=1e-7 MaxStep=0.04
    Goal { Name="base" Voltage=0.70 }
  ) { Coupled { Poisson Electron Hole } }
  Plot(FilePrefix="bjt_m0_vbe070")

  NewCurrentPrefix="bjt_m0_collector_"
  Quasistationary(
    InitialStep=0.005 Increment=1.35 Decrement=2.0
    MinStep=1e-7 MaxStep=0.05
    Goal { Name="collector" Voltage=3.0 }
  ) {
    Coupled { Poisson Electron Hole }
    CurrentPlot(Time=(Range=(0 1) Intervals=30))
    Plot(FilePrefix="bjt_m0_vce" Time=(0; 0.0333333333333; 0.333333333333; 1) NoOverWrite)
  }
}
