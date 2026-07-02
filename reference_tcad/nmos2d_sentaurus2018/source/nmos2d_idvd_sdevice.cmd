File {
  Grid    = "nmos2d_msh.tdr"
  Plot    = "nmos2d_idvd_des.tdr"
  Current = "nmos2d_idvd.plt"
  Output  = "nmos2d_idvd.log"
}

Electrode {
  { Name="Source" Voltage=0.0 }
  { Name="Drain"  Voltage=0.0 }
  { Name="Gate" Voltage=2.0 }
  { Name="Body"   Voltage=0.0 }
}

Physics {

  # Low-field mobility, matching the first pn2d IV reference model family.

  Mobility(
    DopingDependence
  )

  # SRH recombination.

  Recombination(
    SRH
  )

  # Intrinsic density model.

  EffectiveIntrinsicDensity(
    OldSlotboom
  )
}

Plot {

  # Electrostatic

  Potential
  ElectricField

  # Carrier density

  eDensity
  hDensity

  # Quasi Fermi

  eQuasiFermi
  hQuasiFermi

  # Current

  eCurrent
  hCurrent
  TotalCurrent

  # Doping

  Doping
  DonorConcentration
  AcceptorConcentration

  # Recombination

  SRHRecombination

  # Mobility diagnostics

  eMobility
  hMobility
}

Math {
  Extrapolate
  RelErrControl
  Digits=5
  Iterations=50
  NotDamped=100
}

Solve {

  Coupled(Iterations=100) { Poisson }

  Coupled(Iterations=100) {
    Poisson Electron Hole
  }

  # Id-Vd output sweep at fixed on-state gate bias.

  Quasistationary(
    InitialStep=1e-3
    MinStep=1e-8
    MaxStep=0.02
    Increment=1.3
    Goal {
      Name="Drain"
      Voltage=0.5
    }
  ) {
    Coupled {
      Poisson Electron Hole
    }
    Plot(FilePrefix="nmos2d_idvd_multibias" Time=(Range=(0 1) Intervals=25) NoOverWrite)
  }
}