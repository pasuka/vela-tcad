File {
  Grid    = "pn2d_msh.tdr"
  Plot    = "pn2d_iv_des.tdr"
  Current = "pn2d_iv.plt"
  Output  = "pn2d_iv.log"
}

Electrode {
  { Name="Anode" Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}

Physics {

  # Low-field mobility

  Mobility(
    DopingDependence
  )

  # SRH recombination

  Recombination(
    SRH
  )

  # Intrinsic density model

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

  # Forward bias sweep

  Quasistationary(
    InitialStep=1e-3
    MinStep=1e-8
    MaxStep=0.025
    Increment=1.3
    Goal {
      Name="Anode"
      Voltage=2.0
    }
  ) {
    Coupled {
      Poisson Electron Hole
    }
    Plot(FilePrefix="pn2d_iv_multibias" Time=(Range=(0 1) Intervals=40) NoOverWrite)
  }
}
