File {
  Grid    = "pn2d_msh.tdr"
  Plot    = "pn2d_0v_des.tdr"
  Current = "pn2d_0v.plt"
  Output  = "pn2d_0v.log"
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
}

Math {
  Extrapolate
  RelErrControl
  Digits=5
  Iterations=50
  NotDamped=100
}

Solve {

  # Poisson initialization

  Coupled(Iterations=100) { Poisson }

  #Thermal equilibrium

  Coupled(Iterations=100) {
    Poisson Electron Hole
  }
}
