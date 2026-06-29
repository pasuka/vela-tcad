File {
  Grid    = "pn2d_msh.tdr"
  Plot    = "pn2d_bv_des.tdr"
  Current = "pn2d_bv.plt"
  Output  = "pn2d_bv.log"
}

Electrode {
  { Name="Anode"   Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}

Physics {
  Mobility(
    DopingDependence
    HighFieldSaturation
  )

  Recombination(
    SRH
    Avalanche(VanOverstraeten)
  )

  EffectiveIntrinsicDensity(
    OldSlotboom
  )
}

Plot {

  # Electrostatic and quasi-Fermi potentials

  Potential
  eQuasiFermi
  hQuasiFermi

  # Carrier density

  eDensity
  hDensity

  # Electric field

  ElectricField
  ElectricField/Vector

  # Current density

  eCurrent
  hCurrent
  TotalCurrent
  eCurrentDensity/Vector
  hCurrentDensity/Vector
  TotalCurrentDensity/Vector

  # Doping and charge

  Doping
  DonorConcentration
  AcceptorConcentration
  SpaceCharge

  # Recombination and avalanche source

  SRHRecombination
  eAlphaAvalanche
  hAlphaAvalanche
  AvalancheGeneration

  # Mobility, velocity, and ionization-integral diagnostics

  eMobility
  hMobility
  eVelocity
  hVelocity
  eIonIntegral
  hIonIntegral
  MeanIonIntegral
}
Math {
  Extrapolate
  RelErrControl
  Digits=5
  Iterations=80
  NotDamped=100
  Method=Blocked
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled(Iterations=100) { Poisson Electron Hole }

  Quasistationary(
    InitialStep=1e-4
    MinStep=1e-10
    MaxStep=0.05
    Increment=1.2
    Decrement=2.0
    Goal { Name="Anode" Voltage=-20.0 }
  ) {
    Coupled { Poisson Electron Hole }
    Plot(FilePrefix="pn2d_bv_multibias" Time=(Range=(0 1) Intervals=400) NoOverWrite)
  }
}
