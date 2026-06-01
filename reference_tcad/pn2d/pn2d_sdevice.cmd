File {
  Grid      = "pn2d_msh.tdr"
  Plot      = "pn2d_des.tdr"
  Current   = "pn2d_iv.plt"
  Output    = "pn2d_des.log"
}

Electrode {
  { Name="Anode"   Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}

Physics {
  Fermi
  EffectiveIntrinsicDensity( OldSlotboom )

  Mobility(
    DopingDep
    HighFieldSaturation
  )

  Recombination(
    SRH
    Auger
  )
}

Plot {

  Doping
  DonorConcentration
  AcceptorConcentration

  eDensity
  hDensity
  Potential
  ElectricField
  eCurrentDensity
  hCurrentDensity
  TotalCurrentDensity

  eQuasiFermi
  hQuasiFermi
  SpaceCharge
  SRHRecombination
  AugerRecombination
}

Math {
  Extrapolate
  RelErrControl
  CNormPrint
  Iterations=30
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled(Iterations=100) { Poisson Electron Hole }

  Quasistationary(
    InitialStep=1e-3
    MinStep=1e-6
    MaxStep=0.02
    Increment=1.4
    Decrement=2.0
    Goal { Name="Anode" Voltage=1.0 }
  ) {
    Coupled { Poisson Electron Hole }
    Plot(FilePrefix="pn2d_forward" Time=(0;0.25;0.5;0.75;1.0))
  }
}
