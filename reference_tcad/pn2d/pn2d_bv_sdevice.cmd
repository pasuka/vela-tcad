File {
  Grid      = "pn2d_msh.tdr"
  Plot      = "pn2d_bv_des.tdr"
  Current   = "pn2d_bv.plt"
  Output    = "pn2d_bv_des.log"
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
    Avalanche( OkutoCrowell )
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
  AvalancheGeneration
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
    MinStep=1e-8
    MaxStep=0.2
    Increment=1.2
    Decrement=2.0
    Goal { Name="Cathode" Voltage=50.0 }
  ) {
    Coupled { Poisson Electron Hole }
    Plot(FilePrefix="pn2d_reverse" Time=(0;0.2;0.4;0.6;0.8;1.0))
  }
}
