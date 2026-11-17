


File{
  Grid      = "n1_msh.tdr"
  Plot      = "run_des.tdr"
  Parameter = "pp5_des.par"
  Current   = "run_des.plt"
  Output    = "run_des.log"
}

Electrode{
  { Name="drain"     Voltage=0.0 }
  { Name="source"    Voltage=0.0 }
  { Name="gate"      Voltage=0.0 Barrier= -0.55 }
  { Name="substrate" Voltage=0.0 }
}

Physics{

  Recombination(
    Avalanche(ElectricField)
  )
  Fermi
}





Plot{
  *--Density and Currents, etc
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  *--Fields and charges
  ElectricField/Vector Potential SpaceCharge

  *--Doping Profiles
  Doping

  *--Generation/Recombination
  SRH Band2Band * Auger
  ImpactIonization eImpactIonization hImpactIonization
  eIonIntegral hIonIntegral MeanIonIntegral eAlphaAvalanche hAlphaAvalanche
}

Math {
  Extrapolate
  Iterations=20
  Notdamped =100
  RelErrControl
  ComputeIonizationIntegrals ( )
  BreakAtIonIntegral(1 1.1)
  Transient=BE
}

Solve {
  *- Build-up of initial solution:
  Coupled(Iterations=100){ Poisson }
  Coupled { Poisson Electron Hole }

  Quasistationary(
    InitialStep=1e-4 Increment=1.41
    MinStep=1e-7 MaxStep=0.025
    Goal{ Name="drain" Voltage=100. }
  ) { Coupled { Poisson } }
}
