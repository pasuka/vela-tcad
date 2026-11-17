


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

  EffectiveIntrinsicDensity( OldSlotboom )
  Mobility(
    DopingDep
    HighFieldsaturation( GradQuasiFermi )
    Enormal
  )
  Recombination(
    SRH( DopingDep )

    Avalanche( GradQuasiFermi )
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
  AvalDerivatives
  ErrRef(Electron)=1.e10
  ErrRef(Hole)=1.e10
  BreakCriteria{ Current(Contact="drain" AbsVal=1.3e-3) }
  Transient=BE
}

Solve {
  *- Build-up of initial solution:
  Coupled(Iterations=100){ Poisson }
  Coupled { Poisson Electron Hole }

  Quasistationary(
    InitialStep=0.002 Increment=1.41
    MinStep=2e-6 MaxStep=0.025
    Goal{ Name="drain" Voltage=6. }
  ) { Coupled { Poisson Electron Hole } }

  Set ("drain" mode current)
  QuasiStationary (
    InitialStep=1e-6 Maxstep=0.05 MinStep=1e-9 Increment=1.41
    Goal { name=drain current=1.3e-3 }
  ) { Coupled { Poisson Electron Hole } }
}
