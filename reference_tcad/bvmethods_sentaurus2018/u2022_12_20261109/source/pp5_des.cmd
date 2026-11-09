


File{
  Grid      = "n1_msh.tdr"
  Plot      = "n5_des.tdr"
  Parameter = "pp5_des.par"
  Current   = "n5_des.plt"
  Output    = "n5_des.log"
}

Electrode{
  { Name="drain"     Voltage=0.0 Resistor= 1e7  }
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
    InitialStep=1e-07 Increment=1.41 
    MinStep=1e-10 MaxStep=0.025
    Goal{ Name="drain" Voltage=100000. }
  ) { Coupled { Poisson Electron Hole } }
}

