#setdep @previous@

File{
  Grid      = "@tdr@"
  Plot      = "@tdrdat@"
  Parameter = "@parameter@"
  Current   = "@plot@"
  Output    = "@log@"
}

Electrode{
  { Name="source"    Voltage=0.0 }
  { Name="drain"     Voltage=0.0 }
  { Name="gate"      Voltage=-0.5 }
  { Name="substrate" Voltage=0.0 }
}

*Example of Thermode section needed for self-heating simulations
#Thermode {
#  { Name="[contact-name];" Temperature=[K]; [other options] }
#  { [other contacts] }
#  ...
#}

Physics{
  eQuantumPotential
  EffectiveIntrinsicDensity( OldSlotboom )     
  Mobility(
    DopingDep
    eHighFieldsaturation( GradQuasiFermi )
    hHighFieldsaturation( GradQuasiFermi )
    Enormal
  )
  Recombination(
    SRH( DopingDep TempDependence )
  )           
}

*Examples of restricted Physics sections
#Physics(Material="[material name]"){
#  [models active only in the given material]
#}
#Physics(Region="[region name]"){
#  [models active only in the given region]
#}
#Physics(RegionInterface="[region name]/[region name]") {
#  [models active only in the given region interface]
#}
#Physics(MaterialInterface="[material name]/[material name]") {
#  [models active only in the given material interface]
#}

Plot{
  *--Density and Currents, etc
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  eMobility/Element hMobility/Element
  eVelocity hVelocity
  eQuasiFermi hQuasiFermi
  
  *--Temperature 
  eTemperature hTemperature Temperature
  
  *--Fields and charges
  ElectricField/Vector Potential SpaceCharge
  
  *--Doping Profiles
  Doping DonorConcentration AcceptorConcentration
  
  *--Generation/Recombination
  SRH Band2Band Auger
  AvalancheGeneration eAvalancheGeneration hAvalancheGeneration
  
  *--Driving forces
  eGradQuasiFermi/Vector hGradQuasiFermi/Vector
  eEparallel hEparallel eENormal hENormal
  
  *--Band structure/Composition
  BandGap 
  BandGapNarrowing
  Affinity
  ConductionBand ValenceBand
  eQuantumPotential hQuantumPotential
}

Math {
  Extrapolate               * (not needed for transient sweeps)
  #Avalderivatives          * (only if Avalanche models are active)
  RelErrControl
  Digits=5                  * (default)
  ErrRef(electron)=1.e10    * (default)
  ErrRef(hole)=1.e10        * (default)
  Iterations=20
  Notdamped=100
}

CurrentPlot {
  ElectricField/Vector( (0.01 0.025) )
}

Solve {
  *- Build-up of initial solution:
  NewCurrentPrefix="init_"
  Coupled(Iterations=100){ Poisson eQuantumPotential }
  Coupled{ Poisson Electron Hole eQuantumPotential }
  
  Save ( FilePrefix = "n@node@_init" )
  
  *- Bias drain to target bias
  Quasistationary(
    InitialStep=0.01 MinStep=1e-5 MaxStep=0.1
    Goal{ Name="drain" Voltage= 0.1  }
  ) { Coupled { Poisson Electron Hole eQuantumPotential } }
  
  NewCurrentPrefix="IdVgsLin_"
  *- Gate voltage sweep
  Quasistationary(
    InitialStep=1e-3 Increment=1.3 Decrement=2 MaxStep=0.1 MinStep=1e-5
    Goal{ Name="gate" Voltage= 2.2 }
  ) { Coupled { Poisson Electron Hole eQuantumPotential }
    *Plot at a set of given values of the t variable
    Plot( -Loadable Fileprefix="n@node@_inter" NoOverWrite Time= (0.3; 0.6) )
    *Plot at regular intervals
    #Plot(Fileprefix="n@node@_inter" NoOverWrite Time=(Range=(0 1) Intervals=6))
    *I-V calculated at given values of the t variable
    #CurrentPlot(Time=(0.1; 0.5; 0.8))
    *I-V calculated at regular intervals
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  } *end gate voltage sweep
  
  Plot ( FilePrefix = "n@node@_Lin" )
  
  Load ( FilePrefix = "n@node@_init" )
  
  NewCurrentPrefix=""
  *- Bias drain to target bias
  Quasistationary(
    InitialStep=0.01 MinStep=1e-5 MaxStep=0.1
    Goal{ Name="drain" Voltage= 1.1  }
  ) { Coupled { Poisson Electron Hole eQuantumPotential } }
  
  NewCurrentPrefix="IdVgsSat_"
  *- Gate voltage sweep
  Quasistationary(
    InitialStep=1e-3 Increment=1.3 Decrement=2 MaxStep=0.1 MinStep=1e-5
    Goal{ Name="gate" Voltage= 2.2 }
  ) { Coupled { Poisson Electron Hole eQuantumPotential }
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  } *end gate voltage sweep   
  
  Plot ( FilePrefix = "n@node@_Sat" )
  
  *- Example of parameter (lattice temperature) sweep
  #Quasistationary(
  #InitialStep=1.e-3 Increment=1.5
  #MaxStep=0.05 MinStep=1.e-7
  #Goal { Model="DeviceTemperature" Parameter="Temperature" Value=400.} 
  #) { Coupled { Poisson Electron Hole eQuantumPotential } }
  
}*end Solve
