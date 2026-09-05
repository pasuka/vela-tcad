
*-- DD






File{
   Grid      = "n1_msh.tdr"
   Plot      = "n13_des.tdr"
   Parameter = "pp13_des.par"
   Current   = "n13_des.plt"
   Output    = "n13_des.log"
}

Electrode{
   { Name="source"    Voltage=0.0 }
   { Name="drain"     Voltage=0.0 }
   { Name="gate"      Voltage=1. }
   { Name="substrate" Voltage=0.0 }
}

Thermode{
  { Name="substrate" Temperature=300 SurfaceResistance=5e-4 }
  { Name="drain" Temperature=300 SurfaceResistance=1e-3 }
  { Name="source" Temperature=300 SurfaceResistance=1e-3 }
}


Physics{

   * DriftDiffusion
   eQuantumPotential
   Fermi
   EffectiveIntrinsicDensity( OldSlotboom )
   Mobility(
      DopingDep
      eHighFieldsaturation( 		 GradQuasiFermi )
      hHighFieldsaturation( GradQuasiFermi )
      Enormal
   )
   Recombination(
      SRH( DopingDep TempDependence )
   )
}

Plot{
*--Density and Currents, etc
   eDensity hDensity
   TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
   eMobility hMobility
   eVelocity hVelocity
   eQuasiFermi hQuasiFermi

*--Temperature
   eTemperature Temperature * hTemperature

*--Fields and charges
   ElectricField/Vector Potential SpaceCharge

*--Doping Profiles
   Doping DonorConcentration AcceptorConcentration

*--Generation/Recombination
   SRH Band2Band * Auger
   ImpactIonization eImpactIonization hImpactIonization

*--Driving forces
   eGradQuasiFermi/Vector hGradQuasiFermi/Vector
   eEparallel hEparallel eENormal hENormal

*--Band structure/Composition
   BandGap
   BandGapNarrowing
   Affinity
   ConductionBand ValenceBand
   eQuantumPotential
}

Math {
   Extrapolate
   Iterations= 20
   Notdamped= 100
   Method= Blocked
   SubMethod= Pardiso
}

Solve {
   *- Build-up of initial solution:
   NewCurrentPrefix="init_"
   Coupled(Iterations=100){ Poisson eQuantumPotential }
   Coupled{ Poisson Electron Hole eQuantumPotential }

   *-  drain bias sweep
   NewCurrentPrefix="IdVds_"
   Quasistationary(
      InitialStep=1e-3 MinStep=1e-5 MaxStep=1
      Goal{ Name="drain" Voltage= 2. }
   ) { Coupled { Poisson Electron Hole eQuantumPotential }
       CurrentPlot(Time=(Range=(0 1) Intervals=20))
     }

   System("rm init_n13_des.plt")
}
