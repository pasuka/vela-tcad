File {
 Grid="gate_msh.tdr" Parameter="models.par"
 Plot="workfunction_des.tdr" Current="workfunction.plt" Output="workfunction.log"
}
Electrode {
 { Name="substrate" Voltage=0 }
 { Name="gate" Voltage=0 WorkFunction=4.6110036106760788 }
}
Physics { Temperature=300 }
Plot { Potential eDensity hDensity eQuasiFermi hQuasiFermi ElectricField/Vector Doping }
Math { RelErrControl Digits=8 Iterations=100 NotDamped=100 Method=Blocked }
Solve {
 Coupled(Iterations=100) { Poisson }
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.01 MaxStep=0.1 MinStep=1e-8
   Goal { Name="gate" Voltage=0.39481494677204887 }) { Coupled { Poisson Electron Hole } }
 Quasistationary(InitialStep=0.025 MaxStep=0.025 MinStep=1e-8
   Goal { Name="gate" Voltage=0.4348149467720489 }) {
   Coupled { Poisson Electron Hole }
   Plot(FilePrefix="workfunction_scan" Time=(Range=(0 1) Intervals=40) NoOverWrite)
 }
}
