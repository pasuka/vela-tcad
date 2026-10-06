File {
 Grid="gate_msh.tdr" Parameter="models.par"
 Plot="barrier_positive_des.tdr" Current="barrier_positive.plt" Output="barrier_positive.log"
}
Electrode {
 { Name="substrate" Voltage=0 }
 { Name="gate" Voltage=0 Barrier=0.55 }
}
Physics { Temperature=300 }
Plot { Potential eDensity hDensity eQuasiFermi hQuasiFermi ElectricField/Vector Doping }
Math { RelErrControl Digits=8 Iterations=100 NotDamped=100 Method=Blocked }
Solve {
 Coupled(Iterations=100) { Poisson }
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.01 MaxStep=0.1 MinStep=1e-8
   Goal { Name="gate" Voltage=0.94481494677204891 }) { Coupled { Poisson Electron Hole } }
 Quasistationary(InitialStep=0.025 MaxStep=0.025 MinStep=1e-8
   Goal { Name="gate" Voltage=0.98481494677204895 }) {
   Coupled { Poisson Electron Hole }
   Plot(FilePrefix="barrier_positive_scan" Time=(Range=(0 1) Intervals=40) NoOverWrite)
 }
}
