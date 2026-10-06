File {
 Grid="gate_msh.tdr" Parameter="models.par"
 Plot="barrier_negative_des.tdr" Current="barrier_negative.plt" Output="barrier_negative.log"
}
Electrode {
 { Name="substrate" Voltage=0 }
 { Name="gate" Voltage=0 Barrier=-0.3 }
}
Physics { Temperature=300 }
Plot { Potential eDensity hDensity eQuasiFermi hQuasiFermi ElectricField/Vector Doping }
Math { RelErrControl Digits=8 Iterations=100 NotDamped=100 Method=Blocked }
Solve {
 Coupled(Iterations=100) { Poisson }
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.01 MaxStep=0.1 MinStep=1e-8
   Goal { Name="gate" Voltage=0.094814946772048894 }) { Coupled { Poisson Electron Hole } }
 Quasistationary(InitialStep=0.025 MaxStep=0.025 MinStep=1e-8
   Goal { Name="gate" Voltage=0.13481494677204889 }) {
   Coupled { Poisson Electron Hole }
   Plot(FilePrefix="barrier_negative_scan" Time=(Range=(0 1) Intervals=40) NoOverWrite)
 }
}
