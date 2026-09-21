File {
 Grid="pn2d_msh.tdr" Parameter="models.par"
 Plot="reverse_des.tdr" Current="reverse.plt" Output="reverse.log"
}
Electrode { { Name="Anode" Voltage=0 } { Name="Cathode" Voltage=0 } }
Physics { Temperature=300 AreaFactor=1
 Mobility(DopingDependence)
 # Recombination disabled
 EffectiveIntrinsicDensity(OldSlotboom)
}
Plot { Potential ElectricField/Vector eDensity hDensity eQuasiFermi hQuasiFermi
 eCurrent/Vector hCurrent/Vector TotalCurrent/Vector Doping DonorConcentration AcceptorConcentration
 SRHRecombination eMobility hMobility EffectiveIntrinsicDensity BandGap BandGapNarrowing
}
Math { Extrapolate RelErrControl Digits=8 Iterations=80 NotDamped=100 }
Solve {
 Coupled(Iterations=100) { Poisson }
 Coupled(Iterations=100) { Poisson Electron Hole }
 Plot(FilePrefix="reverse_p0p00")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.1 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p10")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.2 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p20")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.3 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p30")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.4 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p40")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.5 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p50")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.6 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p60")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.7 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p70")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.8 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p80")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-0.9 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m0p90")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=-1 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="reverse_m1p00")
}
