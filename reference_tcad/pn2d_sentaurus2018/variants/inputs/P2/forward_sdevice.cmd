File {
 Grid="pn2d_msh.tdr" Parameter="models.par"
 Plot="forward_des.tdr" Current="forward.plt" Output="forward.log"
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
 Plot(FilePrefix="forward_p0p00")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.02 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p02")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.04 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p04")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.06 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p06")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.08 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p08")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.1 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p10")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.12 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p12")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.14 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p14")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.16 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p16")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.18 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p18")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.2 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p20")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.22 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p22")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.24 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p24")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.26 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p26")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.28 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p28")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.3 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p30")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.32 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p32")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.34 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p34")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.36 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p36")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.38 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p38")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.4 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p40")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.42 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p42")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.44 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p44")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.46 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p46")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.48 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p48")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.5 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p50")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.52 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p52")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.54 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p54")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.56 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p56")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.58 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p58")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.6 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p60")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.62 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p62")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.64 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p64")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.66 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p66")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.68 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p68")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.7 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p70")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.72 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p72")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.74 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p74")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.76 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p76")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.78 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p78")
 Quasistationary(InitialStep=0.05 MinStep=1e-8 MaxStep=1 Increment=1.3 Goal { Name="Anode" Voltage=0.8 }) {
  Coupled { Poisson Electron Hole }
 }
 Plot(FilePrefix="forward_p0p80")
}
