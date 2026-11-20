* Linear and saturation branches of the SingleDevice SDevice deck.
File {
  Grid = "device.grd"
  Parameter = "models.par"
  Plot = "singledevice_des.tdr"
  Current = "singledevice.plt"
  Output = "singledevice.log"
}
Electrode {
  { Name="source" Voltage=0.0 }
  { Name="drain" Voltage=0.0 }
  { Name="gate" Voltage=-0.5 }
  { Name="substrate" Voltage=0.0 }
}
Physics {
  Temperature = 300
  eQuantumPotential
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(DopingDep eHighFieldsaturation(GradQuasiFermi)
                    hHighFieldsaturation(GradQuasiFermi) Enormal)
  Recombination(SRH(DopingDep TempDependence))
}
Plot {
  Potential
}
Math {
  Extrapolate
  RelErrControl
  Digits=5
  ErrRef(electron)=1e10
  ErrRef(hole)=1e10
  Iterations=20
  Notdamped=100
}
Solve {
  NewCurrentPrefix="init_"
  Coupled(Iterations=100) { Poisson eQuantumPotential }
  Coupled { Poisson Electron Hole eQuantumPotential }
  Save(FilePrefix="singledevice_init")
  Quasistationary(
    InitialStep=0.01 MinStep=1e-5 MaxStep=0.1
    Goal { Name="drain" Voltage=0.1 }
  ) { Coupled { Poisson Electron Hole eQuantumPotential } }
  NewCurrentPrefix="IdVgsLin_"
  Quasistationary(
    InitialStep=1e-3 Increment=1.3 Decrement=2 MinStep=1e-5 MaxStep=0.1
    Goal { Name="gate" Voltage=2.2 }
  ) {
    Coupled { Poisson Electron Hole eQuantumPotential }
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  }
  NewCurrentPrefix=""
  Load(FilePrefix="singledevice_init")
  Quasistationary(
    InitialStep=0.01 MinStep=1e-5 MaxStep=0.1
    Goal { Name="drain" Voltage=1.1 }
  ) { Coupled { Poisson Electron Hole eQuantumPotential } }
  NewCurrentPrefix="IdVgsSat_"
  Quasistationary(
    InitialStep=1e-3 Increment=1.3 Decrement=2 MinStep=1e-5 MaxStep=0.1
    Goal { Name="gate" Voltage=2.2 }
  ) {
    Coupled { Poisson Electron Hole eQuantumPotential }
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  }
}
