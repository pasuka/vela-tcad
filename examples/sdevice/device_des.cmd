File {
  Grid="mesh.grd" Parameter="models.par"
  Output="run" Current="current" Plot="final"
}
Electrode {
  {Name="anode" Voltage=0}
  {Name="cathode" Voltage=0}
}
Math { NumberOfThreads=1 ExitOnFailure }
Solve {
  Coupled {Poisson}
  Coupled {Poisson Electron Hole}
  Quasistationary(InitialStep=0.1 MinStep=0.001 MaxStep=0.25 Increment=2 Decrement=2
    Goal {Name="anode" Voltage=1}) {
    Coupled {Poisson Electron Hole}
  }
}
