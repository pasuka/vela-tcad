* Isothermal DG adaptation of ../source/pp13_des.cmd for the DF-ISE input.
File {
  Grid = "n1_msh.grd"
  Parameter = "pp13_des.par"
  Plot = "n13_des.tdr"
  Current = "n13_des.plt"
  Output = "n13_des.log"
}
Electrode {
  { Name="source" Voltage=0.0 }
  { Name="drain" Voltage=0.0 }
  { Name="gate" Voltage=1.0 }
  { Name="substrate" Voltage=0.0 }
}
Physics {
  Temperature = 300
  eQuantumPotential
  Fermi
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
  Iterations=20
  Notdamped=100
  Method=Blocked
  SubMethod=Pardiso
}
Solve {
  NewCurrentPrefix="init_"
  Coupled(Iterations=100) { Poisson eQuantumPotential }
  Coupled { Poisson Electron Hole eQuantumPotential }
  NewCurrentPrefix="IdVds_"
  Quasistationary(
    InitialStep=1e-3 MinStep=1e-5 MaxStep=1
    Goal { Name="drain" Voltage=2. }
  ) {
    Coupled { Poisson Electron Hole eQuantumPotential }
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  }
}
