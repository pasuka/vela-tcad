* Two-dimensional n-type silicon Schottky diode translated from the Charon
* schottky_contacts/n_diode CVFEM-DD benchmark.

Electrode {
  {
    Name="anode" Voltage=0.0 Schottky Workfunction=4.75
    * Charon uses An=250 and Ap=130 A/(cm^2 K^2).  At 300 K, using the
    * corresponding silicon band-edge densities of states, these are the
    * thermionic surface velocities in cm/s.
    eRecVelocity=4.91e6 hRecVelocity=2.74e6
  }
  { Name="cathode" Voltage=0.0 }
}

File {
  Grid="schottky_n_msh.tdr"
  Plot="schottky_n_des.tdr"
  Current="schottky_n_iv"
  Output="schottky_n_des.log"
}

Physics {
  Recombination(SRH)
}

Plot {
  Potential ElectricField/Vector SpaceCharge
  eDensity hDensity eCurrent/Vector hCurrent/Vector
  eQuasiFermi hQuasiFermi
  DonorConcentration AcceptorConcentration Doping
  ConductionBandEnergy ValenceBandEnergy BandGap
  SRH eMobility hMobility
}

Math {
  Extrapolate
  RelErrControl
  Digits=6
  ErrRef(Electron)=1.0e8
  ErrRef(Hole)=1.0e8
  Iterations=25
  NotDamped=50
  NumberOfThreads=4
  ExitOnFailure
}

Solve {
  Coupled(Iterations=100 LineSearchDamping=1.0e-3) { Poisson }
  Coupled(Iterations=100) { Poisson Electron Hole }

  NewCurrentPrefix="forward_"
  Quasistationary(
    InitialStep=0.01 Increment=1.35 Decrement=2.0
    MaxStep=0.05 MinStep=1.0e-6
    Goal { Name="anode" Voltage=1.0 }
  ) {
    Coupled { Poisson Electron Hole }
  }
}
