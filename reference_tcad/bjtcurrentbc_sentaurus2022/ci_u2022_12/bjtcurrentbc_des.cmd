Electrode{
  { name="base"      voltage=0 }
  { name="emitter"   voltage=0 }
  { name="collector" voltage=0 }
}

File{
  Grid    = "bjt.grd"
  Parameter = "models.par"
  Plot    = "native.tdr"
  Current = "native.plt"
  Output  = "native.log"
}

Plot { Potential }

Physics{
  Recombination(
    SRH(DopingDep)
    Auger
  )
  Mobility( DopingDep HighFieldSaturation)
  EffectiveIntrinsicDensity( OldSlotboom )
}

Math{
  Extrapolate
  RelErrControl
  Digits=5
  ErrReff(electron)= 1.0e7
  ErrReff(hole)    = 1.0e7
  Iterations=20
  Notdamped=100
}

Solve {

  Coupled(Iterations=100){ Poisson }
  Coupled{ Poisson Electron Hole }

  Quasistationary (
    InitialStep=0.01 Increment=1.3
    MaxStep =0.2 MinStep = 1e-6
    Goal{ Name="base" Voltage=0.5 }
  ){ Coupled {Poisson Electron Hole } }

  Set( "base" mode current )

  Quasistationary (
    InitialStep=0.01 Increment=1.3
    MaxStep =0.2 MinStep = 1e-6
    Goal{ Name="base" Current=1e-8 }
  ){ Coupled {Poisson Electron Hole } }


  NewCurrentPrefix="single_"
  Quasistationary (
    InitialStep=0.01 Increment=1.3
    MaxStep=0.2 MinStep=1e-6
    Goal { Name="collector" Voltage=1.0 }
  ){
    Coupled { Poisson Electron Hole }
    CurrentPlot(Time=(Range=(0 1) Intervals=20))
  }
}
