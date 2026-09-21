Electrode {
	{ Name= "drain"     Voltage= 0.0 hRecVelocity=  1.93E6 }
	{ Name= "gate"      Voltage= 0.0 Material=  "PolySi"(N) }
	{ Name= "source"    Voltage= 0.0 hRecVelocity=  1.93E6 }
	{ Name= "substrate" Voltage= 0.0}
}

Thermode {
	{ Name= "th_lat"  Temperature= 300 SurfaceResistance= 0.005}
}

File {
	Grid=       "n1_fps.tdr"
	Parameters= "sdevice.par"
	Output=     "n4_des.log"
	Current=    "n4_des.plt"
	Plot=       "n4_des.tdr"
}

Physics {Fermi}

Physics (Material= "Silicon") {
	eQuantumPotential(density) hQuantumPotential(density)
	Mobility(
		HighFieldSaturation 
		Enormal (IALMob(AutoOrientation))
	)
	EffectiveIntrinsicDensity(OldSlotboom)
	Recombination(
	SRH( DopingDependence TempDependence )
	Auger
	)
}


Plot {
	eDensity hDensity
	TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
	eMobility hMobility
	eVelocity hVelocity
	eQuasiFermi hQuasiFermi
	eTemperature hTemperature Temperature 
	ElectricField/Vector Potential SpaceCharge
	Doping DonorConcentration AcceptorConcentration
	SRH Band2Band 
	AvalancheGeneration eAvalancheGeneration hAvalancheGeneration
	eGradQuasiFermi/Vector hGradQuasiFermi/Vector
	eEparallel hEparallel eENormal hENormal
	BandGap 
	BandGapNarrowing
	Affinity
	ConductionBand ValenceBand
	eBarrierTunneling hBarrierTunneling * BarrierTunneling
	eTrappedCharge  hTrappedCharge
	eGapStatesRecombination hGapStatesRecombination
	eDirectTunnel hDirectTunnel
}

Math {  Digits=6
	Extrapolate
	Notdamped= 100
	Iterations= 25
	ExitOnFailure
	ErrRef(Electron)=  1e8
	ErrRef(Hole)    =  1e8
	
	*CNormPrint
	RefDens_eGradQuasiFermi_EparallelToInterface= 1e12
	RefDens_hGradQuasiFermi_EparallelToInterface= 1e12
}

Solve {
*- Creating initial guess:
	Coupled(Iterations= 100 LineSearchDamping= 1e-4){ Poisson } 
	Coupled { Poisson Electron Hole }

	Quasistationary (
		Initialstep= 0.01 Increment= 1.35
		MaxStep= 0.4 Minstep= 1.e-4
		Goal { Name= "gate" Voltage= 4.0}
	){ Coupled { Poisson } }
	Save(FilePrefix= "n4_Vg1")

	Quasistationary (
		Initialstep= 0.01 Increment= 1.35
		MaxStep= 0.4 Minstep= 1.e-4
		Goal { Name= "gate" Voltage= 8.0}
	){ Coupled { Poisson } }
	Save(FilePrefix= "n4_Vg2")

	Load(FilePrefix= "n4_Vg1")
	Coupled { Poisson Electron Hole Temperature }
	NewCurrentPrefix= "IdVd_Vg1_"
	Quasistationary (
		Initialstep= 2.5e-4 Increment= 1.35
		Minstep= 1.e-6 MaxStep= 0.05 
		Goal { Name= "drain" Voltage= 40.0 }
	){ Coupled { Poisson Electron Hole Temperature }
		CurrentPlot( Time= (Range= (0 1) Intervals= 30) )
        Plot(FilePrefix="field_vg4" Time=(Range=(0 1) Intervals=30) NoOverWrite)
	}

	Plot(FilePrefix= "n4_Vg1")
	NewCurrentPrefix= "IdVd_Vg2_"
	Load(FilePrefix= "n4_Vg2")
	Coupled { Poisson Electron Hole Temperature }
	NewCurrentPrefix= "IdVd_Vg2_"
	Quasistationary (
		Initialstep= 2.5e-4 Increment= 1.35
		Minstep= 1.e-6 MaxStep= 0.05
		Goal { Name= "drain" Voltage= 40.0 }
	){ Coupled { Poisson Electron Hole Temperature }
		CurrentPlot( Time= (Range= (0 1) Intervals= 30)  )
        Plot(FilePrefix="field_vg8" Time=(Range=(0 1) Intervals=30) NoOverWrite)
	}
}

