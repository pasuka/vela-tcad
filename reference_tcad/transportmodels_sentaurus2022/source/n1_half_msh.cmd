Title ""

Controls {
}

IOControls {
	EnableSections
}

Definitions {
	Constant "Const.Substrate" {
		Species = "BoronActiveConcentration"
		Value = 1.5e+18
	}
	Constant "Const.Gate" {
		Species = "ArsenicActiveConcentration"
		Value = 1e+20
	}
	AnalyticalProfile "Impl.SDprof" {
		Species = "ArsenicActiveConcentration"
		Function = Gauss(PeakPos = 0, PeakVal = 6e+20, ValueAtDepth = 1.5e+18, Depth = 0.05)
		LateralFunction = Gauss(Factor = 0.4)
	}
	AnalyticalProfile "Impl.SDextprof" {
		Species = "ArsenicActiveConcentration"
		Function = Gauss(PeakPos = 0, PeakVal = 2e+20, ValueAtDepth = 1.5e+18, Depth = 0.0125)
		LateralFunction = Gauss(Factor = 0.25)
	}
	Refinement "Global_def" {
		MaxElementSize = ( 0.125 0.0925 )
		MinElementSize = ( 0.025 0.0185 )
	}
	Refinement "Active_def" {
		MaxElementSize = ( 0.00625 0.025 )
		MinElementSize = ( 0.00125 0.005 )
	}
	Refinement "Channel_def" {
		MaxElementSize = ( 0.003125 0.00625 )
		MinElementSize = ( 0.000625 0.00125 )
	}
	Refinement "RD_def_0" {
		MaxElementSize = ( 1.1012 0.185 0 )
		MinElementSize = ( 0.001 0.001 0.001 )
		RefineFunction = MaxTransDiff(Variable = "DopingConcentration",Value = 1)
	}
	Refinement "refintdef_0" {
		MaxElementSize = 100
		MinElementSize = 0.0002
		RefineFunction = MaxLenInt(Interface("R.Polygate","R.Gateox"), Value=0.0002, factor=2, UseRegionNames)
	}
	Refinement "refintdef_1" {
		MaxElementSize = 100
		MinElementSize = 0.0001
		RefineFunction = MaxLenInt(Interface("R.Substrate","R.Gateox"), Value=0.0001, factor=1.6, UseRegionNames)
	}
}

Placements {
	Constant "PlaceCD.Substrate" {
		Reference = "Const.Substrate"
		EvaluateWindow {
			Element = region ["R.Substrate"]
		}
	}
	Constant "PlaceCD.Gate" {
		Reference = "Const.Gate"
		EvaluateWindow {
			Element = region ["R.Polygate"]
		}
	}
	AnalyticalProfile "Impl.Drain" {
		Reference = "Impl.SDprof"
		ReferenceElement {
			Element = Line [(0 0.37) (0 0.085)]
			Direction = positive
		}
	}
	AnalyticalProfile "Impl.DrainExt" {
		Reference = "Impl.SDextprof"
		ReferenceElement {
			Element = Line [(0 0.37) (0 0.025)]
			Direction = positive
		}
	}
	Refinement "Global_plc" {
		Reference = "Global_def"
		RefineWindow = Rectangle [(-100 -100) (100 100)]
	}
	Refinement "Active_plc" {
		Reference = "Active_def"
		RefineWindow = Rectangle [(0 0) (0.06 0.185)]
	}
	Refinement "Channel_plc" {
		Reference = "Channel_def"
		RefineWindow = Rectangle [(0 0) (0.02 0.045)]
	}
	Refinement "RD_plc_0" {
		Reference = "RD_def_0"
		RefineWindow = material ["Silicon"]
	}
	Refinement "refintplc_0" {
		Reference = "refintdef_0"
		RefineWindow = region ["R.Polygate"]
	}
	Refinement "refintplc_1" {
		Reference = "refintdef_1"
		RefineWindow = region ["R.Substrate"]
	}
}
