Title "Untitled"

Controls {
}

IOControls {
	EnableSections
}

Definitions {
	Constant "P_Doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+17
	}
	Constant "N_Doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+17
	}
	Refinement "GlobalMesh" {
		MaxElementSize = ( 0.08 0.08 )
		MinElementSize = ( 0.01 0.01 )
	}
	Refinement "JunctionMesh" {
		MaxElementSize = ( 0.02 0.02 )
		MinElementSize = ( 0.002 0.002 )
	}
}

Placements {
	Constant "P_Doping_Reg" {
		Reference = "P_Doping"
		EvaluateWindow {
			Element = region ["R.Si"]
		}
	}
	Constant "N_Doping_Reg" {
		Reference = "N_Doping"
		EvaluateWindow {
			Element = region ["R.NRegion"]
		}
	}
	Refinement "GlobalMesh_Reg" {
		Reference = "GlobalMesh"
		RefineWindow = region ["R.Si"]
	}
	Refinement "JunctionMesh_Placement" {
		Reference = "JunctionMesh"
		RefineWindow = Rectangle [(0.85 0) (1.15 0.5)]
	}
}

