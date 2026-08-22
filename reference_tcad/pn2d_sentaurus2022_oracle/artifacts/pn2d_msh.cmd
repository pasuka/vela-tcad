Title "Untitled"

Controls {
}

IOControls {
	EnableSections
}

Definitions {
	Constant "P.Doping" {
		Species = "BoronActiveConcentration"
		Value = 1e+17
	}
	Constant "N.Doping" {
		Species = "PhosphorusActiveConcentration"
		Value = 1e+17
	}
	Refinement "Global.Mesh" {
		MaxElementSize = ( 0.05 0.05 )
		MinElementSize = ( 0.01 0.01 )
	}
	Refinement "Junction.Mesh" {
		MaxElementSize = ( 0.01 0.02 )
		MinElementSize = ( 0.005 0.005 )
	}
}

Placements {
	Constant "P.Place" {
		Reference = "P.Doping"
		EvaluateWindow {
			Element = Rectangle [(0 0) (1 0.5)]
		}
	}
	Constant "N.Place" {
		Reference = "N.Doping"
		EvaluateWindow {
			Element = Rectangle [(1 0) (2 0.5)]
		}
	}
	Refinement "Global.Mesh.Place" {
		Reference = "Global.Mesh"
		RefineWindow = region ["R.Si"]
	}
	Refinement "Junction.Mesh.Place" {
		Reference = "Junction.Mesh"
		RefineWindow = Rectangle [(0.9 0) (1.1 0.5)]
	}
}

