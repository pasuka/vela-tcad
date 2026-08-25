(sde:clear)
(sdegeo:create-rectangle
  (position 0 0 0) (position 1 1 0) "Silicon" "R.Si")
(sdegeo:create-rectangle
  (position 1 0 0) (position 2 1 0) "Oxide" "R.Ox")

(sdegeo:define-contact-set "Left")
(sdegeo:set-contact (find-edge-id (position 0 0.5 0)) "Left")
(sdegeo:define-contact-set "Right")
(sdegeo:set-contact (find-edge-id (position 2 0.5 0)) "Right")

(sdedr:define-constant-profile
  "Si.Background" "PhosphorusActiveConcentration" 1e15)
(sdedr:define-constant-profile-region
  "Si.Background.Place" "Si.Background" "R.Si")

(sdedr:define-refeval-window
  "Full.Window" "Rectangle" (position 0 0 0) (position 2 1 0))
(sdedr:define-refeval-window
  "Polygon.Window" "Polygon"
  (list
    (position 0.2 0.2 0)
    (position 0.8 0.2 0)
    (position 0.6 0.8 0)
    (position 0.2 0.7 0)
    (position 0.2 0.2 0)))
(sdedr:define-refeval-window
  "Line.Window" "Line" (position 0.2 0.5 0) (position 0.8 0.5 0))

(sdedr:define-refinement-size "Material.Ref" 0.15 0.15 0.04 0.04)
(sdedr:define-refinement-material
  "Material.Ref.Place" "Material.Ref" "Silicon")

(sdedr:define-refinement-size "Polygon.Ref" 0.08 0.04 0.02 0.01)
(sdedr:define-refinement-placement
  "Polygon.Ref.Place" "Polygon.Ref" "Polygon.Window")

(sdedr:define-refinement-size "Line.Ref" 0.04 0.04 0.01 0.01)
(sdedr:define-refinement-placement
  "Line.Ref.Place" "Line.Ref" "Line.Window")

(sdedr:define-refinement-size "Interface.Ref" 0.15 0.15 0.02 0.02)
(sdedr:define-refinement-placement
  "Interface.Ref.Place" "Interface.Ref" "Full.Window")
(sdedr:define-refinement-function
  "Interface.Ref" "MaxLenInt" "R.Si" "R.Ox" 0.03 1.5
  "DoubleSide" "UseRegionNames")

(sde:build-mesh "snmesh" "" "m3_refinement")
