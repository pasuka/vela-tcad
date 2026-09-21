;==========================================================
; 2D PN Junction
; Sentaurus Structure Editor
;
; Purpose:
;   Generate a simple 2D PN junction
;
; Unit:
;   Length = um
;   Doping = cm^-3
;==========================================================

(sde:clear)

;----------------------------------------------------------
; Device geometry
;----------------------------------------------------------

(define L 2.0)       ; total device length
(define H 0.5)       ; device height
(define XJ 1.0)      ; PN junction position

;----------------------------------------------------------
; Create silicon region
;----------------------------------------------------------

(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon"
  "R.Si"
)

;----------------------------------------------------------
(sdegeo:insert-vertex (position 0.0 0.125 0.0))
(sdegeo:insert-vertex (position 0.0 0.1875 0.0))
(sdegeo:insert-vertex (position 0.0 0.3125 0.0))
(sdegeo:insert-vertex (position 0.0 0.375 0.0))
; Define contacts
;----------------------------------------------------------

(sdegeo:define-contact-set
  "Anode"
  4.0
  (color:rgb 1 0 0)
  "##"
)

(sdegeo:define-contact-set
  "Cathode"
  4.0
  (color:rgb 0 0 1)
  "##"
)

; Left electrode

(sdegeo:set-current-contact-set "Anode")

(sdegeo:define-2d-contact (find-edge-id (position 0.0 0.0625 0.0)) "Anode")
(sdegeo:define-2d-contact (find-edge-id (position 0.0 0.15625 0.0)) "Anode")
(sdegeo:define-2d-contact (find-edge-id (position 0.0 0.25 0.0)) "Anode")
(sdegeo:define-2d-contact (find-edge-id (position 0.0 0.34375 0.0)) "Anode")
(sdegeo:define-2d-contact (find-edge-id (position 0.0 0.4375 0.0)) "Anode")

; Right electrode

(sdegeo:set-current-contact-set "Cathode")

(sdegeo:define-2d-contact
  (find-edge-id
    (position L (/ H 2.0) 0.0)
  )
  "Cathode"
)

;----------------------------------------------------------
; Doping windows
;----------------------------------------------------------

(sdedr:define-refeval-window
  "P.Window"
  "Rectangle"
  (position 0.0 0.0 0.0)
  (position XJ H 0.0)
)

(sdedr:define-refeval-window
  "N.Window"
  "Rectangle"
  (position XJ 0.0 0.0)
  (position L H 0.0)
)

;----------------------------------------------------------
; Constant doping profiles
;----------------------------------------------------------

(sdedr:define-constant-profile
  "P.Doping"
  "BoronActiveConcentration"
  1e+17
)

(sdedr:define-constant-profile
  "N.Doping"
  "PhosphorusActiveConcentration"
  1e+17
)

;----------------------------------------------------------
; Apply doping
;----------------------------------------------------------

(sdedr:define-constant-profile-placement
  "P.Place"
  "P.Doping"
  "P.Window"
)

(sdedr:define-constant-profile-placement
  "N.Place"
  "N.Doping"
  "N.Window"
)

;----------------------------------------------------------
; Global mesh
;----------------------------------------------------------

(sdedr:define-refinement-size
  "Global.Mesh"
  0.05 0.05
  0.01 0.01
)

(sdedr:define-refinement-region
  "Global.Mesh.Place"
  "Global.Mesh"
  "R.Si"
)

;----------------------------------------------------------
; Junction mesh refinement
;----------------------------------------------------------

(sdedr:define-refeval-window
  "Junction.Window"
  "Rectangle"
  (position 0.9 0.0 0.0)
  (position 1.1 H 0.0)
)

(sdedr:define-refinement-size
  "Junction.Mesh"
  0.01 0.02
  0.005 0.005
)

(sdedr:define-refinement-placement
  "Junction.Mesh.Place"
  "Junction.Mesh"
  "Junction.Window"
)

;----------------------------------------------------------
; Build mesh
;----------------------------------------------------------


(sdedr:define-refeval-window "Edge0" "Rectangle" (position 0.0 0.1 0.0) (position 0.05 0.15 0.0))
(sdedr:define-refinement-size "EdgeMesh0" 0.01 0.00625 0.005 0.003125)
(sdedr:define-refinement-placement "EdgePlace0" "EdgeMesh0" "Edge0")

(sdedr:define-refeval-window "Edge1" "Rectangle" (position 0.0 0.1625 0.0) (position 0.05 0.2125 0.0))
(sdedr:define-refinement-size "EdgeMesh1" 0.01 0.00625 0.005 0.003125)
(sdedr:define-refinement-placement "EdgePlace1" "EdgeMesh1" "Edge1")

(sdedr:define-refeval-window "Edge2" "Rectangle" (position 0.0 0.2875 0.0) (position 0.05 0.3375 0.0))
(sdedr:define-refinement-size "EdgeMesh2" 0.01 0.00625 0.005 0.003125)
(sdedr:define-refinement-placement "EdgePlace2" "EdgeMesh2" "Edge2")

(sdedr:define-refeval-window "Edge3" "Rectangle" (position 0.0 0.35 0.0) (position 0.05 0.4 0.0))
(sdedr:define-refinement-size "EdgeMesh3" 0.01 0.00625 0.005 0.003125)
(sdedr:define-refinement-placement "EdgePlace3" "EdgeMesh3" "Edge3")

(sde:build-mesh
  "snmesh"
  ""
  "pn2d"
)
