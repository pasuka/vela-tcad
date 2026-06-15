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

(sdegeo:define-2d-contact
  (find-edge-id
    (position 0.0 (/ H 2.0) 0.0)
  )
  "Anode"
)

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
  1e17
)

(sdedr:define-constant-profile
  "N.Doping"
  "PhosphorusActiveConcentration"
  1e17
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

(sde:build-mesh
  "snmesh"
  ""
  "pn2d"
)
