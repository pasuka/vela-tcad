;==========================================================
; 2D NMOS Id-Vd reference structure
; Sentaurus Structure Editor
;
; Purpose:
;   Generate a simple rectangular 2D NMOS cross-section for
;   Sentaurus/Vela Drain Id-Vd comparison.
;
; Unit:
;   Length = um
;   Doping = cm^-3
;==========================================================

(sde:clear)

; Aggregate silicon marker for import/test inventory: "R.Si".

;----------------------------------------------------------
; Device geometry
;----------------------------------------------------------

(define L 2.0)
(define HSi 0.4)
(define Tox 0.05)
(define XSource 0.35)
(define XDrain 1.65)

; Source, body/channel, and drain are separate rectangles so that
; contacts and region-average doping remain unambiguous after import.

(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position 0.35 0.4 0.0)
  "Silicon"
  "R.Source"
)

(sdegeo:create-rectangle
  (position 0.35 0.0 0.0)
  (position 1.65 0.4 0.0)
  "Silicon"
  "R.Body"
)

(sdegeo:create-rectangle
  (position 1.65 0.0 0.0)
  (position 2.0 0.4 0.0)
  "Silicon"
  "R.Drain"
)

(sdegeo:create-rectangle
  (position 0.35 0.4 0.0)
  (position 1.65 0.45 0.0)
  "Oxide"
  "R.Ox"
)

;----------------------------------------------------------
; Define contacts
;----------------------------------------------------------

(sdegeo:define-contact-set "Source" 4.0 (color:rgb 0 0.6 0) "##")
(sdegeo:define-contact-set "Drain"  4.0 (color:rgb 0 0 1) "##")
(sdegeo:define-contact-set "Gate"   4.0 (color:rgb 0.8 0.6 0) "##")
(sdegeo:define-contact-set "Body"   4.0 (color:rgb 1 0 0) "##")

(sdegeo:set-current-contact-set "Source")
(sdegeo:define-2d-contact (find-edge-id (position 0.175 0.4 0.0)) "Source")

(sdegeo:set-current-contact-set "Drain")
(sdegeo:define-2d-contact (find-edge-id (position 1.825 0.4 0.0)) "Drain")

(sdegeo:set-current-contact-set "Gate")
(sdegeo:define-2d-contact (find-edge-id (position 1.0 0.45 0.0)) "Gate")

(sdegeo:set-current-contact-set "Body")
(sdegeo:define-2d-contact (find-edge-id (position 1.0 0.0 0.0)) "Body")

;----------------------------------------------------------
; Doping windows
;----------------------------------------------------------

(sdedr:define-refeval-window
  "P.Body.Window"
  "Rectangle"
  (position 0.35 0.0 0.0)
  (position 1.65 0.4 0.0)
)

(sdedr:define-refeval-window
  "N.Source.Window"
  "Rectangle"
  (position 0.0 0.0 0.0)
  (position 0.35 0.4 0.0)
)

(sdedr:define-refeval-window
  "N.Drain.Window"
  "Rectangle"
  (position 1.65 0.0 0.0)
  (position 2.0 0.4 0.0)
)

;----------------------------------------------------------
; Constant doping profiles
;----------------------------------------------------------

(sdedr:define-constant-profile
  "P.Body.Doping"
  "BoronActiveConcentration"
  1e17
)

(sdedr:define-constant-profile
  "N.Source.Doping"
  "PhosphorusActiveConcentration"
  1e17
)

(sdedr:define-constant-profile
  "N.Drain.Doping"
  "PhosphorusActiveConcentration"
  1e17
)

(sdedr:define-constant-profile-placement
  "P.Body.Place"
  "P.Body.Doping"
  "P.Body.Window"
)

(sdedr:define-constant-profile-placement
  "N.Source.Place"
  "N.Source.Doping"
  "N.Source.Window"
)

(sdedr:define-constant-profile-placement
  "N.Drain.Place"
  "N.Drain.Doping"
  "N.Drain.Window"
)

;----------------------------------------------------------
; Mesh refinement
;----------------------------------------------------------

(sdedr:define-refinement-size
  "Global.Mesh"
  0.05 0.04
  0.01 0.01
)

(sdedr:define-refinement-region "Source.Mesh.Place" "Global.Mesh" "R.Source")
(sdedr:define-refinement-region "Body.Mesh.Place"   "Global.Mesh" "R.Body")
(sdedr:define-refinement-region "Drain.Mesh.Place"  "Global.Mesh" "R.Drain")
(sdedr:define-refinement-region "Ox.Mesh.Place"     "Global.Mesh" "R.Ox")

(sdedr:define-refeval-window
  "Channel.Window"
  "Rectangle"
  (position 0.30 0.34 0.0)
  (position 1.70 0.45 0.0)
)

(sdedr:define-refinement-size
  "Channel.Mesh"
  0.02 0.01
  0.005 0.0025
)

(sdedr:define-refinement-placement
  "Channel.Mesh.Place"
  "Channel.Mesh"
  "Channel.Window"
)

;----------------------------------------------------------
; Build mesh
;----------------------------------------------------------

(sde:build-mesh
  "snmesh"
  ""
  "nmos2d"
)