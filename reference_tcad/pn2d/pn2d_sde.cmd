; pn2d_sde.cmd
; 2D PN junction diode structure for calibration

(sde:clear)

; -----------------------------
; Parameters
; -----------------------------
(define L  2.0)     ; um
(define H  0.5)     ; um
(define Xj 1.0)     ; junction position, um

(define Na 1e17)    ; P doping, cm-3
(define Nd 1e17)    ; N doping, cm-3

; -----------------------------
; Geometry
; -----------------------------
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L   H   0.0)
  "Silicon"
  "R.Si")

; -----------------------------
; Doping regions
; left: P, right: N
; -----------------------------
(sdedr:define-constant-profile "P_Doping" "BoronActiveConcentration" Na)
(sdedr:define-constant-profile-region "P_Doping_Reg" "P_Doping" "R.Si")

(sdedr:define-constant-profile "N_Doping" "PhosphorusActiveConcentration" Nd)

; N region overwrite on right side
(sdegeo:create-rectangle
  (position Xj 0.0 0.0)
  (position L  H   0.0)
  "Silicon"
  "R.NRegion")

(sdedr:define-constant-profile-region "N_Doping_Reg" "N_Doping" "R.NRegion")

; -----------------------------
; Contacts
; -----------------------------
(sdegeo:define-contact-set "Anode"   4.0 (color:rgb 1 0 0) "##")
(sdegeo:define-contact-set "Cathode" 4.0 (color:rgb 0 0 1) "##")

; left boundary = Anode
(sdegeo:set-current-contact-set "Anode")
(sdegeo:define-2d-contact
  (find-edge-id (position 0.0 (/ H 2.0) 0.0))
  "Anode")

; right boundary = Cathode
(sdegeo:set-current-contact-set "Cathode")
(sdegeo:define-2d-contact
  (find-edge-id (position L (/ H 2.0) 0.0))
  "Cathode")

; -----------------------------
; Mesh refinement
; -----------------------------
(sdedr:define-refinement-size "GlobalMesh"
  0.08 0.08
  0.01 0.01)

(sdedr:define-refinement-region "GlobalMesh_Reg" "GlobalMesh" "R.Si")

(sdedr:define-refinement-window "JunctionWindow"
  "Rectangle"
  (position (- Xj 0.15) 0.0 0.0)
  (position (+ Xj 0.15) H   0.0))

(sdedr:define-refinement-size "JunctionMesh"
  0.02 0.02
  0.002 0.002)

(sdedr:define-refinement-placement
  "JunctionMesh_Placement"
  "JunctionMesh"
  "JunctionWindow")

; -----------------------------
; Save structure and mesh
; -----------------------------
(sde:build-mesh "snmesh" "" "pn2d")
