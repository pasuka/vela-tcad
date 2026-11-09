(sdegeo:set-default-boolean "ABA")
(sdegeo:create-rectangle
  (position 0.0 0.0 0.0) (position 2.0 1.0 0.0) "Silicon" "R.P")
(sdegeo:create-rectangle
  (position 1.0 0.0 0.0) (position 3.0 1.0 0.0) "Silicon" "R.N")

(sdegeo:insert-vertex (position 0.0 0.5 0.0))
(sdegeo:define-contact-set "Anode" 4 (color:rgb 1 0 0) "##")
(sdegeo:set-contact (find-edge-id (position 0.0 0.25 0.0)) "Anode")
(sdegeo:define-contact-set "Cathode" 4 (color:rgb 0 0 1) "##")
(sdegeo:set-contact (find-edge-id (position 3.0 0.5 0.0)) "Cathode")

(sdedr:define-constant-profile
  "Background" "PhosphorusActiveConcentration" 1.0e14)
(sdedr:define-constant-profile-material
  "Background.Place" "Background" "Silicon")
(sdedr:define-constant-profile
  "P.Profile" "BoronActiveConcentration" 1.0e16)
(sdedr:define-constant-profile-region
  "P.Place" "P.Profile" "R.P")
(sdedr:define-constant-profile
  "N.Profile" "PhosphorusActiveConcentration" 2.0e16)
(sdedr:define-constant-profile-region
  "N.Place" "N.Profile" "R.N" 0.0 "LocalReplace")

(sdedr:define-refeval-window
  "N.Baseline" "Line" (position 1.5 0.8 0.0) (position 2.5 0.8 0.0))
(sdedr:define-gaussian-profile
  "N.Gaussian" "ArsenicActiveConcentration"
  "PeakPos" 0.2 "PeakVal" 1.0e18
  "ValueAtDepth" 1.0e16 "Depth" 0.4 "Gauss" "Factor" 1.0)
(sdedr:define-analytical-profile-placement
  "N.Gaussian.Place" "N.Gaussian" "N.Baseline"
  "Positive" "NoReplace" "NoEval" "R.N" 0.0 "region")

; Offline Gmsh ignores the requested Sentaurus mesh build.  The licensed VM
; fixture runner uses it to freeze a T-2022.03-SP2 mesh oracle for this exact
; owned input without maintaining a separate commercial-only deck.
(sde:build-mesh "snmesh" "" "m2_core")
