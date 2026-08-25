(sdegeo:set-default-boolean "ABA")

(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position 1.0 1.0 0.0)
  "Silicon" "R.P")

(sdegeo:create-rectangle
  (position 1.5 0.0 0.0)
  (position 2.5 1.0 0.0)
  "Silicon" "R.N")

(sdegeo:define-contact-set "Anode" 4.0 (color:rgb 1 0 0) "##")
(sdegeo:set-contact (find-edge-id (position 0.0 0.5 0.0)) "Anode")
(sdegeo:define-contact-set "Cathode" 4.0 (color:rgb 0 0 1) "##")
(sdegeo:set-contact (find-edge-id (position 2.5 0.5 0.0)) "Cathode")

(sdedr:define-constant-profile
  "P.Profile" "BoronActiveConcentration" 1.0e16)
(sdedr:define-constant-profile-region
  "P.Place" "P.Profile" "R.P" "NoReplace")
(sdedr:define-constant-profile
  "N.Profile" "PhosphorusActiveConcentration" 2.0e16)
(sdedr:define-constant-profile-region
  "N.Place" "N.Profile" "R.N" "NoReplace")
