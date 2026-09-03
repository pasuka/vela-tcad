(sde:clear)
(sdegeo:set-default-boolean "ABA")

; Local-refinement sensitivity variant of bjt_sde.cmd.
; Geometry, contacts, and doping are identical to the accepted baseline.
(define silicon
  (sdegeo:create-rectangle
    (position 0.0 0.0 0.0)
    (position 6.0 2.0 0.0)
    "Silicon" "R.Silicon"))

(sdegeo:insert-vertex (position 1.25 0.0 0.0))
(sdegeo:insert-vertex (position 2.00 0.0 0.0))
(sdegeo:insert-vertex (position 2.75 0.0 0.0))
(sdegeo:insert-vertex (position 4.25 0.0 0.0))

(sdegeo:define-contact-set "base" 4.0 (color:rgb 0.0 0.0 1.0) "##")
(sdegeo:define-2d-contact (find-edge-id (position 1.625 0.0 0.0)) "base")
(sdegeo:define-contact-set "emitter" 4.0 (color:rgb 1.0 0.0 0.0) "##")
(sdegeo:define-2d-contact (find-edge-id (position 3.50 0.0 0.0)) "emitter")
(sdegeo:define-contact-set "collector" 4.0 (color:rgb 0.0 0.6 0.0) "##")
(sdegeo:define-2d-contact (find-edge-id (position 3.00 2.0 0.0)) "collector")

(sdedr:define-constant-profile
  "Doping.BackgroundDonor" "PhosphorusActiveConcentration" 5.0e15)
(sdedr:define-constant-profile-region
  "Place.BackgroundDonor" "Doping.BackgroundDonor" "R.Silicon")

(sdedr:define-refeval-window
  "Ref.Doping" "Line" (position 0.0 0.0 0.0) (position 6.0 0.0 0.0))

(sdedr:define-analytical-profile
  "Doping.BaseBulk" "BoronActiveConcentration"
  "bbxm=3.0;bbxh=1.75;bbym=0.175;bbyh=0.175;bbxc=0.12;bbyc=0.16"
  "6.0e17*exp(-(((abs(x-bbxm)-bbxh+abs(abs(x-bbxm)-bbxh))/(2*bbxc))^2)-(((abs(y-bbym)-bbyh+abs(abs(y-bbym)-bbyh))/(2*bbyc))^2))" 0.0 "General")
(sdedr:define-analytical-profile-placement
  "Place.BaseBulk" "Doping.BaseBulk" "Ref.Doping" "Both" "NoReplace" "Eval")

(sdedr:define-analytical-profile
  "Doping.BaseSurface" "BoronActiveConcentration"
  "bsxm=3.0;bsxh=1.75;bsxc=0.12;bsyc=0.16"
  "4.0e18*exp(-(((abs(x-bsxm)-bsxh+abs(abs(x-bsxm)-bsxh))/(2*bsxc))^2)-((abs(y)/bsyc)^2))" 0.0 "General")
(sdedr:define-analytical-profile-placement
  "Place.BaseSurface" "Doping.BaseSurface" "Ref.Doping" "Both" "NoReplace" "Eval")

(sdedr:define-analytical-profile
  "Doping.Emitter" "PhosphorusActiveConcentration"
  "emxm=3.5;emxh=0.75;emxc=0.1275;emyc=0.17"
  "7.0e19*exp(-(((abs(x-emxm)-emxh+abs(abs(x-emxm)-emxh))/(2*emxc))^2)-((abs(y)/emyc)^2))" 0.0 "General")
(sdedr:define-analytical-profile-placement
  "Place.Emitter" "Doping.Emitter" "Ref.Doping" "Both" "NoReplace" "Eval")

(sdedr:define-analytical-profile
  "Doping.Collector" "PhosphorusActiveConcentration"
  "coym=2.0;coyc=0.27"
  "1.0e19*exp(-((abs(y-coym)/coyc)^2))" 0.0 "General")
(sdedr:define-analytical-profile-placement
  "Place.Collector" "Doping.Collector" "Ref.Doping" "Both" "NoReplace" "Eval")

; Preserve the accepted baseline mesh definitions.
(sdedr:define-refinement-size "RefDef.Global" 0.15 0.10 0.015 0.008)
(sdedr:define-refinement-function
  "RefDef.Global" "DopingConcentration" "MaxTransDiff" 0.8)
(sdedr:define-refinement-region "PlaceRef.Global" "RefDef.Global" "R.Silicon")

(sdedr:define-refeval-window
  "RefWin.Top" "Rectangle" (position 1.0 0.0 0.0) (position 5.0 0.85 0.0))
(sdedr:define-refinement-size "RefDef.Top" 0.05 0.035 0.008 0.004)
(sdedr:define-refinement-placement "PlaceRef.Top" "RefDef.Top" "RefWin.Top")

(sdedr:define-refeval-window
  "RefWin.Bottom" "Rectangle" (position 0.0 1.45 0.0) (position 6.0 2.0 0.0))
(sdedr:define-refinement-size "RefDef.Bottom" 0.10 0.04 0.012 0.004)
(sdedr:define-refinement-placement "PlaceRef.Bottom" "RefDef.Bottom" "RefWin.Bottom")

; Refine the p-base / base-collector depletion band where the weak-current
; discrepancy is concentrated (144/180 tail nodes were at y=0.50-0.75 um).
(sdedr:define-refeval-window
  "RefWin.BaseCollector" "Rectangle" (position 0.75 0.42 0.0) (position 5.25 0.88 0.0))
(sdedr:define-refinement-size "RefDef.BaseCollector" 0.025 0.0175 0.004 0.002)
(sdedr:define-refinement-placement
  "PlaceRef.BaseCollector" "RefDef.BaseCollector" "RefWin.BaseCollector")

; Resolve the lateral emitter-base junction flanks without globally refining
; the heavily doped emitter interior.
(sdedr:define-refeval-window
  "RefWin.EmitterBaseLeft" "Rectangle" (position 2.52 0.0 0.0) (position 2.98 0.55 0.0))
(sdedr:define-refinement-size "RefDef.EmitterBaseLeft" 0.015 0.0125 0.002 0.0015)
(sdedr:define-refinement-placement
  "PlaceRef.EmitterBaseLeft" "RefDef.EmitterBaseLeft" "RefWin.EmitterBaseLeft")

(sdedr:define-refeval-window
  "RefWin.EmitterBaseRight" "Rectangle" (position 4.02 0.0 0.0) (position 4.48 0.55 0.0))
(sdedr:define-refinement-size "RefDef.EmitterBaseRight" 0.015 0.0125 0.002 0.0015)
(sdedr:define-refinement-placement
  "PlaceRef.EmitterBaseRight" "RefDef.EmitterBaseRight" "RefWin.EmitterBaseRight")

(sde:build-mesh "snmesh" "-a -c boxmethod" "bjt")
