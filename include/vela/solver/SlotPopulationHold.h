#pragma once
// Intentionally namespace-neutral: also included inside the runner's private
// namespace. Diagnostic policy only; it does not alter any physical equation.
struct SlotPopulationHoldState {
 bool held=false;
 int consecutiveNonpositive=0;
 int activations=0;
 int releases=0;
 void observe(double densityCm3, double thresholdCm3, double relative, int required) {
  // The caller validates finite inputs and a positive threshold/required count.
  if(held) {
   if(relative>0.) {held=false;consecutiveNonpositive=0;++releases;}
   return;
  }
  consecutiveNonpositive=(1.+relative<=0.)?consecutiveNonpositive+1:0;
  if(densityCm3<thresholdCm3 && consecutiveNonpositive>=required) {
   held=true;++activations;
  }
 }
};
