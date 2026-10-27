#pragma once
// Diagnostic controller only. Physical residuals and final stopping gates belong
// to the caller; raw update RMS deliberately is not a mode-selection input.
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace vela::diagnostics {
enum class DerivativeMode { Bootstrap, Trial, Full };
inline const char* derivativeModeName(DerivativeMode mode) {
    switch(mode) {
    case DerivativeMode::Bootstrap:return "bootstrap";
    case DerivativeMode::Trial:return "trial";
    case DerivativeMode::Full:return "full";
    }
    return "unknown";
}
struct DerivativeSwitchConfig {
    int minimumTrialSteps=2,maximumTrialSteps=4,failureConfirmations=2;
    int bootstrapProbeInterval=2,cooldownSteps=1;
    double minimumModelRatio=.1;
};
struct DerivativeObservation {
    bool valid=false,physicalGates=false;
    double before=0.,after=0.,modelAfter=0.;
};
struct DerivativeQuality {
    bool reliable=false,adequate=false;
    double predicted=0.,actual=0.,noise=0.,ratio=0.;
};
inline DerivativeQuality derivativeQuality(const DerivativeObservation& o,double threshold) {
    DerivativeQuality q;
    if(!o.valid||!std::isfinite(o.before)||!std::isfinite(o.after)||!std::isfinite(o.modelAfter))return q;
    q.predicted=o.before-o.modelAfter;q.actual=o.before-o.after;
    q.noise=64.*std::numeric_limits<double>::epsilon()*
        std::max({std::abs(o.before),std::abs(o.after),std::abs(o.modelAfter),std::numeric_limits<double>::min()});
    q.reliable=q.predicted>q.noise;
    if(q.reliable)q.ratio=q.actual/q.predicted;
    // Near the original physical gates, do not divide near-zero decreases to
    // veto a feasible full step. Final RMS still must be checked by the caller.
    q.adequate=o.physicalGates||(q.reliable&&q.actual>q.noise&&q.ratio>=threshold);
    return q;
}
class DerivativeSwitch {
    DerivativeSwitchConfig cfg_;
    DerivativeMode mode_=DerivativeMode::Bootstrap;
    int bootstrap_=0,cooldown_=0,trial_=0,good_=0,bad_=0;
    void fallback(const std::string& why) {
        mode_=DerivativeMode::Bootstrap;bootstrap_=trial_=good_=bad_=0;
        cooldown_=cfg_.cooldownSteps+1;reason=why;
    }
public:
    std::string reason="bootstrap_start";
    explicit DerivativeSwitch(DerivativeSwitchConfig c={},bool startFull=false):cfg_(c) {
        if(c.minimumTrialSteps<1||c.maximumTrialSteps<c.minimumTrialSteps||c.maximumTrialSteps>100||
           c.failureConfirmations<2||c.bootstrapProbeInterval<1||c.cooldownSteps<0||
           !(c.minimumModelRatio>0.&&c.minimumModelRatio<1.))throw std::invalid_argument("invalid derivative switch configuration");
        if(startFull){mode_=DerivativeMode::Full;reason="carried_full_mode";}
    }
    DerivativeMode mode()const{return mode_;}
    int trialSteps()const{return trial_;}
    int badSteps()const{return bad_;}
    bool requestFull(bool beforePhysicalGates) {
        if(mode_!=DerivativeMode::Bootstrap)return true;
        if(cooldown_>0){reason="bootstrap_cooldown";return false;}
        if(beforePhysicalGates||bootstrap_>=cfg_.bootstrapProbeInterval){
            mode_=DerivativeMode::Trial;trial_=good_=bad_=0;
            reason=beforePhysicalGates?"physical_gates_trigger_restore":"periodic_restore_probe";return true;
        }
        reason="bootstrap_progress";return false;
    }
    void bootstrapApplied(){++bootstrap_;if(cooldown_>0)--cooldown_;}
    // false means reject this candidate before changing the physical state.
    bool fullCandidate(const DerivativeObservation& o) {
        if(!o.valid){fallback("invalid_full_candidate");return false;}
        if(mode_==DerivativeMode::Trial)++trial_;
        const auto q=derivativeQuality(o,cfg_.minimumModelRatio);
        if(q.adequate){++good_;bad_=0;}else{++bad_;good_=0;}
        if(bad_>=cfg_.failureConfirmations){fallback("confirmed_model_failure");return false;}
        if(mode_==DerivativeMode::Trial&&trial_>=cfg_.minimumTrialSteps&&good_>=cfg_.minimumTrialSteps){
            mode_=DerivativeMode::Full;reason="full_mode_confirmed";return true;
        }
        if(mode_==DerivativeMode::Trial&&trial_>=cfg_.maximumTrialSteps){fallback("restore_window_exhausted");return false;}
        reason=q.adequate?"full_step_retained":"observing_model_failure";return true;
    }
};
} // namespace vela::diagnostics
