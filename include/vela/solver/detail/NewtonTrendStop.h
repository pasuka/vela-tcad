#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace vela::detail {
// Failure-only diagnostic policy. The caller must give existing success gates priority.
struct NewtonTrendSample {
    int iteration = 0;
    std::array<double, 6> rows{}; // electron N/E/max, hole N/E/max; passing rows contribute zero
    std::array<double, 6> gates{}; // primary, e/h closure, branch QF, reserved; <=1 passes
    double rawDirection = 0;
    double alpha = 1;
};
class NewtonTrendStop {
    std::vector<NewtonTrendSample> samples_;
    int consecutive_ = 0;
    static double median(double a, double b, double c) {
        std::array<double, 3> x{a,b,c}; std::sort(x.begin(),x.end()); return x[1];
    }
    static bool improved(double before, double after, double fraction) {
        return before > 0 && (before-after)/before >= fraction;
    }
public:
    int consecutiveWindows() const { return consecutive_; }
    std::string observe(const NewtonTrendSample& s, bool qualified = false) {
        if (qualified) { consecutive_=0; return {}; }
        bool finite = std::isfinite(s.rawDirection) && std::isfinite(s.alpha);
        for (double v:s.rows) finite = finite && std::isfinite(v) && v>=0;
        for (double v:s.gates) finite = finite && std::isfinite(v) && v>=0;
        if (!finite) return "nonfinite_trend_metric";
        samples_.push_back(s);
        const int k = static_cast<int>(samples_.size())-1;
        // Early tiny-step stop: three consecutive accepted steps, no >=1% improvement
        // in any row-quality metric or any still-failing gate at each transition.
        if (k>=3) {
            bool tiny=true;
            for (int t=k-2;t<=k;++t) {
                const auto& a=samples_[t-1]; const auto& b=samples_[t];
                tiny = tiny && b.alpha<=1.0/1024;
                for (int j=0;j<6;++j) {
                    if (improved(a.rows[j],b.rows[j],.01)) tiny=false;
                    if (b.gates[j]>1 && improved(a.gates[j],b.gates[j],.01)) tiny=false;
                }
            }
            if (tiny) return "tiny_step_stagnation";
        }
        // Six samples spanning five updates: median(first3) vs median(last3).
        if (s.iteration<12 || k<5) { consecutive_=0; return {}; }
        bool stalled=true;
        for (int j=0;j<6;++j) {
            const double a=median(samples_[k-5].rows[j],samples_[k-4].rows[j],samples_[k-3].rows[j]);
            const double b=median(samples_[k-2].rows[j],samples_[k-1].rows[j],samples_[k].rows[j]);
            if (improved(a,b,.01)) stalled=false;
            const double ga=median(samples_[k-5].gates[j],samples_[k-4].gates[j],samples_[k-3].gates[j]);
            const double gb=median(samples_[k-2].gates[j],samples_[k-1].gates[j],samples_[k].gates[j]);
            // Include current value so a newly failing gate is not hidden by the median.
            if (s.gates[j]>1 && improved(ga,gb,.01)) stalled=false;
        }
        const double da=median(samples_[k-5].rawDirection,samples_[k-4].rawDirection,samples_[k-3].rawDirection);
        const double db=median(samples_[k-2].rawDirection,samples_[k-1].rawDirection,samples_[k].rawDirection);
        int small=0; for(int t=k-2;t<=k;++t) if(samples_[t].alpha<=1.0/128) ++small;
        const bool auxiliary=!improved(da,db,.10) || small>=2;
        consecutive_ = stalled && auxiliary ? consecutive_+1 : 0;
        return consecutive_>=3 ? "stagnation_detected" : "";
    }
};
}
