#pragma once
#include <Eigen/QR>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace vela::experimental {
struct ElectrothermalLocalWeights {
    Eigen::VectorXd weights;
    unsigned degree=0;
    double amplification=0.;
};
// Anchor the fit at the latest accepted state; fit changes on normalized bias.
// Limit noise amplification, downgrading a quadratic fit before giving up.
inline ElectrothermalLocalWeights electrothermalLocalWeights(const std::vector<double>& bias,double target) {
    if(bias.size()<3 || bias.size()>4 || !std::isfinite(target))
        throw std::invalid_argument("Local prediction requires three or four history points");
    for(std::size_t i=0;i<bias.size();++i)
        if(!std::isfinite(bias[i]) || (i && !(bias[i]>bias[i-1])))
            throw std::invalid_argument("Local prediction history must increase");
    const double span=bias.back()-bias.front(),q=(target-bias.back())/span;
    if(!(span>1e-12 && q>0. && q<=3.))throw std::invalid_argument("Local prediction span or reach invalid");
    const Eigen::Index count=static_cast<Eigen::Index>(bias.size()-1);
    for(unsigned degree:{2u,1u}) {
        Eigen::MatrixXd design(count,degree);Eigen::VectorXd evaluation(degree);
        for(Eigen::Index i=0;i<count;++i) {
            const double z=(bias[i]-bias.back())/span;
            design(i,0)=z;if(degree==2)design(i,1)=z*z;
        }
        evaluation[0]=q;if(degree==2)evaluation[1]=q*q;
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(design);qr.setThreshold(1e-8);
        if(qr.rank()!=degree)continue;
        Eigen::VectorXd weights(count+1);
        weights.head(count)=qr.solve(Eigen::MatrixXd::Identity(count,count)).transpose()*evaluation;
        weights[count]=1.-weights.head(count).sum();
        const double amplification=weights.lpNorm<1>();
        if(weights.allFinite() && amplification<=8.)return {std::move(weights),degree,amplification};
    }
    throw std::invalid_argument("Local prediction would amplify history noise excessively");
}
}
