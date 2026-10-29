// Included in vela::detail after the local Tri3 AD helpers.
// Explicit source-current convention: SG with the statistical correction
// incorporated into the effective band potential, not a generalized-Einstein
// multiplier. Live Fermi populations and their full derivatives are inputs.
template<class S>
inline S phumobSourceLog(const S& x)
{
    const Real v=localAdValue(x);
    if (!(v>0.) || !std::isfinite(v))
        throw PopulationDomainError("SG source requires finite positive populations");
    S result(std::log(v));
    if constexpr (std::is_same_v<S,Tri3LocalForwardDual>)
        for(std::size_t i=0;i<result.derivative.size();++i)
            result.derivative[i]=x.derivative[i]/v;
    return result;
}

template<class S>
inline S phumobSourceSgDrop(S left,S right,S dq,Real vt)
{
    // Returns the particle-continuity edge drop (without mobility or length).
    // Reverse the edge when necessary so expm1 never overflows. Do not discard
    // the derivative at dq=0: a zero current still has a linear response.
    Real sign=1.;
    if(localAdValue(dq)>0.) {std::swap(left,right);dq=-dq;sign=-1.;}
    const S arg=phumobSourceLog(right)-phumobSourceLog(left)+dq;
    return S(sign*vt)*right*localAdBernoulli(arg)*localAdExpm1(dq);
}

template<class S>
inline std::array<S,2> phumobSourceSgVector(
    const DeviceMesh& mesh,const Cell& cell,std::size_t vertex,
    const std::array<S,3>& density,const std::array<S,3>& qf,
    Real vt,bool electron)
{
    const std::size_t a=vertex,b=(a+1)%3,c=(a+2)%3;
    const auto& pa=mesh.getNode(cell.node_ids[a]);
    const auto& pb=mesh.getNode(cell.node_ids[b]);
    const auto& pc=mesh.getNode(cell.node_ids[c]);
    const Real dx1=pb.x-pa.x,dy1=pb.y-pa.y,dx2=pc.x-pa.x,dy2=pc.y-pa.y;
    const Real det=dx1*dy2-dx2*dy1;
    if(std::abs(det)<=1e-300 || !(vt>0.))
        throw std::invalid_argument("SG source requires a nondegenerate Tri3 and positive thermal voltage");
    const Real carrierSign=electron?1.:-1.;
    const auto edge=[&](std::size_t j) {
        return S(-carrierSign)*phumobSourceSgDrop(density[a],density[j],
            (qf[j]-qf[a])*S(carrierSign/vt),vt);
    };
    const S v1=edge(b),v2=edge(c);
    return {(v1*S(dy2)-v2*S(dy1))/S(det),
            (v2*S(dx1)-v1*S(dx2))/S(det)};
}
