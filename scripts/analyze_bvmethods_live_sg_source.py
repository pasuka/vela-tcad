"""Fixed-state SG/avalanche audit. Matching a vector is not element identity.

Native coordinates are um, densities cm^-3, fields V/cm, J/mu is
A/(cm^2)/(cm^2/V/s), alpha cm^-1. No PDE or trajectory qualification.
"""
import math
import numpy as np


def bernoulli(x):
    if abs(x)<1e-5:
        return 1-x/2+x*x/12-x**4/720
    if x>600:
        return x*math.exp(-x)
    if x< -600:
        return -x
    return x/math.expm1(x)


def sg_edge_current(density0, density1, qf_drop, length_cm, sign,
                    vt=.02585199526648486, q=1.602192e-19, g=1.):
    """Cancellation-free current without mobility, along endpoint 0 -> 1."""
    if min(density0,density1,length_cm,vt,q,g)<=0 or sign not in [-1,1]:
        raise ValueError('positive populations, length, constants and g required')
    delta=sign*qf_drop/(vt*g)
    arg=math.log(density1/density0)+delta
    # Algebraically equivalent endpoint forms avoid overflow for large drops.
    if delta>0:
        value=density0*bernoulli(-arg)*(-math.expm1(-delta))
    else:
        value=density1*bernoulli(arg)*math.expm1(delta)
    return -sign*q*vt*g*value/length_cm


def representative_current(points_cm, density, qf, psi, rep, sign,
                           eta=None, boundary_normals=(), is_contact=False):
    points_cm=np.asarray(points_cm,float);density=np.asarray(density,float)
    others=[j for j in range(3) if j!=rep]
    vectors=points_cm[others]-points_cm[rep]
    lengths=np.linalg.norm(vectors,axis=1);directions=vectors/lengths[:,None]
    values=[]
    for j,h in zip(others,lengths):
        lr=math.log(density[j]/density[rep]);g=1.
        if eta is not None and abs(lr)>1e-10:
            g=(eta[j]-eta[rep])/lr
        values.append(sg_edge_current(density[rep],density[j],qf[j]-qf[rep],h,sign,g=g))
    current=np.linalg.solve(directions,values)
    field=-np.linalg.solve(vectors,np.asarray(psi)[others]-psi[rep])
    # At electrode/insulator junctions the electrode vertex retains both
    # current components. This is an observed diagnostic support contract.
    if len(boundary_normals) and not is_contact:
        ns=np.asarray(boundary_normals,float)
        projection=np.eye(2)-np.linalg.pinv(ns)@ns
        current=projection@current;field=projection@field
    magnitude=np.linalg.norm(current)
    force=max(0.,float(field@current/magnitude)) if magnitude else 0.
    return current,force


def van_overstraeten(field, carrier, temperature=300.):
    if carrier not in ['electron','hole'] or temperature<=0:
        raise ValueError('invalid carrier/temperature')
    if field<=0:
        return 0.
    low=field<4e5
    a=7.03e5 if carrier=='electron' else (1.582e6 if low else 6.71e5)
    b=1.231e6 if carrier=='electron' else (2.036e6 if low else 1.693e6)
    gamma=math.tanh(.063/(2*8.6173324e-5*300))/math.tanh(.063/(2*8.6173324e-5*temperature))
    return gamma*a*math.exp(-gamma*b/field)


def match_candidates(observed, candidates, relative_tolerance=1e-8, absolute_tolerance=0.):
    """Return every admissible identity; never turn an argmin into uniqueness.

    Tolerances describe a fingerprint's numerical resolution, not solver gates.
    Force/alpha are deliberately excluded from vector-based identification.
    """
    observed=np.asarray(observed,float)
    if observed.shape!=(2,) or not np.all(np.isfinite(observed)) or not candidates:
        raise ValueError('finite vector and nonempty candidates required')
    errors=np.array([np.linalg.norm(np.asarray(c['current'])-observed) for c in candidates])
    limit=absolute_tolerance+relative_tolerance*np.linalg.norm(observed)
    accepted=[candidates[i]['identity'] for i in np.flatnonzero(errors<=limit)]
    best=int(np.argmin(errors))
    return dict(identities=accepted, unique=len(accepted)==1,
                best_identity=candidates[best]['identity'], absolute_error=float(errors[best]),
                relative_error=float(errors[best]/max(np.linalg.norm(observed),1e-300)))
