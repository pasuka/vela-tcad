"""Private fixed-state diagnostic formulas; not a production mobility profile.

Inputs are cm, V, cm^-3, cm^2/(V s), and V/cm at 300 K. The geometric
projection must come from true transport boundary faces. Partials are branch
local; a step crossing the 1 V/cm HFS cutoff needs separate treatment.
"""
import math
import numpy as np


def surface_hfs(bulk, normal_field, distance_cm, impurity_cm3, drive, electron):
    """Return mobility, d(mu)/d(normal field), and d(mu)/d(HFS drive)."""
    if not all(math.isfinite(x) for x in (bulk, normal_field, distance_cm, impurity_cm3, drive)):
        raise ValueError('finite mobility operands required')
    if bulk <= 0 or min(normal_field, distance_cm, impurity_cm3, drive) < 0:
        raise ValueError('positive bulk and nonnegative physical operands required')
    B,C,lam,delta,eta = ((4.75e5,580.*4.641588833612778e-4,.125,5.82e10,5.82e32)
                        if electron else (9.925e4,2947.*4.641588833612778e-4,.0317,2.0546e10,2.0546e32))
    # At zero normal magnitude this is the right derivative. Composing with
    # abs(E.normal) still has a cusp, which must not be certified as smooth.
    low, dn = bulk, -(bulk*1e-4)**2*math.exp(-distance_cm*.01/1e-8)*1e6/B
    if normal_field > 0:
        F=normal_field*100.; damp=math.exp(-distance_cm*.01/1e-8)
        # Both concentration operands are SI; N0=1 cm^-3=1e6 m^-3.
        acoustic2=C*((impurity_cm3*1e6+1e6)/1e6)**lam/F**(1/3)
        acoustic=B/F+acoustic2
        rough=(F/100.)**2/delta+F**3/eta
        si=1/(1/(bulk*1e-4)+damp*(1/acoustic+rough))
        dacoustic=-B/F**2-acoustic2/(3*F)
        drough=2*F/(1e4*delta)+3*F**2/eta
        low=si*1e4
        dn=-si*si*damp*(-dacoustic/acoustic**2+drough)*1e6
    if drive < 1.:
        return low,dn,0.
    beta,vs=(1.109,1.07e7) if electron else (1.213,8.37e6)
    z=low*drive/vs; power=z**beta
    dm=(1+power)**(-1-1/beta)
    return low/(1+power)**(1/beta),dm*dn,-low*low/vs*z**(beta-1)*dm


def triangle_mobility(points_cm, psi, qf, distance_cm, bulk, impurity_cm3,
                      projection, contact, electron):
    """Live local values and exact psi/QF partials at fixed bulk/doping.

    This is a local constitutive candidate, not a full coupled DD Jacobian:
    population/statistics chain rules and source derivatives are not included.
    """
    p=np.asarray(points_cm,float); ps=np.asarray(psi,float); q=np.asarray(qf,float)
    d=np.asarray(distance_cm,float); P=np.asarray(projection,float)
    if p.shape!=(3,2) or P.shape!=(2,2) or any(np.shape(a)!=(3,) for a in (ps,q,d,bulk,impurity_cm3)):
        raise ValueError('triangle and three local operands required')
    if not np.allclose(P@P,P,rtol=0,atol=1e-12) or not np.allclose(P,P.T,rtol=0,atol=1e-12):
        raise ValueError('orthogonal projection required')
    inv=np.linalg.inv(p[1:]-p[0]); G=np.column_stack((-inv[:,0]-inv[:,1],inv[:,0],inv[:,1]))
    # Subtract endpoints before multiplying: common potential/QF offsets
    # must not enter the dot-product cancellation path of a near-zero field.
    E=-inv@(ps[1:]-ps[0]); dg=inv@(d[1:]-d[0]); signed=float(E@dg); normal=abs(signed)
    force=inv@(ps[1:]-ps[0]) if contact else P@inv@(q[1:]-q[0]); drive=float(np.linalg.norm(force))
    dnormal=-np.sign(signed)*dg@G
    ddrive=force@(G if contact else P@G)/drive if drive else np.zeros(3)
    values=[]; jp=[]; jq=[]
    for k in range(3):
        mu,dn,df=surface_hfs(bulk[k],normal,d[k],impurity_cm3[k],drive,electron)
        values.append(mu);jp.append(dn*dnormal+(df*ddrive if contact else 0));jq.append(np.zeros(3) if contact else df*ddrive)
    return np.array(values),np.array(jp),np.array(jq)
