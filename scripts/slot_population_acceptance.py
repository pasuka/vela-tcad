"""Opt-in SLOT population-driver acceptance; not a global solver default."""
import math


def accept(rho, raw_error, *, closure, boundary, linear, mapping, positive,
           policy="primary_with_guards"):
    """Use the unprojected physical update error, never the damped step."""
    if policy not in ("primary_only", "primary_with_guards"):
        raise ValueError("unknown SLOT population acceptance policy")
    if not all(math.isfinite(x) and x >= 0 for x in (rho, raw_error)):
        return False
    return bool((rho <= 1e-10 or raw_error < 1)
                and (policy == "primary_only" or (closure and boundary))
                and linear and mapping and positive)


def readback_audit(rho, readback_rho, rhs_converged):
    """A saved-state consistency check, not a second nonlinear RHS gate."""
    finite = all(math.isfinite(x) and x >= 0 for x in (rho, readback_rho))
    tolerance = max(1e-15, 1e-6 * abs(rho)) if finite else 0.
    return dict(rhs_pass=bool(rhs_converged), residual=readback_rho,
                consistency_tolerance=tolerance,
                consistent=bool(finite and abs(rho-readback_rho) <= tolerance))
