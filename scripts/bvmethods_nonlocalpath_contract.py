"""Independent N0 source/state checks; no model implementation or fitted data."""
import math

Q_C = 1.602176634e-19
H_J_S = 6.62607015e-34
KB_J_K = 1.380649e-23
M0_KG = 9.1093837015e-31


def indirect_parameters(A_cm3_s, B_V_cm, gap_reference_eV, phonon_eV, mass_ratio=0., *, constants=None):
    """Independent manual (515)-(516) conversion with CODATA 2018 constants.

    The reference gap is an explicit hypothesis, not inferred from a biased
    effective gap. This diagnostic is not a production NonLocalPath model.
    A/B normalization uses F0=1 V/cm. Returned masses are in units of m0 and
    Cpath in J^2 cm/kg. A zero ratio selects manual (506)-(507).
    Alternative version controls must supply the complete explicit constant
    set; this does not replace any production global physical constants.
    """
    values = [A_cm3_s, B_V_cm, gap_reference_eV, phonon_eV, mass_ratio]
    if any(not math.isfinite(x) for x in values) or min(values[:4]) <= 0 or mass_ratio < 0:
        raise ValueError("positive finite A/B/gap/phonon and nonnegative ratio required")
    physical = dict(q_C=Q_C, h_J_s=H_J_S, kb_J_K=KB_J_K, m0_kg=M0_KG)
    if constants is not None:
        if constants.keys() != physical.keys():
            raise ValueError("an explicit complete physical constant set is required")
        physical = dict(constants)
    if any(not math.isfinite(x) or x <= 0 for x in physical.values()):
        raise ValueError("physical constants must be finite and positive")
    q, h, kb, m0 = (physical[k] for k in ['q_C', 'h_J_s', 'kb_J_K', 'm0_kg'])
    eg, ep = gap_reference_eV * q, phonon_eV * q
    mr = (3 * q * h * (100 * B_V_cm) /
          (2**3.5 * math.pi * eg**1.5))**2
    if mass_ratio == 0:
        if 2 * mr >= m0:
            raise ValueError("automatic valence mass is nonpositive")
        mc = 1 / (1 / (2 * mr) + 1 / m0)
        mv = 1 / (1 / (2 * mr) - 1 / m0)
    else:
        mc, mv = mr * (1 + 1 / mass_ratio), mr * (1 + mass_ratio)
    phonons = 1 / math.expm1(ep / (kb * 300))
    c_si = (A_cm3_s * 1e6 * 2**(21/4) * h**2.5 * mr**1.25 * ep * eg**1.75 /
            ((mc * mv)**1.5 * (1 + 2 * phonons) * (q * 100)**2.5))
    return {"m_r_m0": mr / m0, "m_c_m0": mc / m0,
            "m_v_m0": mv / m0, "Cpath_J2_cm_kg": c_si / .01,
            "reference_gap_eV": gap_reference_eV, "phonon_eV": phonon_eV, "constants": physical}


def generation_current_per_um(rates_cm3_s, areas_um2):
    """Integrate a volume generation rate over a 2D device of 1 um thickness.

    The caller must identify the actual nodal integration support. This function
    does not infer boxes from triangles or claim that a TDR nodal sum is a source
    integral. Signed generation/recombination is retained.
    """
    rates, areas = list(rates_cm3_s), list(areas_um2)
    if len(rates) != len(areas) or not rates:
        raise ValueError("matching nonempty rate and area arrays required")
    if any(not math.isfinite(x) for x in rates + areas) or any(x < 0 for x in areas):
        raise ValueError("finite rates and nonnegative areas required")
    return Q_C * 1e-12 * math.fsum(g * a for g, a in zip(rates, areas))


def state_identity(reference, candidate):
    """Require the complete native labeled vector, not just selected weak rows.

    This establishes equality at the precision of the supplied representation;
    it cannot certify hidden precision of a formatted native debug vector.
    """
    if not reference or reference.keys() != candidate.keys():
        raise ValueError("incomplete or different native variable mapping")
    if any(not math.isfinite(x) for x in [*reference.values(), *candidate.values()]):
        raise ValueError("nonfinite state")
    changed = sum(reference[k] != candidate[k] for k in reference)
    return {"entries": len(reference), "identical": changed == 0,
            "changed": changed,
            "maximum_absolute_difference": max(abs(reference[k] - candidate[k]) for k in reference)}


def temperature_identity(actual_values_K, requested_K):
    """Inspect physical temperature, allowing four binary64 serialization ULPs.

    A normalized 400/300 saved value may return as 400.00000000000006 K.
    This representation allowance does not accept a different physical state.
    """
    values = list(actual_values_K)
    if not values or any(not math.isfinite(x) or x <= 0 for x in values + [requested_K]):
        raise ValueError("positive finite physical temperatures required")
    allowance = 4 * math.ldexp(1., math.frexp(requested_K)[1] - 53)
    return {"requested_K": requested_K, "minimum_K": min(values), "maximum_K": max(values),
            "serialization_allowance_K": allowance,
            "identical": all(abs(x - requested_K) <= allowance for x in values)}


def pair_balance(electron_integrals, hole_integrals, absolute_pair_events):
    """Normalize a signed pair ledger by absolute events, preserving tiny rates."""
    e, h, events = list(electron_integrals), list(hole_integrals), list(absolute_pair_events)
    if not e or not h or any(not math.isfinite(x) for x in e + h + events):
        raise ValueError("finite nonempty carrier integrals required")
    if any(x < 0 for x in events):
        raise ValueError("absolute event weights cannot be negative")
    error = abs(math.fsum(e) - math.fsum(h))
    denominator = math.fsum(events)
    if denominator == 0:
        if any(x != 0 for x in e + h):
            raise ValueError("nonzero source without pair events")
        return 0.0
    return error / denominator


def weighted_source_difference(reference, candidate, areas):
    """Return absolute and relative L1 including zero/tail samples.

    No arbitrary unit floor is used. A nonzero candidate against zero reference
    has no relative metric and cannot be reported as a zero relative error.
    """
    r, c, w = list(reference), list(candidate), list(areas)
    if not r or not len(r) == len(c) == len(w):
        raise ValueError("matching nonempty arrays required")
    if any(not math.isfinite(x) for x in r + c + w) or any(x < 0 for x in w):
        raise ValueError("finite samples and nonnegative support required")
    denominator = math.fsum(abs(x) * a for x, a in zip(r, w))
    error = math.fsum(abs(x - y) * a for x, y, a in zip(r, c, w))
    return {"weighted_absolute_l1": error, "reference_absolute_l1": denominator,
            "relative_l1": error / denominator if denominator else (0.0 if error == 0 else None),
            "maximum_absolute_difference": max(abs(x - y) for x, y in zip(r, c)),
            "all_samples_included": True}
