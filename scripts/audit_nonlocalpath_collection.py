"""Independent geometry/unit checks for interpreted NonLocalPath observations.

This module contains no native binary-layout decoder, fitted source coefficient,
solver dispatch or claim of native window/integration identity. Callers retain
the raw observations and their version/hash provenance separately.
"""
import math
import statistics


def centroid_chord(positions, valence):
    """Independent Tri3 centroid chord along minus the P1 band gradient.

    Coordinates and energies must each use one consistent unit. This is a
    local geometry primitive, without boundary tracing or source qualification.
    The returned barycentric gradients permit deterministic face traversal.
    """
    if (len(positions)!=3 or any(len(p)!=2 for p in positions) or len(valence)!=3 or
            not all(math.isfinite(x) for p in positions for x in p) or
            not all(map(math.isfinite,valence))):
        raise ValueError("finite Tri3 coordinates and band values required")
    u=[positions[1][k]-positions[0][k] for k in (0,1)]
    v=[positions[2][k]-positions[0][k] for k in (0,1)]
    det=u[0]*v[1]-u[1]*v[0]
    if not math.isfinite(det) or det==0:
        raise ValueError("nondegenerate triangle required")
    g1=[v[1]/det,-v[0]/det];g2=[-u[1]/det,u[0]/det]
    gradients=[[-g1[k]-g2[k] for k in (0,1)],g1,g2]
    grad=[g1[k]*(valence[1]-valence[0])+g2[k]*(valence[2]-valence[0]) for k in (0,1)]
    norm=math.hypot(*grad)
    if not math.isfinite(norm) or norm==0:
        raise ValueError("nonzero finite source gradient required")
    direction=[-x/norm for x in grad]
    slopes=[math.fsum(g[k]*direction[k] for k in (0,1)) for g in gradients]
    lower=max((-1/(3*s) for s in slopes if s>0),default=-math.inf)
    upper=min((-1/(3*s) for s in slopes if s<0),default=math.inf)
    if not math.isfinite(lower) or not math.isfinite(upper) or not lower<0<upper:
        raise ValueError("unrepresentable source chord")
    return dict(direction=direction, gradient_norm=norm, area=abs(det)/2,
                centroid=[math.fsum(p[k] for p in positions)/3 for k in (0,1)],
                gradients=gradients, slopes=slopes, lower_distance=lower, upper_distance=upper)


def recursive_midpoint_samples(evaluate, lower, upper, *, maximum_width,
                               minimum_width, relative_tolerance,
                               absolute_tolerance, maximum_depth):
    """Compare signed parent/child midpoint estimates, retaining both children.

    Explicit controls define a diagnostic quadrature, not native defaults.
    Resource or minimum-width stops retain the estimate and mark it unresolved;
    callers must reject that result for actual RHS use. The scalar error test
    does not certify collection-vector errors or oscillatory integrands.
    """
    if (not all(math.isfinite(v) for v in (lower,upper,maximum_width,minimum_width,
            relative_tolerance,absolute_tolerance)) or not lower<upper or
            not math.isfinite(upper-lower) or not 0<=minimum_width<maximum_width or
            relative_tolerance<0 or absolute_tolerance<0 or
            type(maximum_depth) is not int or not 0<=maximum_depth<=40):
        raise ValueError("finite ordered window and explicit quadrature controls required")
    samples=[];result=dict(evaluations=0,width_splits=0,error_splits=0,
                           unresolved_intervals=0,maximum_depth=0)
    def value(a,b):
        v=evaluate(a,b);result["evaluations"]+=1
        if not math.isfinite(v):raise ValueError("nonfinite midpoint integral")
        return v
    def visit(a,b,parent,depth):
        mid=a+(b-a)/2
        if not a<mid<b:
            samples.append((a,b,parent));result["unresolved_intervals"]+=1;return
        left,right=value(a,mid),value(mid,b)
        error=abs((left+right)-parent)
        budget=absolute_tolerance*(b-a)/(upper-lower)+relative_tolerance*abs(parent)
        forced=b-a>maximum_width
        inaccurate=error>budget
        result["maximum_depth"]=max(result["maximum_depth"],depth)
        if forced or inaccurate:
            if depth>=maximum_depth or b-a<=minimum_width:
                result["unresolved_intervals"]+=1
            else:
                result["width_splits" if forced else "error_splits"]+=1
                visit(a,mid,left,depth+1);visit(mid,b,right,depth+1);return
        samples.extend([(a,mid,left),(mid,b,right)])
    visit(lower,upper,value(lower,upper),0)
    result.update(samples=samples,integral=math.fsum(v for a,b,v in samples),
                  resolved=result["unresolved_intervals"]==0)
    return result


def fraction_barrier_integrals(distances_m, conduction_J, valence_J, energy_J,
                               conduction_mass_kg, valence_mass_kg, hbar_J_s):
    """Audit a monotone barrier using P1 band fractions and averaged scales.

    z=(E-Ev)/(Ec-Ev) is interpolated per segment. The sqrt(2*m*gap)/hbar
    scales are averaged at that segment's original endpoints, including for
    partial turning segments. This representation differs from P1 barrier
    energy when the gap varies; it does not replace the C++ N3 kernel.
    Coincident, identical region boundary knots are retained and skipped.
    The result contains integrals only, without a source or RHS qualification.
    """
    return _fraction_barrier_integrals(distances_m,conduction_J,valence_J,energy_J,
                                      conduction_mass_kg,valence_mass_kg,hbar_J_s,False)


def piecewise_fraction_barrier_integrals(distances_m, conduction_J, valence_J, energy_J,
                                         conduction_mass_kg, valence_mass_kg, hbar_J_s):
    """Integrate each clipped P1 fraction leg, including interior reversals.

    This mathematical control keeps positive physical segment lengths. It
    does not select a physical path, qualify multiple turning points, or
    replace the strictly monotone audit contract above.
    """
    return _fraction_barrier_integrals(distances_m,conduction_J,valence_J,energy_J,
                                      conduction_mass_kg,valence_mass_kg,hbar_J_s,True)


def _fraction_barrier_integrals(distances_m,conduction_J,valence_J,energy_J,
                                conduction_mass_kg,valence_mass_kg,hbar_J_s,piecewise):
    n = len(distances_m)
    if (n < 2 or len(conduction_J) != n or len(valence_J) != n or
            not all(math.isfinite(x) for x in [*distances_m, *conduction_J, *valence_J, energy_J]) or
            not all(math.isfinite(x) and x > 0 for x in
                    [conduction_mass_kg, valence_mass_kg, hbar_J_s])):
        raise ValueError("finite matched SI knots and positive masses/hbar required")
    gaps = [c-v for c, v in zip(conduction_J, valence_J)]
    if not all(math.isfinite(g) and g > 0 for g in gaps):
        raise ValueError("positive finite band gaps required")
    z = [(energy_J-v)/g for v, g in zip(valence_J, gaps)]
    if not all(map(math.isfinite, z)) or z[0] > 0 or z[-1] < 1:
        raise ValueError("profile must bracket both turning fractions")
    cross = conduction_mass_kg/(conduction_mass_kg+valence_mass_kg)
    result = dict(actionV=0., actionC=0., inverseV=0., inverseC=0.)
    for i in range(n-1):
        span = distances_m[i+1]-distances_m[i]
        if span == 0 and conduction_J[i] == conduction_J[i+1] and valence_J[i] == valence_J[i+1]:
            continue
        if not math.isfinite(span) or span <= 0 or (not piecewise and z[i+1] <= z[i]):
            raise ValueError("ordered distance and strictly increasing fraction required")
        delta=abs(z[i+1]-z[i])
        if delta==0:
            value=z[i]
            if value in (0.,1.):raise ValueError("singular constant turning segment")
            if not 0<value<1:continue
            carrier,mass=("V",valence_mass_kg) if value<cross else ("C",conduction_mass_kg)
            fraction=value if carrier=="V"else 1-value
            k=(math.sqrt(2*mass*gaps[i])+math.sqrt(2*mass*gaps[i+1]))/(2*hbar_J_s)
            result["action"+carrier]+=span*k*math.sqrt(fraction)
            result["inverse"+carrier]+=span/k/math.sqrt(fraction)
            continue
        scale = span/delta
        zlow,zhigh=min(z[i],z[i+1]),max(z[i],z[i+1])
        for carrier, mass, low, high in [
                ("V", valence_mass_kg, max(0., zlow), min(cross, zhigh)),
                ("C", conduction_mass_kg, max(cross, zlow), min(1., zhigh))]:
            if low >= high:
                continue
            k = (math.sqrt(2*mass*gaps[i])+math.sqrt(2*mass*gaps[i+1]))/(2*hbar_J_s)
            if carrier == "C":
                low, high = 1-high, 1-low
            a, b = math.sqrt(low), math.sqrt(high)
            result["action"+carrier] += k*scale*(2/3)*(high-low)*(low+a*b+high)/(a+b)
            result["inverse"+carrier] += scale/k*2*(high-low)/(a+b)
    if not all(math.isfinite(v) and v > 0 for v in result.values()):
        raise ValueError("unrepresentable or incomplete barrier integrals")
    return result


def first_crossing_partition_weights(energies, lower, upper):
    """Collect each energy at the endpoints of its first downward crossing.

    Bisect the actual crossing segment's band energies. Retain original point
    identity, including an uphill endpoint preceding a new low; never sort the
    physical path or normalize overlapping bins. Outside the path's available
    range use its first point or first minimum. This is a diagnostic mapping,
    not qualification of the caller's physical turning-point rule.
    """
    if (not energies or not all(map(math.isfinite,energies)) or
            not math.isfinite(lower) or not math.isfinite(upper) or
            not lower<upper or not math.isfinite(upper-lower)):
        raise ValueError("finite band profile and ordered energy window required")
    width=upper-lower;weights=[0.]*len(energies)
    weights[0]=max(0.,upper-max(lower,energies[0]))/width
    record=energies[0];lowest=0
    for i,(a,b) in enumerate(zip(energies,energies[1:])):
        lo=max(lower,b);hi=min(upper,a,record);mid=a/2+b/2
        if lo<hi:
            weights[i]+=max(0.,hi-max(lo,mid))/width
            weights[i+1]+=max(0.,min(hi,mid)-lo)/width
        if b<record:record=b;lowest=i+1
    weights[lowest]+=max(0.,min(upper,record)-lower)/width
    if abs(math.fsum(weights)-1)>128*math.ulp(1.):
        raise ValueError("incomplete first-crossing energy partition")
    return weights


def bounded_energy_partition_weights(energies, lower, upper):
    """Audit midpoint energy bins with nonmonotonicity away from the window.

    Preserve every point's original identity; require the clipped intervals
    to cover the entire window exactly once. Overlaps/gaps fail explicitly.
    This is a collection control, without physical path qualification.
    """
    if (not energies or not all(map(math.isfinite,energies)) or energies[0]<energies[-1] or
            not math.isfinite(lower) or not math.isfinite(upper) or not lower<upper or
            not math.isfinite(upper-lower)):
        raise ValueError("finite band profile and ordered window required")
    mid=[a/2+b/2 for a,b in zip(energies,energies[1:])]
    intervals=[(max(lower,lo),min(upper,hi)) for lo,hi in zip(mid+[-math.inf],[math.inf]+mid)]
    cursor=lower;tolerance=128*math.ulp(1.)*max(abs(lower),abs(upper),upper-lower)
    for lo,hi in sorted((lo,hi)for lo,hi in intervals if lo<hi):
        if abs(lo-cursor)>tolerance:raise ValueError("energy bins overlap or leave a gap")
        cursor=hi
    if abs(cursor-upper)>tolerance:raise ValueError("incomplete energy bin coverage")
    return [max(0.,hi-lo)/(upper-lower)for lo,hi in intervals]


def energy_partition_weights(descending_energies, lower, upper):
    """Overlap a finite energy window with nearest-band-point intervals.

    Each interval boundary is the arithmetic midpoint of adjacent band
    energies. End intervals extend to infinity. Equal-energy adjacent points
    retain their separate, coincident support; no path point is removed.
    This is an independent collection hypothesis, not a source kernel or a
    claim that nonmonotone bands, physical turning points or RHS are qualified.
    All energies and bounds must use the same explicit energy unit.
    """
    if (not descending_energies or
            not all(math.isfinite(x) for x in descending_energies) or
            not math.isfinite(lower) or not math.isfinite(upper) or
            not lower < upper or not math.isfinite(upper-lower) or
            any(a < b for a, b in zip(descending_energies, descending_energies[1:]))):
        raise ValueError("finite descending band energies and a positive finite window required")
    midpoints = [a/2+b/2 for a, b in zip(descending_energies, descending_energies[1:])]
    high = [math.inf] + midpoints
    low = midpoints + [-math.inf]
    return [max(0., min(upper, hi)-max(lower, lo))/(upper-lower)
            for hi, lo in zip(high, low)]


def scatter_collection_coefficients(point_supports, coefficients_per_m_s):
    """Reconstruct nodal collection from signed coefficients on path points.

    This is a linear support audit, not an autonomous energy-window model.
    Each support is a sequence of (global node, barycentric weight) pairs.
    No weight renormalization, clipping or small-source floor is applied.
    """
    if len(point_supports) != len(coefficients_per_m_s):
        raise ValueError("one signed coefficient per path point required")
    terms = {}
    for support, coefficient in zip(point_supports, coefficients_per_m_s):
        if not support or not math.isfinite(coefficient):
            raise ValueError("nonempty support and finite coefficient required")
        if any(not isinstance(i, int) or i < 0 or not math.isfinite(w) or w < 0
               for i, w in support):
            raise ValueError("nonnegative finite barycentric weights and valid nodes required")
        if abs(math.fsum(w for i, w in support) - 1) > 128 * math.ulp(1.):
            raise ValueError("support weights must sum to one without normalization")
        for node, weight in support:
            terms.setdefault(node, []).append(coefficient * weight)
    nodal = {node: math.fsum(values) for node, values in sorted(terms.items())}
    total = math.fsum(coefficients_per_m_s)
    absolute = math.fsum(abs(x) for x in coefficients_per_m_s)
    difference = math.fsum(nodal.values()) - total
    return {"nodal_integrals_per_m_s": nodal, "coefficient_total_per_m_s": total,
            "absolute_coefficient_sum_per_m_s": absolute,
            "balance_difference_per_m_s": difference,
            "relative_balance": abs(difference) / absolute if absolute else 0.,
            "energy_kernel_qualified": False}


def audit_pair_centroids(positions_m, hole_integrals, electron_integrals,
                         expected_displacement_m):
    """Audit P1 collection moments using signed particle line-source integrals.

    Inputs are nodal integrals in m^-1 s^-1, not nodal rate densities. The
    displacement must come from independent bands/geometry, never a fit to
    these centroids. A zero signed total has no centroid; its raw moment and
    absolute source remain available. This does not qualify an energy kernel.
    """
    if (len(expected_displacement_m) != 2 or
            not all(math.isfinite(x) for x in expected_displacement_m) or
            any(len(p) != 2 or not all(math.isfinite(x) for x in p)
                for p in positions_m)):
        raise ValueError("finite SI positions and independent displacement required")
    blocks = {}
    for name, source in [("hole", hole_integrals), ("electron", electron_integrals)]:
        if any(not isinstance(i, int) or i < 0 or i >= len(positions_m) or
               not math.isfinite(value) for i, value in source.items()):
            raise ValueError("finite signed nodal integrals and valid nodes required")
        total = math.fsum(source.values())
        moment = [math.fsum(value * positions_m[i][k] for i, value in source.items())
                  for k in (0, 1)]
        blocks[name] = {"signed_total_per_m_s": total,
                        "absolute_nodal_sum_per_m_s": math.fsum(abs(x) for x in source.values()),
                        "first_moment_per_s": moment,
                        "centroid_m": [x / total for x in moment] if total else None}
    h, e = blocks["hole"], blocks["electron"]
    residual = [math.fsum([e["first_moment_per_s"][k], -h["first_moment_per_s"][k],
                          -expected_displacement_m[k] * h["signed_total_per_m_s"]])
                for k in (0, 1)]
    actual = ([e["centroid_m"][k] - h["centroid_m"][k] for k in (0, 1)]
              if h["centroid_m"] is not None and e["centroid_m"] is not None else None)
    return {"blocks": blocks,
            "signed_pair_difference_per_m_s": e["signed_total_per_m_s"] - h["signed_total_per_m_s"],
            "centroid_displacement_m": actual,
            "first_moment_residual_per_s": residual,
            "energy_kernel_qualified": False}


def density_normalization(physical_cm3, printed_density):
    """Identify units from populations only, never from a source/RHS fit."""
    ratios = []
    zeros = 0
    for node, value in printed_density.items():
        physical = physical_cm3[node]
        if not math.isfinite(value) or not math.isfinite(physical) or value < 0 or physical < 0:
            raise ValueError("finite nonnegative populations required")
        if value == 0 or physical == 0:
            if value != physical:
                raise ValueError("unmatched zero population")
            zeros += 1
        else:
            ratios.append(physical / value)
    if not ratios:
        raise ValueError("nonzero populations required to establish units")
    scale = statistics.median(ratios)
    spread = max(abs(value / scale - 1) for value in ratios)
    if spread > 1e-12:
        raise ValueError("populations do not define a common printed normalization")
    return {"concentration_cm3_per_unit": scale, "maximum_relative_spread": spread,
            "nonzero_samples": len(ratios), "matched_zero_samples": zeros,
            "basis": "physical population divided by printed variable"}


def bind_source_rhs(on, off, physical_integrals_per_cm_s, concentration_cm3):
    """Check +G in native RHS1 after an exact printed full-state identity gate.

    physical_integrals_per_cm_s maps block name -> global node -> G integral
    in cm^-1 s^-1. No q, arbitrary unit floor or inferred source factor is used.
    The diagnostic denominator is the absolute source sum, NOT the frozen full
    RHS denominator of all assembly terms. Small differences can therefore be
    limited by cancellation in the native assembler or printed RHS values.
    """
    if (on["names"] != off["names"] or not on["variables"] or
            on["variables"] != off["variables"]):
        raise ValueError("complete matching printed variables required")
    keys = set(on["variables"])
    if set(on["rhs1"]) != keys or set(off["rhs1"]) != keys:
        raise ValueError("incomplete native RHS mapping")
    if not math.isfinite(concentration_cm3) or concentration_cm3 <= 0:
        raise ValueError("explicit positive concentration normalization required")
    result = {}
    selected = set()
    for name, integrals in physical_integrals_per_cm_s.items():
        block = on["names"].index(name)
        rows = []
        for key in sorted(keys):
            if key[0] != block:
                continue
            expected = integrals[key[1]] / concentration_cm3
            observed = on["rhs1"][key] - off["rhs1"][key]
            if not all(map(math.isfinite, (expected, observed, on["rhs1"][key], off["rhs1"][key]))):
                raise ValueError("nonfinite source/RHS value")
            rows.append({"node": key[1], "expected": expected, "observed": observed,
                         "absolute_error": abs(observed - expected)})
        if not rows:
            raise ValueError("empty carrier block")
        selected.add(block)
        denominator = math.fsum(abs(row["expected"]) for row in rows)
        error = math.fsum(row["absolute_error"] for row in rows)
        result[name] = {"rows": rows, "absolute_source_sum": denominator,
                        "absolute_error_sum": error,
                        "source_only_relative_L1": error / denominator if denominator else (0 if error == 0 else None)}
    outside = [(key, on["rhs1"][key] - off["rhs1"][key]) for key in sorted(keys) if key[0] not in selected]
    return {"blocks": result, "full_printed_state_identical": True,
            "outside_carrier_changed_rows": sum(value != 0 for key, value in outside),
            "outside_carrier_absolute_change": math.fsum(abs(value) for key, value in outside),
            "full_rhs_qualified": False}


def audit_centroid_path(positions_m, triangles, valence_J, path, length_scale_m, thermal_voltage_V, charge_C):
    """Validate a decoded ray against independently supplied mesh/bands/units.

    Each point provides cell, local_vertices, weights and native arclength t.
    Native ABI/addresses are deliberately outside this mathematical contract.
    A large returned mismatch is evidence against a hypothesis, not an excuse
    to refit the supplied units. Reflected rays need separate straightness use.
    """
    if (len(positions_m) != len(valence_J) or
            not all(math.isfinite(x) for p in positions_m for x in p) or
            not all(math.isfinite(x) for x in valence_J) or
            any(len(p) != 2 for p in positions_m) or
            not all(math.isfinite(x) and x > 0 for x in (length_scale_m, thermal_voltage_V, charge_C))):
        raise ValueError("finite SI mesh/bands and positive explicit units required")
    seed = path["seed_cell"]
    vertices = triangles[seed]
    if len(vertices) != 3 or len(set(vertices)) != 3 or min(vertices) < 0 or max(vertices) >= len(positions_m):
        raise ValueError("valid Tri3 seed required")
    p0, p1, p2 = [positions_m[i] for i in vertices]
    u = [p1[k] - p0[k] for k in (0, 1)]
    v = [p2[k] - p0[k] for k in (0, 1)]
    det = u[0] * v[1] - u[1] * v[0]
    if det == 0:
        raise ValueError("degenerate source triangle")
    a, b = valence_J[vertices[1]] - valence_J[vertices[0]], valence_J[vertices[2]] - valence_J[vertices[0]]
    gradient = [(a * v[1] - b * u[1]) / det, (u[0] * b - v[0] * a) / det]
    expected_ray = [-x / charge_C * length_scale_m / thermal_voltage_V for x in gradient]
    ray = path["ray_gradient"][:2]
    norm = math.hypot(*ray)
    if norm == 0 or not all(map(math.isfinite, ray)):
        raise ValueError("nonzero finite ray required")
    direction = [x / norm for x in ray]
    centroid = [math.fsum(positions_m[i][k] for i in vertices) / 3 for k in (0, 1)]
    h = max(math.dist(positions_m[i], positions_m[j]) for i in vertices for j in vertices)
    points = []
    for point in path["points"]:
        cell = triangles[point["cell"]]
        local, weights = point["local_vertices"], point["weights"]
        if (len(local) != len(weights) or not weights or
                not all(math.isfinite(w) and w >= 0 for w in weights) or
                abs(math.fsum(weights) - 1) > 128 * math.ulp(1.0) or
                any(i < 0 or i >= len(cell) for i in local) or
                not math.isfinite(point["t"])):
            raise ValueError("invalid interpreted point support")
        xy = [math.fsum(w * positions_m[cell[i]][k] for i, w in zip(local, weights)) for k in (0, 1)]
        points.append((point["t"], xy, point["cell"]))
    if (len(points) < 2 or points[0][2] != seed or points[1][2] != seed or
            not points[0][0] < 0 < points[1][0]):
        raise ValueError("source faces must bracket the origin in the seed cell")
    duplicate_boundaries = 0
    coordinate_scale = max(h, *(abs(x) for _, xy, _ in points for x in xy))
    for previous, current in zip(points, points[1:]):
        if current[0] < previous[0]:
            raise ValueError("path arclength must not decrease")
        if current[0] == previous[0]:
            # Two regions may each retain the same physical boundary point.
            # Accept its duplicate representation, never a zero-length jump.
            if (current[2] == previous[2] or
                    math.dist(current[1], previous[1]) > 128 * math.ulp(coordinate_scale)):
                raise ValueError("equal arclength requires the same point in distinct cells")
            duplicate_boundaries += 1
    t0, t1 = points[0][0], points[1][0]
    alpha = -t0 / (t1 - t0)
    anchor = [(1-alpha) * points[0][1][k] + alpha * points[1][1][k] for k in (0, 1)]
    return {"points": len(points), "duplicate_boundary_points": duplicate_boundaries,
            "centroid_error_over_cell_diameter": math.dist(anchor, centroid) / h,
            "ray_gradient_relative_error": math.hypot(*(ray[k]-expected_ray[k] for k in (0, 1))) / norm,
            "straightness_error_over_cell_diameter": max(math.dist(xy, [centroid[k]+t*length_scale_m*direction[k] for k in (0, 1)]) / h for t, xy, cell in points),
            "arclength_error_over_cell_diameter": max(abs(math.dist(points[i][1],points[i-1][1])-(points[i][0]-points[i-1][0])*length_scale_m) / h for i in range(1,len(points))),
            "native_window_mapping_qualified": False}
