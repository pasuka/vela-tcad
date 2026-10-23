"""Read-only fixed-state support checks; no solver or model defaults."""
import math
import numpy as np


def edge_divergence(node_count, pairs, flux):
    """Conservative incidence sum with explicit endpoint orientation."""
    if node_count <= 0 or len(pairs) != len(flux):
        raise ValueError('invalid incidence dimensions')
    result = np.zeros(node_count)
    for (i, j), value in zip(pairs, flux):
        if i == j or min(i, j) < 0 or max(i, j) >= node_count or not math.isfinite(value):
            raise ValueError('invalid edge or flux')
        result[i] += value
        result[j] -= value
    return result


def original_vertex_values(rows, coordinates, native_coordinates):
    """Map original mesh IDs, retaining and counting contact-copy evidence."""
    values = {}
    for index, value in rows:
        if index in values or not math.isfinite(value):
            raise ValueError('duplicate or nonfinite native scalar')
        values[index] = value
    for i, xy in enumerate(coordinates):
        if i not in values or i not in native_coordinates:
            raise ValueError('missing original vertex')
        if not np.allclose(xy, native_coordinates[i], rtol=0, atol=1e-13):
            raise ValueError('coordinate mapping mismatch')
    return np.array([values[i] for i in range(len(coordinates))]), len(values)-len(coordinates)


def element_box_mobility(nodes, mobility, measures):
    """A hypothesis using an independently identified cell's vertex measures."""
    if len(nodes) != 3 or len(set(nodes)) != 3 or len(mobility) != 3:
        raise ValueError('three distinct vertices required')
    weights = [measures[i] for i in nodes]
    if any(not math.isfinite(v) or v < 0 for v in [*weights, *mobility]) or sum(weights) <= 0:
        raise ValueError('invalid mobility or measure')
    return math.fsum(v*w for v, w in zip(mobility, weights))/math.fsum(weights)


def e2_cell_integrals(points_um, psi, measures_um2, A_SI=3.4e23, B_SI=2.26e9):
    """Original E2 constants, linear cell field, explicit source quadrature."""
    points = np.asarray(points_um, float)
    if points.shape != (3, 2) or len(psi) != 3 or len(measures_um2) != 3:
        raise ValueError('triangle and three source measures required')
    operands = np.r_[points.ravel(), psi, measures_um2, A_SI, B_SI]
    if not np.all(np.isfinite(operands)) or min(measures_um2) < 0 or A_SI < 0 or B_SI <= 0:
        raise ValueError('invalid source data')
    vectors = (points[1:]-points[0])*1e-6
    if np.linalg.det(vectors) == 0:
        raise ValueError('degenerate cell')
    field = np.linalg.norm(np.linalg.solve(vectors, np.asarray(psi)[1:]-psi[0]))
    rate_cm3_s = A_SI*field*field*math.exp(-B_SI/field)/1e6 if field else 0.
    return rate_cm3_s*np.asarray(measures_um2)


def e2_nodal_integrals(fields_V_cm, measures_um2, A_cm, B_V_cm):
    """Explicit-parameter E2 integral and vector-field partials, fixed boxes.

    Field recovery and its potential derivatives are deliberately not supplied
    here. A caller must identify that operator independently.
    """
    E=np.asarray(fields_V_cm,float);V=np.asarray(measures_um2,float)
    if E.ndim!=2 or E.shape[1]!=2 or V.shape!=(len(E),):
        raise ValueError('two-component fields and matching measures required')
    if not np.all(np.isfinite(E)) or not np.all(np.isfinite(V)) or np.any(V<0) or not math.isfinite(A_cm+B_V_cm) or A_cm<0 or B_V_cm<=0:
        raise ValueError('invalid E2 operands')
    F=np.linalg.norm(E,axis=1);rate=np.zeros(len(E));J=np.zeros_like(E);active=F>0
    rate[active]=A_cm*F[active]**2*np.exp(-B_V_cm/F[active])
    J[active]=((rate[active]*(2/F[active]+B_V_cm/F[active]**2)/F[active])*V[active])[:,None]*E[active]
    return rate*V,J
