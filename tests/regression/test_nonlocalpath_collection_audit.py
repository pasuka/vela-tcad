import copy
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
import audit_nonlocalpath_collection as collection_audit
from audit_nonlocalpath_collection import (density_normalization, bind_source_rhs,
                                          audit_centroid_path, audit_pair_centroids,
                                          scatter_collection_coefficients)


class CollectionAuditTests(unittest.TestCase):
    def test_first_crossing_collection_retains_real_segment_endpoint_identity(self):
        self.assertTrue(hasattr(collection_audit,"first_crossing_partition_weights"))
        r=collection_audit.first_crossing_partition_weights([9.,5.,7.,1.],0.,9.)
        self.assertEqual(r,[2/9,2/9,1/9,4/9])
        # The rising node at energy 7 still receives the upper half of the
        # subsequent first-crossing segment. Sorting or deleting it is wrong.
        self.assertGreater(r[2],0.)

    def test_first_crossing_reduces_to_monotone_bins_and_preserves_gauge_units(self):
        fn=collection_audit.first_crossing_partition_weights
        for band in [[9.,5.,1.],[9.,5.,5.,1.],[5.]]:
            for a,b in [(2.,8.),(-3.,12.),(4.,6.)]:
                actual=fn(band,a,b);expected=collection_audit.energy_partition_weights(band,a,b)
                for x,y in zip(actual,expected):self.assertAlmostEqual(x,y)
                self.assertAlmostEqual(math.fsum(actual),1.)
        self.assertEqual(fn([25.,21.,23.,17.],16.,25.),fn([9.,5.,7.,1.],0.,9.))
        self.assertEqual(fn([27.,15.,21.,3.],0.,27.),fn([9.,5.,7.,1.],0.,9.))

    def test_piecewise_fraction_integrates_reversed_and_constant_interior_legs(self):
        self.assertTrue(hasattr(collection_audit,"piecewise_fraction_barrier_integrals"))
        fn=collection_audit.piecewise_fraction_barrier_integrals
        # Constant gap=4; fractions 0, .3, .2, .2, 1 include an interior
        # reversal and a plateau. Every segment contributes positive length.
        z=[0.,.3,.2,.2,1.];ev=[-4*x for x in z];ec=[v+4 for v in ev]
        r=fn([0.,1.,2.,3.,4.],ec,ev,0.,1.,1.,1.)
        k=math.sqrt(8);cross=.5
        expected={key:0. for key in r}
        for a,b in zip(z,z[1:]):
            for carrier in ['V','C']:
                lo,hi=sorted((a,b));lo=max(lo,0. if carrier=='V'else cross);hi=min(hi,cross if carrier=='V'else 1.)
                if lo>=hi and a!=b:continue
                if a==b:
                    if (a<cross)!=(carrier=='V'):continue
                    v=a if carrier=='V'else 1-a
                    action=k*math.sqrt(v);inverse=1/(k*math.sqrt(v))
                else:
                    if carrier=='C':lo,hi=1-hi,1-lo
                    action=k*2/3*(hi**1.5-lo**1.5)/abs(b-a)
                    inverse=2/k*(math.sqrt(hi)-math.sqrt(lo))/abs(b-a)
                expected['action'+carrier]+=action;expected['inverse'+carrier]+=inverse
        for key in r:self.assertAlmostEqual(r[key],expected[key])
        ordinary=collection_audit.fraction_barrier_integrals([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.)
        self.assertEqual(fn([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.),ordinary)

    def test_nonmonotone_energy_bins_require_a_complete_nonoverlapping_window(self):
        self.assertTrue(hasattr(collection_audit,"bounded_energy_partition_weights"))
        fn=collection_audit.bounded_energy_partition_weights
        self.assertEqual(fn([12.,10.,11.,8.,4.,0.],1.,3.),[0.,0.,0.,0.,.5,.5])
        for energies in [[9.,1.,12.,0.],[0.,1.]]:
            with self.assertRaises(ValueError):fn(energies,2.,6.)

    def test_recursive_midpoint_keeps_children_and_signed_integrals(self):
        self.assertTrue(hasattr(collection_audit, "recursive_midpoint_samples"))
        fn = collection_audit.recursive_midpoint_samples
        for sign in [1., -1.]:
            r = fn(lambda a,b: sign*(b-a)*((a+b)/2), 0., 1.,
                   maximum_width=1., minimum_width=0., relative_tolerance=1e-6,
                   absolute_tolerance=0., maximum_depth=12)
            self.assertEqual(len(r["samples"]), 2)
            self.assertEqual(r["integral"], sign*.5)
            self.assertTrue(r["resolved"])

    def test_recursive_midpoint_resolves_curvature_and_forced_width(self):
        fn = collection_audit.recursive_midpoint_samples
        r = fn(lambda a,b: (b-a)*math.exp((a+b)/2), 0., 1.,
               maximum_width=.2, minimum_width=0., relative_tolerance=1e-5,
               absolute_tolerance=0., maximum_depth=20)
        self.assertLess(abs(r["integral"]-math.expm1(1.)), 1e-5)
        self.assertGreater(r["error_splits"], 0)
        self.assertGreater(r["width_splits"], 0)
        self.assertEqual(math.fsum(b-a for a,b,v in r["samples"]),1.)

    def test_recursive_midpoint_never_calls_resource_failure_a_zero(self):
        fn = collection_audit.recursive_midpoint_samples
        r = fn(lambda a,b: (b-a)*((a+b)/2)**2, 0., 1.,
               maximum_width=.01, minimum_width=0., relative_tolerance=1e-12,
               absolute_tolerance=0., maximum_depth=1)
        self.assertFalse(r["resolved"])
        self.assertGreater(r["unresolved_intervals"], 0)
        self.assertGreater(r["integral"], 0.)
        with self.assertRaises(ValueError):
            fn(lambda a,b: math.nan, 0., 1., maximum_width=1., minimum_width=0.,
               relative_tolerance=1e-5, absolute_tolerance=0., maximum_depth=10)

    def test_centroid_chord_is_geometric_and_rotation_covariant(self):
        self.assertTrue(hasattr(collection_audit, "centroid_chord"))
        fn = collection_audit.centroid_chord
        r=fn([[0.,0.],[3.,0.],[0.,3.]], [0.,-3.,0.])
        self.assertEqual(r["direction"],[1.,0.])
        self.assertAlmostEqual(r["lower_distance"],-1.)
        self.assertAlmostEqual(r["upper_distance"],1.)
        self.assertAlmostEqual(r["area"],4.5)
        q=fn([[7.,-4.],[7.,-1.],[4.,-4.]], [13.,10.,13.])
        self.assertEqual(q["direction"],[0.,1.])
        self.assertAlmostEqual(q["lower_distance"],r["lower_distance"])
        self.assertAlmostEqual(q["upper_distance"],r["upper_distance"])
        with self.assertRaises(ValueError): fn([[0.,0.],[1.,0.],[2.,0.]],[0.,1.,2.])

    def test_fraction_barrier_constant_gap_recovers_parabolic_legs(self):
        self.assertTrue(hasattr(collection_audit, "fraction_barrier_integrals"),
                        "transformed-coordinate barrier audit is missing")
        r = collection_audit.fraction_barrier_integrals([0.,2.],[5.,-3.],[1.,-7.],0.,2.,3.,1.)
        z = 2/5
        self.assertAlmostEqual(r["actionV"], math.sqrt(24)*2/3*z**1.5)
        self.assertAlmostEqual(r["actionC"], 4*2/3*(1-z)**1.5)
        self.assertAlmostEqual(r["inverseV"], 2*math.sqrt(z)/math.sqrt(24))
        self.assertAlmostEqual(r["inverseC"], 2*math.sqrt(1-z)/4)

    def test_fraction_barrier_variable_gap_matches_independent_regularized_quadrature(self):
        r = collection_audit.fraction_barrier_integrals([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.)
        # Direct quadrature after z=u^2 removes the reciprocal singularity.
        # Here z is the explicitly linear fraction, rather than the rational
        # fraction obtained by interpolating the original energies first.
        def simpson(f, upper):
            n = 40; h = upper/n
            return h/3*math.fsum(f(i*h)*(1 if i in [0,n] else 4 if i%2 else 2)
                                for i in range(n+1))
        kv = (math.sqrt(8)+math.sqrt(18))/2; kc = 2*kv
        self.assertAlmostEqual(r["actionV"], simpson(lambda u:2*kv*u*u,math.sqrt(.8)))
        self.assertAlmostEqual(r["inverseV"], simpson(lambda u:2/kv,math.sqrt(.8)))
        self.assertAlmostEqual(r["actionC"], simpson(lambda u:2*kc*u*u,math.sqrt(.2)))
        self.assertAlmostEqual(r["inverseC"], simpson(lambda u:2/kc,math.sqrt(.2)))

    def test_fraction_barrier_preserves_energy_gauge(self):
        fn = collection_audit.fraction_barrier_integrals
        self.assertEqual(fn([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.),
                         fn([0.,1.],[20.,16.],[16.,7.],16.,4.,1.,1.))

    def test_fraction_barrier_units_scale_action_and_inverse_separately(self):
        fn = collection_audit.fraction_barrier_integrals
        base = fn([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.)
        for result, sa, si in [(fn([3.,5.],[4.,0.],[0.,-9.],0.,4.,1.,1.),2.,2.),
                               (fn([0.,1.],[16.,0.],[0.,-36.],0.,4.,1.,1.),2.,.5),
                               (fn([0.,1.],[4.,0.],[0.,-9.],0.,16.,4.,1.),2.,.5),
                               (fn([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,2.),.5,2.)]:
            for key in base:
                self.assertAlmostEqual(result[key]/base[key],sa if key.startswith("action") else si)

    def test_fraction_barrier_identical_region_boundary_duplicate_is_retained(self):
        fn = collection_audit.fraction_barrier_integrals
        self.assertEqual(fn([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,1.),
                         fn([0.,0.,1.],[4.,4.,0.],[0.,0.,-9.],0.,4.,1.,1.))

    def test_fraction_barrier_rejects_unresolved_or_nonmonotone_profiles(self):
        fn = collection_audit.fraction_barrier_integrals
        bad = [([],[],[],0.,4.,1.,1.),([0.,1.],[4.,0.],[0.,-9.],0.,0.,1.,1.),
               ([0.,1.],[4.,0.],[0.,-9.],0.,4.,1.,math.nan),
               ([0.,1.],[4.,0.],[0.,-9.],10.,4.,1.,1.),
               ([0.,1.],[4.,-9.],[0.,-9.],0.,4.,1.,1.),
               ([0.,0.,1.],[4.,5.,0.],[0.,0.,-9.],0.,4.,1.,1.),
               ([0.,.5,1.],[4.,4.,0.],[0.,1.,-9.],0.,4.,1.,1.)]
        for args in bad:
            with self.assertRaises(ValueError): fn(*args)

    def test_energy_window_inside_one_nearest_band_bin_collects_at_one_point(self):
        self.assertTrue(hasattr(collection_audit, "energy_partition_weights"),
                        "independent energy-bin overlap audit is missing")
        weights = collection_audit.energy_partition_weights([9., 5., 1.], 4., 6.)
        self.assertEqual(weights, [0., 1., 0.])

    def test_energy_bin_crossing_splits_by_overlap_not_endpoint_interpolation(self):
        self.assertEqual(collection_audit.energy_partition_weights([9.,5.,1.],2.,8.),
                         [1/6, 4/6, 1/6])
        self.assertEqual(collection_audit.energy_partition_weights([9.,5.,1.],6.,8.),
                         [.5, .5, 0.])

    def test_energy_bin_weights_preserve_gauge_and_unit_scaling(self):
        expected = collection_audit.energy_partition_weights([9.,5.,1.],2.,8.)
        self.assertEqual(collection_audit.energy_partition_weights([27.,15.,3.],6.,24.), expected)
        self.assertEqual(collection_audit.energy_partition_weights([25.,21.,17.],18.,24.), expected)

    def test_energy_bin_endpoints_and_repeated_band_points_retain_support(self):
        self.assertEqual(collection_audit.energy_partition_weights([9.,5.,5.,1.],4.,6.),
                         [0., .5, .5, 0.])
        self.assertEqual(collection_audit.energy_partition_weights([5.],-20.,30.),[1.])
        self.assertEqual(collection_audit.energy_partition_weights([9.,5.,1.],10.,11.),[1.,0.,0.])

    def test_energy_bin_weights_do_not_clip_or_reorder_invalid_inputs(self):
        for bands, low, high in [([],0.,1.),([1.,2.],0.,1.),([math.nan],0.,1.),
                                 ([1.],1.,1.),([1.],2.,1.),([1.],0.,math.inf)]:
            with self.assertRaises(ValueError):
                collection_audit.energy_partition_weights(bands,low,high)

    def test_energy_bin_signed_collection_keeps_pair_totals_and_spatial_moments(self):
        weights = collection_audit.energy_partition_weights([9.,5.,1.],2.,8.)
        for rate in [12., -12., 1e-200]:
            result = scatter_collection_coefficients([[(i,1.)] for i in range(3)],
                                                     [rate*w for w in weights])
            self.assertAlmostEqual(result["coefficient_total_per_m_s"]/rate, 1.)
            moment = audit_pair_centroids([[0.,0.],[1.,0.],[2.,0.]],
                                         result["nodal_integrals_per_m_s"],{},[0.,0.])
            self.assertAlmostEqual(moment["blocks"]["hole"]["centroid_m"][0], 1.)

    def test_path_coefficients_scatter_to_five_nodes_without_fixed_thirds(self):
        supports = [[(0,.25),(1,.75)],[(1,.5),(2,.5)],[(3,.2),(4,.8)]]
        result = scatter_collection_coefficients(supports,[4.,2.,5.])
        self.assertEqual(result["nodal_integrals_per_m_s"],{0:1.,1:4.,2:1.,3:1.,4:4.})
        self.assertEqual(result["coefficient_total_per_m_s"],11.)
        self.assertEqual(result["relative_balance"],0.)
        self.assertFalse(result["energy_kernel_qualified"])

    def test_signed_coefficients_keep_cancellation_and_small_values(self):
        supports = [[(0,1.)],[(1,1.)]]
        result = scatter_collection_coefficients(supports,[1.,-1.])
        self.assertEqual(result["absolute_coefficient_sum_per_m_s"],2.)
        self.assertEqual(result["coefficient_total_per_m_s"],0.)
        result = scatter_collection_coefficients(supports,[1e-200,-1e-200])
        self.assertEqual(result["nodal_integrals_per_m_s"],{0:1e-200,1:-1e-200})

    def test_collection_moments_preserve_interpolated_coordinates(self):
        positions = [[0.,0.],[2.,0.],[0.,2.]]
        result = scatter_collection_coefficients([[(0,.25),(1,.25),(2,.5)]],[4.])
        audit = audit_pair_centroids(positions,result["nodal_integrals_per_m_s"],{},[0.,0.])
        self.assertEqual(audit["blocks"]["hole"]["centroid_m"],[.5,1.])

    def test_collection_does_not_repair_invalid_supports(self):
        for supports, values in [([[(0,.4),(1,.4)]],[1.]),
                                 ([[(0,1.)]],[]),([[]],[0.]),
                                 ([[(-1,1.)]],[1.]),([[(0,1.)]],[math.nan]),
                                 ([[(0,1.1),(1,-.1)]],[1.])]:
            with self.assertRaises(ValueError): scatter_collection_coefficients(supports,values)

    def test_pair_moments_preserve_signed_generation_and_recombination(self):
        p = [[0.,0.],[2.,0.],[4.,0.],[6.,0.]]
        for sign in [1., -1., 1e-100]:
            result = audit_pair_centroids(p,{0:sign,1:sign},{2:sign,3:sign},[4.,0.])
            self.assertEqual(result["centroid_displacement_m"],[4.,0.])
            self.assertEqual(result["first_moment_residual_per_s"],[0.,0.])
            self.assertEqual(result["signed_pair_difference_per_m_s"],0.)
            self.assertFalse(result["energy_kernel_qualified"])

    def test_pair_moment_rotation_translation_and_unit_scaling(self):
        p = [[7.,-4.],[7.,-2.],[7.,0.],[7.,2.]]
        result = audit_pair_centroids(p,{0:1.,1:1.},{2:1.,3:1.},[0.,4.])
        self.assertEqual(result["centroid_displacement_m"],[0.,4.])
        scaled = audit_pair_centroids([[x*2,y*2] for x,y in p],{0:3.,1:3.},
                                     {2:3.,3:3.},[0.,8.])
        self.assertEqual(scaled["first_moment_residual_per_s"],[0.,0.])

    def test_cancelled_source_retains_moment_without_centroid(self):
        result = audit_pair_centroids([[0.,0.],[1.,0.]],{0:1.,1:-1.},{},[0.,0.])
        self.assertIsNone(result["blocks"]["hole"]["centroid_m"])
        self.assertEqual(result["blocks"]["hole"]["absolute_nodal_sum_per_m_s"],2.)
        self.assertEqual(result["first_moment_residual_per_s"],[1.,0.])
        self.assertIsNone(result["centroid_displacement_m"])

    def test_wrong_endpoint_shift_is_exposed_and_bad_moments_rejected(self):
        result = audit_pair_centroids([[0.,0.],[2.,0.]],{0:1.},{1:1.},[1.,0.])
        self.assertEqual(result["first_moment_residual_per_s"],[1.,0.])
        for p,h,e,d in [([[0.,0.]],{1:1.},{},[0.,0.]),
                        ([[0.,0.]],{0:math.inf},{},[0.,0.]),
                        ([[0.,math.nan]],{},{},[0.,0.])]:
            with self.assertRaises(ValueError): audit_pair_centroids(p,h,e,d)

    def test_units_from_population_not_source(self):
        result = density_normalization({0: 2e12, 1: 4e12, 2: 0}, {0: 2., 1: 4., 2: 0.})
        self.assertEqual(result["concentration_cm3_per_unit"], 1e12)
        self.assertEqual(result["matched_zero_samples"], 1)
        self.assertEqual(result["maximum_relative_spread"], 0)

    def test_inconsistent_or_missing_population_units_rejected(self):
        for a, b in [({0: 0}, {0: 1}), ({0: 1}, {0: 0}), ({0: 0}, {0: 0}),
                     ({0: 1, 1: 2}, {0: 1, 1: 1}), ({0: math.inf}, {0: 1})]:
            with self.assertRaises(ValueError): density_normalization(a, b)

    def debug(self):
        return {"names": ["psi", "eDensity", "load"],
                "variables": {(0, 0): 0., (1, 0): 2., (2, 0): 7.},
                "rhs1": {(0, 0): 0., (1, 0): 0., (2, 0): 0.}}

    def test_signed_particle_line_source_maps_to_rhs_without_q(self):
        off = self.debug()
        for rate in [10., -10., 1e-100]:
            on = copy.deepcopy(off);on["rhs1"][1, 0] = rate / 1e12
            result = bind_source_rhs(on, off, {"eDensity": {0: rate}}, 1e12)
            self.assertEqual(result["blocks"]["eDensity"]["source_only_relative_L1"], 0)
            self.assertEqual(result["outside_carrier_changed_rows"], 0)
            self.assertFalse(result["full_rhs_qualified"])

    def test_full_state_includes_noncarrier_variables(self):
        off = self.debug();on = copy.deepcopy(off);on["variables"][2, 0] += 1
        with self.assertRaises(ValueError): bind_source_rhs(on, off, {"eDensity": {0: 0}}, 1e12)
        on = copy.deepcopy(off);del on["rhs1"][0, 0]
        with self.assertRaises(ValueError): bind_source_rhs(on, off, {"eDensity": {0: 0}}, 1e12)

    def test_zero_source_cannot_hide_a_nonzero_rhs(self):
        off = self.debug();on = copy.deepcopy(off);on["rhs1"][1, 0] = 1e-100;on["rhs1"][2, 0] = 3
        result = bind_source_rhs(on, off, {"eDensity": {0: 0}}, 1e12)
        self.assertIsNone(result["blocks"]["eDensity"]["source_only_relative_L1"])
        self.assertEqual(result["outside_carrier_absolute_change"], 3)

    def geometry(self):
        return ([[0.,0.],[3.,0.],[0.,3.]], [[0,1,2]], [0.,-3.,0.],
                {"seed_cell": 0, "ray_gradient": [1.,0.], "points": [
                    {"cell":0,"local_vertices":[0,2],"weights":[2/3,1/3],"t":-1.},
                    {"cell":0,"local_vertices":[1,2],"weights":[2/3,1/3],"t":1.}]})

    def test_ray_anchor_is_independently_the_triangle_centroid(self):
        result = audit_centroid_path(*self.geometry(),1.,1.,1.)
        for key in ["centroid_error_over_cell_diameter","ray_gradient_relative_error",
                    "straightness_error_over_cell_diameter","arclength_error_over_cell_diameter"]:
            self.assertLess(result[key], 1e-15)
        self.assertFalse(result["native_window_mapping_qualified"])

    def test_rotation_and_energy_gauge_preserve_geometry(self):
        points, cells, energies, path = self.geometry()
        points = [[-p[1]+7, p[0]-4] for p in points]
        energies = [e+13 for e in energies];path["ray_gradient"] = [0.,1.]
        result = audit_centroid_path(points,cells,energies,path,1.,1.,1.)
        self.assertLess(result["centroid_error_over_cell_diameter"],1e-15)
        self.assertLess(result["ray_gradient_relative_error"],1e-15)

    def test_wrong_units_are_reported_not_fitted(self):
        result = audit_centroid_path(*self.geometry(),2.,1.,1.)
        self.assertGreater(result["straightness_error_over_cell_diameter"],.1)
        self.assertEqual(result["ray_gradient_relative_error"],1)

    def test_region_boundary_duplicate_preserves_physical_path(self):
        points, cells, energies, path = self.geometry()
        cells.append(cells[0][:])
        duplicate = copy.deepcopy(path["points"][-1]);duplicate["cell"] = 1
        path["points"].append(duplicate)
        result = audit_centroid_path(points,cells,energies,path,1.,1.,1.)
        self.assertEqual(result["duplicate_boundary_points"],1)
        self.assertLess(result["arclength_error_over_cell_diameter"],1e-15)
        for mutation in [lambda p:p.update(cell=0),
                         lambda p:p.update(weights=[.5,.5]),
                         lambda p:p.update(t=.5)]:
            bad = copy.deepcopy(path);mutation(bad["points"][-1])
            with self.assertRaises(ValueError):
                audit_centroid_path(points,cells,energies,bad,1.,1.,1.)

    def test_invalid_support_and_origin_are_rejected(self):
        for mutation in [lambda p:p["points"][0].update(weights=[.5,.4]),
                         lambda p:p["points"][0].update(local_vertices=[0,3]),
                         lambda p:p["points"][0].update(t=0.)]:
            args = self.geometry();mutation(args[3])
            with self.assertRaises(ValueError): audit_centroid_path(*args,1.,1.,1.)


if __name__ == "__main__": unittest.main()
