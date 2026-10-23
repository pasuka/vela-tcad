import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from audit_sentaurus_density_updates import audit_update, independent_rhs, vector_difference


class UpdateAuditTests(unittest.TestCase):
    def test_negative_update_is_local_and_positive_small_values_survive(self):
        old={1:1e-30,2:1e-30,3:1e-30};delta={1:-2e-30,2:-.999e-30,3:2e-30}
        actual={1:1e-32,2:old[2]+delta[2],3:3e-30}
        r=audit_update(old,delta,actual)
        self.assertEqual(r['nonpositive_raw_nodes'],[1])
        self.assertEqual(r['positive_below_one_percent_count'],1)
        self.assertLess(r['max_reconstruction_relative_error'],1e-15)

    def test_small_positive_update_is_not_a_one_percent_floor(self):
        r=audit_update({1:1.0},{1:-.999},{1:.01})
        self.assertGreater(r['max_reconstruction_relative_error'],8)

    def test_reject_bad_mapping_and_nonpositive_old(self):
        with self.assertRaises(ValueError):audit_update({1:1},{2:1},{1:2})
        with self.assertRaises(ValueError):audit_update({1:0},{1:1},{1:1})

    def test_contact_rhs_is_counted_once(self):
        values=independent_rhs({(0,1):4.,(1,233):3.,(1,234):3.})
        self.assertEqual(sum(v*v for v in values),25.)
        with self.assertRaises(ValueError):independent_rhs({(1,233):3.,(1,234):4.})

    def test_rhs_difference_has_no_unit_floor(self):
        self.assertAlmostEqual(vector_difference({(3,1):1e-30},{(3,1):2e-30}),1/3)


if __name__=='__main__':unittest.main()
