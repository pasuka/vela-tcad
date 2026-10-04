import sys,unittest,tempfile
from pathlib import Path
from decimal import Decimal
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import simplemos_srh_cloud_20260926 as s

class Audit(unittest.TestCase):
    def test_signed_obtuse_shares_conserve_area_without_clipping(self):
        mesh=dict(nodes=[dict(id=i,x=x,y=y) for i,(x,y) in enumerate(((0,0),(1,0),(.2,.1)))],triangles=[dict(node_ids=[0,1,2],region_id=0)],regions=[dict(id=0,material='Si')])
        a=s.signed_si(mesh)
        self.assertLess(min(a),0)
        self.assertAlmostEqual(sum(a)/1e-12,.05,places=13)
        mesh['regions'][0]['material']='SiO2'
        self.assertTrue((s.signed_si(mesh)==0).all())
    def test_state_response_retains_low_coordinates_below_ulp(self):
        with tempfile.TemporaryDirectory() as td:
            row=dict(node_id=0,packed_potential_scale_V=1.,electron_qf_reference_V=.5,hole_qf_reference_V=-.5)
            for k in ('psi','electron_qf_increment','hole_qf_increment'):
                row['packed_'+k]=1.;row['packed_'+k+'_low']=2**-60
            s.h.csvout(Path(td)/'state.csv',[row])
            v=s.state_vector(dict(dest=td),dict(count=1,free_si=[0]))
            self.assertGreater(v[0],Decimal(1));self.assertEqual(float(v[0]),1.)
            self.assertGreater(v[1],Decimal('1.5'))
if __name__=='__main__':unittest.main()
