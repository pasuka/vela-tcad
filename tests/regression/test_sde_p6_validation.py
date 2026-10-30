"""Independent cut, source-sizing and TDX stored-field qualification checks."""
from pathlib import Path
import tempfile,sys,unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from sde_native_metrics import cut_roots,cuts,multibox_spacing,source_point,sizing
from sde_tdx_validation import validate_pair


class GeometryMetricsTests(unittest.TestCase):
    def test_cut_interpolates_linear_zero_and_empty_cut_cannot_pass(self):
        model={'vertices':[(0.,0.),(1.,0.),(0.,1.)],'cells':[('R.Silicon',[0,1,2])]}
        field=lambda m,r,n:(0.,0.,m['vertices'][n][1]-.25)
        self.assertEqual([.25],cut_roots(model,'x',.2,'R.Silicon',field))
        case={'contract':{'doping':{'kind':'bjt'}}}
        empty={'vertices':[],'cells':[]}
        result=cuts(case,{'mesh_control':{}},empty,empty,field)
        self.assertEqual(3,len(result));self.assertTrue(all(not r['passed'] for r in result))

    def test_signed_multibox_grades_from_the_authored_side(self):
        size={'ratio_x':-2.,'ratio_y':1.,'min_x':.01,'max_x':.2,'min_y':.02,'max_y':.3}
        window={'lower_left':[0.,0.],'upper_right':[1.,1.]}
        self.assertEqual([.01,.3],multibox_spacing(size,window,1.,.5))
        self.assertAlmostEqual(.11,multibox_spacing(size,window,.9,.5)[0])
        self.assertEqual((.2,.3),source_point({'mesh_control':{'operations':[{'axis':'x'}]}},-.2,.3))

    def test_coarse_local_spacing_is_rejected(self):
        ir={'geometry':{'bodies':[]},'mesh_control':{'refinement_sizes':[{'name':'S','max_x':.1,'max_y':.1}],
            'refinement_placements':[{'name':'P','refinement':'S','target_kind':'region','target':'R'}]}}
        report=sizing(ir,{'vertices':[(0.,0.),(1.,0.)]}, {(0,1):{'R'}})
        self.assertFalse(report[0]['passed'])


class TdxStoredFieldTests(unittest.TestCase):
    def test_regional_storage_preserves_small_net_even_if_carriers_round(self):
        # A self-owned triangle: text rounding loses the carrier difference,
        # while TDX correctly retains the separately stored net field.
        fields={'DonorConcentration':{'values':[3e18+512]*3},
                'AcceptorConcentration':{'values':[3e18]*3},
                'DopingConcentration':{'values':[512.]*3}}
        bundle={'nodes':[{'vela_id':i,'x':x,'y':y} for i,(x,y) in enumerate([(0.,0.),(1.,0.),(0.,1.)])],
                'regions':[{'id':1,'name':'R'}],'materials':[{'id':1,'name':'Si'}],
                'triangles':[{'region_id':1,'material_id':1,'node_ids':[0,1,2]}],
                'contacts':[],'interfaces':[],'lines':[],'fields':fields,'provenance':{'material_aliases':{'Silicon':'Si'}}}
        grid='''Vertices (3) {
0 0
1 0
0 1
}
Edges (3) {
0 1
1 2
2 0
}
Elements (1) {
2 0 1 2
}
Region ("R") {
material = Silicon
Elements (1) {
0
}
}
'''
        def dataset(name,value):
            return f'Dataset ("{name}") {{\nvalidity = [ "R" ]\nValues (3) {{\n{value} {value} {value}\n}}\n}}\n'
        data=dataset('DonorConcentration','3e18')+dataset('AcceptorConcentration','3e18')+dataset('DopingConcentration','512')
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'readback';path.with_suffix('.grd').write_text(grid);path.with_suffix('.dat').write_text(data)
            result=validate_pair(bundle,path)
            self.assertEqual('passed',result['status'])
            self.assertEqual(0.,result['maximum_stored_field_relative_errors']['DopingConcentration'])
            path.with_suffix('.dat').write_text(data.replace('512 512 512','0 0 0'))
            with self.assertRaisesRegex(ValueError,'stored DopingConcentration'):validate_pair(bundle,path)


if __name__=='__main__':unittest.main()
