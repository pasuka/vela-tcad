"""Check SI units and fail-closed boundaries of fixed-state IALMob probes."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--geometry',type=Path,required=True)
    parser.add_argument('--element',type=Path,required=True)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        path=Path(tmp)/'input.json'
        def run(runner,data,success=True):
            path.write_text(json.dumps(data),encoding='utf-8')
            r=subprocess.run([str(runner),str(path)],capture_output=True,text=True,timeout=30)
            if success:
                assert r.returncode==0,r.stderr
                return json.loads(r.stdout)
            assert r.returncode!=0 and not r.stdout
        geom=dict(coordinates_m=[[0,0],[1e-6,0],[1e-6,1e-6],[.6e-6,.2e-6]],
            segments=[dict(node0=0,node1=1,semiconductor_point_m=[.5e-6,.5e-6]),
                      dict(node0=1,node1=2,semiconductor_point_m=[.5e-6,.5e-6])],
            crystal_x=[1,0,0],crystal_y=[0,1,0])
        r=run(args.geometry,geom)['nodes']
        assert r[3]['orientation_family']==110 and math.isclose(r[3]['distance_m'],.2e-6,rel_tol=1e-12)
        run(args.geometry,dict(geom,native_orientation=[100]*4),False)
        run(args.geometry,dict(geom,crystal_y=[1,0,0]),False)
        vertex=dict(potential_V=0.,electron_qf_V=0.,hole_qf_V=0.,donors_m3=0.,acceptors_m3=0.,
            electrons_m3=0.,holes_m3=0.,orientation_family=100,electron_response_m3_per_V=0.,hole_response_m3_per_V=0.)
        cell=dict(cell_id=17,coordinates_m=[[0,0],[1e-7,0],[0,1e-7]],interface_distance_m=[0,0,1e-7],
            vertex_measure_m2=[1e-15,2e-15,2e-15],partial_boundary_layer=False,boundary_tangent=[0,0],
            touches_effective_electrode=False,vertices=[dict(vertex) for _ in range(3)])
        data=dict(temperature_K=300.,electron_parameters_cm={'100':{'mumax':1000.}},
            hole_parameters_cm={'100':{'mumax':500.}},elements_SI=[cell])
        r=run(args.element,data)['elements'][0]
        assert r['cell_id']==17
        assert math.isclose(r['electron_m2_per_Vs'],.1,rel_tol=1e-14)
        assert math.isclose(r['hole_m2_per_Vs'],.05,rel_tol=1e-14)
        assert r['electron_derivative_per_V']==[0.]*9
        run(args.element,dict(data,temperature_K=350.),False)
        cell['vertices'][0]['orientation_family']=111
        run(args.element,data,False)
        cell['vertices'][0]['orientation_family']=100
        cell['vertices'][0]['electrons_m3']=-1.
        run(args.element,data,False)
    print('IALMob geometry/element probe SI and scope checks passed')


if __name__=='__main__':
    main()
