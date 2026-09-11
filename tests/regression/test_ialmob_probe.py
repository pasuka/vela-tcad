"""Check the explicit units and fail-closed scope of the local mobility tool."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--runner',type=Path,required=True)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as temporary:
        path=Path(temporary)/'input.json'
        state=dict(id='intrinsic',donors_m3=0.,acceptors_m3=0.,electrons_m3=0.,holes_m3=0.,
                   normalField_V_per_m=0.,interfaceDistance_m=0.,crystal_normal=[1.,1.,0.])
        config=dict(temperature_K=300.,carrier='electron',parameters_cm={'mumax':1000.},states_SI=[state])
        def run():
            path.write_text(json.dumps(config),encoding='utf-8')
            return subprocess.run([str(args.runner),str(path)],text=True,encoding='utf-8',capture_output=True,timeout=30)
        result=run();assert result.returncode==0,result.stderr
        point=json.loads(result.stdout)['results'][0]
        assert math.isclose(point['mobility_m2_per_Vs'],.1,rel_tol=1e-14)
        assert point['id']=='intrinsic' and point['orientation_family']==110
        assert point['coulomb_m2_per_Vs'] is None
        config['derivatives']=True
        result=run();assert result.returncode==0,result.stderr
        point=json.loads(result.stdout)['results'][0]
        assert point['mobility_derivatives_SI']==[0.]*6
        # Unsupported temperatures/options must not silently select a surrogate.
        config['temperature_K']=350.
        result=run();assert result.returncode!=0 and not result.stdout and '300 K' in result.stderr
        config['temperature_K']=300.;config['parameters_cm']['tcoulomb']=.01
        result=run();assert result.returncode!=0 and not result.stdout and 'tcoulomb' in result.stderr
        del config['parameters_cm']['tcoulomb']
        state['electrons_m3']=-1.
        result=run();assert result.returncode!=0 and not result.stdout
        state['electrons_m3']=0.;state['electron_cm3']=1e16
        result=run();assert result.returncode!=0 and not result.stdout and 'electron_cm3' in result.stderr
    print('IALMob local probe units and rejection checks passed')


if __name__=='__main__':
    main()
