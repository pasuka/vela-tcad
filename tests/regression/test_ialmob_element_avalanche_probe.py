"""Check explicit source support, SI units and scope of the avalanche probe."""
import argparse
import copy
import json
import math
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--runner', type=Path, required=True)
    args = parser.parse_args()
    vertex = dict(potential_V=0., electron_qf_V=0., hole_qf_V=0.,
                  donors_m3=1e22, acceptors_m3=0., electrons_m3=1e20,
                  holes_m3=1e19, orientation_family=100,
                  electron_response_m3_per_V=1e20/.025852,
                  hole_response_m3_per_V=1e19/.025852)
    cell = dict(cell_id=17, coordinates_m=[[0, 0], [1e-8, 0], [0, 1e-8]],
                interface_distance_m=[1e-9, 1e-9, 2e-9],
                vertex_measure_m2=[1e-17, 2e-17, 2e-17],
                partial_boundary_layer=False, boundary_tangent=[1, 0],
                touches_effective_electrode=False,
                avalanche_mobility_vertex=1,
                avalanche_source_measure_m2=[1e-17, 1.5e-17, 2.5e-17],
                vertices=[dict(vertex) for _ in range(3)])
    cell['vertices'][1]['electron_qf_V'] = .1
    cell['vertices'][2]['hole_qf_V'] = .1
    data = dict(temperature_K=300., reference_density_m3=0.,
                electron_parameters_cm={'100': {'mumax': 1000.}},
                hole_parameters_cm={'100': {'mumax': 500.}}, elements_SI=[cell])
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp)/'input.json'

        def run(value, success=True):
            path.write_text(json.dumps(value), encoding='utf-8')
            proc = subprocess.run([str(args.runner), str(path)], capture_output=True,
                                  text=True, timeout=30)
            if not success:
                assert proc.returncode != 0 and not proc.stdout
                return None
            assert proc.returncode == 0, proc.stderr
            result = json.loads(proc.stdout)
            assert result['schema'] == 'vela.ialmob.element_avalanche_probe.v1'
            assert len(result['elements']) == 1
            return result['elements'][0]

        result = run(data)
        assert result['cell_id'] == 17
        for carrier, density, a, b in [('electron', 1e20, 7.03e7, 1.231e8),
                                       ('hole', 1e19, 1.582e8, 2.036e8)]:
            row = result[carrier]
            assert math.isclose(row['density_m3'], density, rel_tol=1e-14)
            assert math.isclose(row['current_field_V_per_m'], 1e7, rel_tol=1e-14)
            assert math.isclose(row['alpha_per_m'], a*math.exp(-b/1e7), rel_tol=1e-14)
            g = row['generation_m3_per_s']
            assert math.isclose(g, row['alpha_per_m']*density*
                                row['mobility_m2_per_Vs']*1e7, rel_tol=1e-14)
            derivatives = row['generation_derivative_per_V']
            assert len(derivatives) == 9 and all(math.isfinite(x) for x in derivatives)
            assert math.isclose(sum(x['value'] for x in row['source_per_m_s']),
                                g*5e-17, rel_tol=1e-14)
        moved = copy.deepcopy(data)
        moved['elements_SI'][0]['avalanche_source_measure_m2'].reverse()
        redistributed = run(moved)
        for carrier in ['electron', 'hole']:
            assert redistributed[carrier]['generation_m3_per_s'] == result[carrier]['generation_m3_per_s']
            assert redistributed[carrier]['generation_derivative_per_V'] == result[carrier]['generation_derivative_per_V']
            assert redistributed[carrier]['source_per_m_s'] == list(reversed(result[carrier]['source_per_m_s']))
        run(dict(data, temperature_K=350.), False)
        run(dict(data, reference_density_m3=1e18), False)
        run(dict(data, unrecognized_model=True), False)
        missing = copy.deepcopy(data)
        del missing['elements_SI'][0]['avalanche_source_measure_m2']
        run(missing, False)
    print('IALMob avalanche probe SI, explicit support and scope checks passed')


if __name__ == '__main__':
    main()
