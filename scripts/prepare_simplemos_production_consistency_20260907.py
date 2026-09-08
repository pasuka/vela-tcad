"""Read-only coefficient-path audit on sealed joint/native states; no new DC."""
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import difflib
import json
from pathlib import Path
import subprocess
import time
import analyze_simplemos_joint_eight_point_20260907 as prior

a=prior.a; d=prior.d; p=prior.p; b=prior.b
LOCAL=p.REPO/'build-release/simplemos_production_consistency_20260907'
OUT=p.REPO/'reference_tcad/simplemos_sentaurus2022/production_consistency_20260907'
SOURCES=['src/equation/CoupledDDAssembler.cpp','src/equation/DDAssembler.cpp','src/equation/PoissonAssembler.cpp',
    'src/post/ContactCurrent.cpp','src/tools/vela_example_runner.cpp','src/solver/NewtonSolver.cpp','src/simulation/DCSweep.cpp',
    'include/vela/equation/AssemblerUtils.h','include/vela/equation/CoupledDDAssembler.h','include/vela/solver/NewtonSolver.h',
    'src/mesh/BoxGeometryBuilder.cpp','docs/config_schema.md']


def prepare():
    a.verify(prior.OUT/'validation_evidence.json'); a.verify(b.OUT/'validation_evidence.json')
    LOCAL.mkdir(parents=True,exist_ok=False)
    files=[Path(__file__).resolve(),prior.OUT/'validation_evidence.json',b.OUT/'validation_evidence.json',b.RUNNER]
    files += [p.REPO/x for x in SOURCES]
    for src,copyfile in [('src/equation/CoupledDDAssembler.cpp',b.LOCAL/'CoupledDDAssembler.cpp'),('src/post/ContactCurrent.cpp',b.old.LOCAL/'ContactCurrent.cpp')]:
        dest=LOCAL/(Path(src).stem+'.runtime_source.diff')
        dest.write_text(''.join(difflib.unified_diff((p.REPO/src).read_text(encoding='utf-8').splitlines(True),copyfile.read_text().splitlines(True),fromfile=src,tofile=a.rel(copyfile))),encoding='utf-8')
        files += [dest,copyfile]
    # The active residual body is exactly current production source. Remaining
    # overlays are constructor coefficients/disabled charge hooks, diagnostics,
    # and Jacobian-only changes. Do not claim production Jacobian qualification.
    raw=(p.REPO/'src/equation/CoupledDDAssembler.cpp').read_text(encoding='utf-8'); overlay=(b.LOCAL/'CoupledDDAssembler.cpp').read_text()
    start='VectorXd CoupledDDAssembler::residualImpl('; end='CoupledDDResidualBlockNorms CoupledDDAssembler::'
    def residual_body(text):
        begin=text.index(start); boundary=text.index('\n}\n',begin)+3
        return text[begin:boundary]
    assert residual_body(raw)==residual_body(overlay)
    cases=[]
    for c in prior.cases():
        c=copy.deepcopy(c); c['probes']=[]; joint=Path(c['root'])/'joint/replacement/state.csv'
        for mode in ('legacy','builtin','joint','native_state_joint'):
            root=LOCAL/c['key']/mode; cfg=a.read(Path(c['base'])/'config.json'); cfg.pop('output_state_file'); cfg['solver'].pop('local_update_diagnostics')
            cfg['solver']['carrier_row_convergence']['mode']='report'; cfg['solver']['global_continuity_closure']['mode']='off'
            cfg['state_file']=c['native_initial'] if mode=='native_state_joint' else str(joint)
            assert 'region_resolved_interface_assembly' not in cfg['solver'] or not cfg['solver']['region_resolved_interface_assembly']
            if mode=='builtin':cfg['solver']['region_resolved_interface_assembly']=dict(transport_edge_coupling=True,
                poisson_electron_transport_node_volume=True,poisson_hole_transport_node_volume=True,poisson_dopant_transport_node_volume=True)
            strength=1. if mode in ('joint','native_state_joint') else 0.
            probes=[('functional','terminal_current_functional_probe'),('edges','sg_edge_flux_probe'),('terms','newton_carrier_term_probe')]
            if mode in ('builtin','joint'):probes.append(('mobility','edge_mobility_probe'))
            for name,kind in probes:
                deck=copy.deepcopy(cfg); deck['simulation_type']=kind; deck['contact']='drain'
                if name=='functional':deck.update(residual_output_csv=str(root/'residual.csv'),contact_edge_output_csv=str(root/'contact_edges.csv'))
                else:deck['output_csv']=str(root/(name+'.csv'))
                if name=='terms':deck['carrier_term_probe']={'solved_equation_terms':True}
                dest=root/(name+'.json'); a.write(dest,deck); files.append(dest)
                c['probes'].append(dict(mode=mode,strength=strength,config=str(dest)))
        cases.append(c); files += [joint,Path(c['native_initial']),Path(c['base'])/'config.json']
        files += [Path(c['ratios'])/n for n in ('edges.txt','nodes.txt')]
    a.write(OUT/'contract.json',dict(status='frozen_before_fixed_state_audit',cases=cases,total_probes=sum(len(c['probes']) for c in cases),
        source_files=SOURCES,production_changes=False,new_DC=0,acceptance_changes=False,
        harness='Sealed stable runner; production residualImpl text is identical. Native coefficient overlay is disabled for legacy/builtin. No charge/vector/mobility-step overlay. No nonlinear solve or claim that production Jacobian includes isolated fixes.',
        modes=dict(legacy='Current default coefficients at qualified joint state',builtin='Only transport_edge_coupling plus all three Poisson transport-volume flags at same state',joint='Sealed native Si coefficients and volumes at qualified joint state',native_state_joint='Same joint operator at native-coherent state'),
        gates=dict(port_relative=1e-8,source_replay_relative=1e-8,edge_sum_roundoff_eps=128,native_qf_SG_port_relative=1e-6),
        local_nodes=[1000,1009,320,324,338,1091,1092],
        interpretation='Native SG replay is calibrated by terminal sums. Local conductance/drop decomposition is an algebraic attribution, not an independently calibrated continuity perturbation or causal Id prediction. Native Poisson replay is not an exported native row residual.'))
    files.append(OUT/'contract.json'); d.matrix.freeze(OUT/'freeze.json',files)
    print('Frozen',sum(len(c['probes']) for c in cases),'read-only probes on eight states; production sources retained.',flush=True)


def run():
    a.verify(OUT/'freeze.json')
    def one(c):
        for job in c['probes']:
            cfg=Path(job['config']); status=cfg.with_suffix('.status.json')
            if status.exists():continue
            env=b.env(cfg.parent,Path(c['ratios']),job['strength'],job['strength'])
            for key in list(env):
                if key.startswith(('VELA_DIAGNOSTIC_','VELA_VALIDATE_','VELA_MINORITY_','VELA_SIMPLEMOS_LINEAR_')):env.pop(key)
            started=time.monotonic(); r=subprocess.run([str(b.RUNNER),'--config',str(cfg),'--log','off'],env=env,capture_output=True,text=True)
            cfg.with_suffix('.stdout.txt').write_text(r.stdout); cfg.with_suffix('.stderr.txt').write_text(r.stderr)
            s=json.loads(r.stdout.strip().splitlines()[-1]); s.update(exit_code=r.returncode,elapsed_seconds=time.monotonic()-started)
            a.write(status,s); assert r.returncode==0,(job,s)
        print('Exported physical paths',c['key'],flush=True)
    with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(one,a.read(OUT/'contract.json')['cases']))
    d.matrix.freeze(OUT/'probe_evidence.json',[OUT/'freeze.json']+[f for f in LOCAL.rglob('*') if f.is_file()])


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=('prepare','run')); globals()[parser.parse_args().action]()
