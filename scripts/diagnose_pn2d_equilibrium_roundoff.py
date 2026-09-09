#!/usr/bin/env python3
"""Frozen exported-state test of equilibrium quasi-Fermi roundoff sensitivity."""
import argparse
import csv
from pathlib import Path
from pn2d_variants import DEFAULT_ROOT, REPO, read_json, write_json, run_cmd, sha
from compare_pn2d_variants import node_field, rows
from audit_genius_bjt_conservative_sections import section_current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--id',default='D1')
    parser.add_argument('--mesh',default='M2')
    parser.add_argument('--fields-root',type=Path)
    args=parser.parse_args();root=args.root.resolve();source=root/args.id/args.mesh
    out=root/'diagnostics/equilibrium_roundoff'/f'{args.id}_{args.mesh}';out.mkdir(parents=True,exist_ok=True)
    mesh=read_json(source/'inputs/mesh.json');ids=[n['id'] for n in mesh['nodes']]
    archived=source/'initial_native_double/fields/forward_p0p00/fields'
    fields=args.fields_root or (archived if archived.exists() else source/'fields/forward_p0p00/fields')
    native={name:node_field(fields/f'{field}_region0.csv',ids) for name,field in
            [('psi','ElectrostaticPotential'),('phin','eQuasiFermiPotential'),('phip','hQuasiFermiPotential'),('electrons_m3','eDensity'),('holes_m3','hDensity')]}
    result={'status':'diagnostic_only','scope':'Vela production flux on frozen exported native equilibrium, never an independent solve or native internal-state replay',
            'state_policy':'Identical psi/n/p in both probes; only replace exported equilibrium quasi-Fermi roundoff by exact zero in the second probe.',
            'id':args.id,'mesh':args.mesh,'runner_sha256':sha(REPO/'build-release/vela_example_runner.exe'),
            'input_field_sha256':{p.name:sha(p) for p in fields.glob('*.csv')},
            'native_qfn_max_abs_V':max(map(abs,native['phin'])),'native_qfp_max_abs_V':max(map(abs,native['phip'])),'probes':[]}
    for label in ['exported_qf','exact_equilibrium_qf']:
        state=out/f'{label}.csv'
        with state.open('w',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['node_id',*native])
            for i,identifier in enumerate(ids):
                writer.writerow([identifier,*[(0 if label=='exact_equilibrium_qf' else value[i]) if key in ['phin','phip'] else value[i]*(1e6 if key.endswith('_m3') else 1) for key,value in native.items()]])
        cfg=read_json(source/'probe.json')
        cfg.update(simulation_type='sg_edge_flux_probe',state_file=str(state),output_csv=str(out/f'{label}_edges.csv'))
        for key in ['mesh_file','node_doping_file','materials_file']: cfg[key]=str(source/cfg[key])
        for contact in cfg['contacts']: contact['bias']=0
        config=out/f'{label}.json';write_json(config,cfg)
        rc=run_cmd([REPO/'build-release/vela_example_runner.exe','--config',config],REPO,out/f'{label}.log')
        if rc: raise RuntimeError(f'frozen diagnostic failed: {label}')
        edges=rows(Path(cfg['output_csv']))
        result['probes'].append({'name':label,'returncode':rc,'state_sha256':sha(state),'config_sha256':sha(config),
            'sections':[{'x_um':x,**section_current(edges,axis='x',cut_um=x)} for x in read_json(root/'contract.json')['conservation']['sections_x_um']]})
    write_json(out/'result.json',result)


if __name__=='__main__': main()
