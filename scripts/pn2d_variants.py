#!/usr/bin/env python3
"""Generate and execute the prospective PN2D incremental validation campaign.

All binary artifacts stay in an ignored build directory. No reference state is
ever used to initialize a Vela solve. Each SDevice target is a separate Goal,
so exported comparison voltages are solved endpoints, not interpolated plots.
"""
from __future__ import annotations

import argparse
import csv
import copy
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

from generate_pn2d_config import render_named_template
from convert_tcad_export import load_mesh
from sentaurus_import import parse_quoted_list, parse_values_block
from sentaurus_parameter_ir import parse_parameter_ir, _strip_comments

REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / 'reference_tcad/pn2d_sentaurus2018'
SPEC = FIXTURE / 'variants/experiments.json'
CONTRACT = FIXTURE / 'variants/acceptance.json'
DEFAULT_ROOT = REPO / 'build-release/pn2d_variants'


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def biases(branch):
    return [round(i * branch['step_V'], 10) for i in range(round(branch['stop_V'] / branch['step_V']) + 1)]


def validate_spec(spec):
    if spec['geometry_um']!={'length':2.0,'height':0.5,'junction_x':1.0} or spec['temperature_K']!=300:
        raise ValueError('This campaign supports only the frozen straight 300 K PN geometry')
    if spec['area_factor']!=1 or spec['depth_um']!=1:
        raise ValueError('The native A/um mapping requires AreaFactor=1 and depth=1 um')
    if spec['constant_mobility_cm2_Vs']!={'electron':1417.0,'hole':470.5} or spec['srh_lifetime_s']!={'electron':1e-5,'hole':3e-6}:
        raise ValueError('Manifest physics constants must match the explicit shipped models.par laws')
    expected=['P0','P1','P2','P3','D1','D2','G1','G2']
    if [e['id'] for e in spec['experiments']]!=expected:
        raise ValueError('Expected exactly eight ordered unique configurations')
    for e in spec['experiments']:
        if not set(e['overrides'])<=set(spec['baseline']) or len(e['overrides'])!=(0 if e['id']=='P0' else 1):
            raise ValueError('Each variant must change exactly one declared factor')


def tag(v):
    return ('m' if v < 0 else 'p') + f'{abs(v):.2f}'.replace('.', 'p')


def material(bgn):
    # Freeze the shipped Silicon DOS and gap laws at 300 K, not fitted ni.
    T = 300.0
    eg0 = 1.16964
    eg = eg0 - 4.73e-4 * T*T / (636.0 + T)
    me = ((6 * 0.1905 * eg0 / eg)**2 * 0.9163)**(1/3)
    numerator = 0.443587 + 0.003609528*T + 0.0001173515*T*T + 1.263218e-6*T**3 + 3.025581e-9*T**4
    denominator = 1 + 0.004683382*T + 0.0002286895*T*T + 7.469271e-7*T**3 + 1.727481e-9*T**4
    nc, nv = 2.5094e19 * me**1.5, 2.5094e19 * numerator / denominator
    effective_gap = eg + (-0.01595 if bgn else 0.0)
    ni = math.sqrt(nc*nv) * math.exp(-effective_gap/(2*8.617333262145e-5*T))
    return {'materials': [{'name': 'Si', 'eps_r': 11.7, 'ni': ni, 'mun': 1417.0, 'mup': 470.5,
                           'bandgap_eV': effective_gap, 'electron_affinity_eV': 4.05,
                           'Nc_m3': nc, 'Nv_m3': nv, 'temperature_K': T}]}


def sde_text(spec, cfg, mesh):
    text = (FIXTURE/'source/pn2d_sde.cmd').read_text(encoding='utf-8')
    text = text.replace('"BoronActiveConcentration"\n  1e17', f'"BoronActiveConcentration"\n  {cfg["NA_cm3"]:.12g}')
    text = text.replace('"PhosphorusActiveConcentration"\n  1e17', f'"PhosphorusActiveConcentration"\n  {cfg["ND_cm3"]:.12g}')
    scale = mesh['scale']
    text = re.sub(r'^  (0\.05 0\.05|0\.01 0\.01|0\.01 0\.02|0\.005 0\.005)$',
                  lambda m: '  '+' '.join(str(float(v)*scale) for v in m[1].split()), text, flags=re.M)
    if 'junction_x_scale' in mesh:
        # Resolve the narrow low-bias SRH peak without changing endpoint or
        # transverse refinement. Preserve M0/M1 evidence, never overwrite it.
        start = text.index('(sdedr:define-refinement-size\n  "Junction.Mesh"')
        end = text.index('\n)', start)
        sx = mesh['junction_x_scale']
        text = text[:start] + f'(sdedr:define-refinement-size\n  "Junction.Mesh"\n  {0.01*sx} {0.02*scale}\n  {0.005*sx} {0.005*scale}' + text[end:]
    if mesh['segmented']:
        endpoints = [0.125, 0.1875, 0.3125, 0.375]
        insertion = '\n'.join(f'(sdegeo:insert-vertex (position 0.0 {y} 0.0))' for y in endpoints) + '\n'
        text = text.replace('; Define contacts', insertion + '; Define contacts')
        start = text.index('(sdegeo:define-2d-contact')
        end = text.index('; Right electrode', start)
        low, high = 0.25*(1-cfg['anode_fraction']), 0.25*(1+cfg['anode_fraction'])
        ys = [0.0, *endpoints, 0.5]
        contacts = '\n'.join(f'(sdegeo:define-2d-contact (find-edge-id (position 0.0 {(a+b)/2} 0.0)) "Anode")'
                             for a,b in zip(ys,ys[1:]) if a >= low-1e-12 and b <= high+1e-12)
        text = text[:start] + contacts + '\n\n' + text[end:]
        edge = ''
        edge_scale=mesh.get('contact_edge_scale',scale)
        for i,y in enumerate(endpoints):
            edge += f'\n(sdedr:define-refeval-window "Edge{i}" "Rectangle" (position 0.0 {y-0.025} 0.0) (position 0.05 {y+0.025} 0.0))\n'
            edge += f'(sdedr:define-refinement-size "EdgeMesh{i}" {0.01*edge_scale} {0.00625*edge_scale} {0.005*edge_scale} {0.003125*edge_scale})\n'
            edge += f'(sdedr:define-refinement-placement "EdgePlace{i}" "EdgeMesh{i}" "Edge{i}")\n'
        text = text.replace('(sde:build-mesh', edge + '\n(sde:build-mesh')
    return text


def sdevice_text(spec, cfg, branch_name, mesh=None):
    branch = spec['branches'][branch_name]
    mesh=mesh or {}
    initial_step=mesh.get('sdevice_goal_initial_step',0.05)
    bits=mesh.get('sentaurus_precision_bits')
    if bits not in (None,80,128): raise ValueError('Unsupported native arithmetic precision')
    precision=f'ExtendedPrecision({bits}) ' if bits else ''
    rhs=mesh.get('sdevice_rhs_min')
    if rhs is not None:
        if not 0 < rhs < 1: raise ValueError('Native RHSMin must be positive and below one')
        precision+=f'RHSMin={rhs:g} '
    mob = 'Mobility(DopingDependence)' if cfg['mobility']=='masetti' else '# No mobility dependence: ConstantMobility parameter law at 300 K'
    srh = 'Recombination(SRH)' if cfg['srh'] else '# Recombination disabled'
    bgn = 'OldSlotboom' if cfg['bgn'] else 'NoBandGapNarrowing'
    text = f'''File {{
 Grid="pn2d_msh.tdr" Parameter="models.par"
 Plot="{branch_name}_des.tdr" Current="{branch_name}.plt" Output="{branch_name}.log"
}}
Electrode {{ {{ Name="Anode" Voltage=0 }} {{ Name="Cathode" Voltage=0 }} }}
Physics {{ Temperature=300 AreaFactor=1
 {mob}
 {srh}
 EffectiveIntrinsicDensity({bgn})
}}
Plot {{ Potential ElectricField/Vector eDensity hDensity eQuasiFermi hQuasiFermi
 eCurrent/Vector hCurrent/Vector TotalCurrent/Vector Doping DonorConcentration AcceptorConcentration
 SRHRecombination eMobility hMobility EffectiveIntrinsicDensity BandGap BandGapNarrowing
}}
Math {{ {precision}Extrapolate RelErrControl Digits=8 Iterations=80 NotDamped=100 }}
Solve {{
 Coupled(Iterations=100) {{ Poisson }}
 Coupled(Iterations=100) {{ Poisson Electron Hole }}
 Plot(FilePrefix="{branch_name}_{tag(0)}")
'''
    for v in biases(branch)[1:]:
        text += f''' Quasistationary(InitialStep={initial_step:g} MinStep=1e-8 MaxStep=1 Increment=1.3 Goal {{ Name="Anode" Voltage={v:.12g} }}) {{
  Coupled {{ Poisson Electron Hole }}
 }}
'''
        if not mesh.get('plot_only_saved_biases') or v in spec['saved_biases_V']:
            text+=f' Plot(FilePrefix="{branch_name}_{tag(v)}")\n'
    return text+'}\n'


def vela_config(spec, cfg, branch):
    c, _ = render_named_template('pn2d_iv')
    c.update(state_format='hdf5', mesh_file='inputs/mesh.json', node_doping_file='inputs/doping.csv', materials_file='materials.json', output_csv=f'{branch}.csv')
    s = c['solver']
    s.update(temperature_K=300.0, taun=1e-5, taup=3e-6, recombination=['srh'] if cfg['srh'] else [],
             bandgap_narrowing='old_slotboom' if cfg['bgn'] else 'none', srh_doping_dependence={'enabled':False})
    s['mobility']['model'] = cfg['mobility']
    # Baseline A/B probe: source-only row equilibration stalls near 1 mV;
    # include flux in left scaling while retaining the exact equations and all
    # local/global convergence gates. This is campaign-local, not a default.
    s['continuity_row_scaling']['flux_fraction'] = 1.0
    c['sweep'].update(bias_points=biases(spec['branches'][branch]), step=spec['branches'][branch]['step_V'], stop=spec['branches'][branch]['stop_V'],
                       vtk_prefix=f'{branch}/state', write_vtk=True, write_state_file=f'{branch}/last.h5',
                       write_state_every_point_prefix=f'{branch}/accepted', stop_on_failure=True)
    c['sweep']['diagnostics'].update(terminal_balance={'enabled':True,'contacts':['Anode','Cathode']},
                                  continuity_balance={'enabled':True,'contacts':['Anode','Cathode']}, transport={'enabled':True})
    c['sweep']['diagnostics']['newton_history']['csv_file'] = f'{branch}_newton.csv'
    return c


def generate(root):
    spec = read_json(SPEC)
    validate_spec(spec)
    root.mkdir(parents=True, exist_ok=True)
    if (root/'contract.json').exists() and sha(root/'contract.json') != sha(CONTRACT):
        raise ValueError('Contract changed: use a new campaign directory; existing evidence cannot be rescored silently')
    shutil.copyfile(CONTRACT, root/'contract.json')
    entries=[]
    for exp in spec['experiments']:
        cfg={**spec['baseline'], **exp['overrides']}
        numerics=spec.get('sentaurus_numerics_by_stage',{}).get(exp['stage'],{})
        for mesh_name, mesh in spec['meshes'].items():
            if mesh_name=='original' and exp['id']!='P0':
                continue
            d=root/exp['id']/mesh_name
            d.mkdir(parents=True, exist_ok=True)
            (d/'pn2d_sde.cmd').write_text(sde_text(spec,cfg,mesh),encoding='utf-8')
            shutil.copyfile(FIXTURE/'source/models.par',d/'models.par')
            write_json(d/'materials.json',material(cfg['bgn']))
            for branch in spec['branches']:
                (d/f'{branch}_sdevice.cmd').write_text(sdevice_text(spec,cfg,branch,{**mesh,**numerics}),encoding='utf-8')
                write_json(d/f'vela_{branch}.json',vela_config(spec,cfg,branch))
            entry={**exp,'parameters':cfg,'mesh':mesh_name,'directory':d.relative_to(root).as_posix(),
                   'input_sha256':{name:sha(d/name) for name in ['pn2d_sde.cmd','models.par','materials.json','forward_sdevice.cmd','reverse_sdevice.cmd','vela_forward.json','vela_reverse.json']},
                   'status':'generated_not_run'}
            if numerics: entry['sentaurus_numerics']=numerics
            write_json(d/'manifest.json',entry)
            entries.append(entry)
    write_json(root/'campaign.json',{'schema':'vela.pn2d.campaign.v1','spec_sha256':sha(SPEC),'contract_sha256':sha(CONTRACT),'entries':entries})
    return entries


def export_review_inputs(root):
    """Save compact generated M0 text inputs for review, not binary outputs."""
    output=FIXTURE/'variants/inputs'
    for e in read_json(root/'campaign.json')['entries']:
        if e['mesh']!='M0': continue
        dest=output/e['id']; dest.mkdir(parents=True,exist_ok=True)
        source=root/e['directory']
        names=['pn2d_sde.cmd','forward_sdevice.cmd','reverse_sdevice.cmd','vela_forward.json','vela_reverse.json','materials.json']
        for name in names: shutil.copyfile(source/name,dest/name)
        write_json(dest/'manifest.json',{**e,'parameter_file':'../../../source/models.par','parameter_sha256':sha(FIXTURE/'source/models.par')})


def audit_parameters(root, dump):
    if dump.resolve()==(FIXTURE/'source/models.par').resolve():
        raise ValueError('A current-version remote parameter dump is required; self-comparison is not reference evidence')
    sections={'Epsilon','Bandgap','OldSlotboom','eDOSMass','hDOSMass','ConstantMobility','DopingDependence','Scharfetter'}
    def values(path):
        # The 2022 dump contains inactive Species/Ionization syntax beyond the
        # shared IR parser's scope. Project only the audited complete blocks;
        # never drop statements from inside an active block to make it parse.
        lines=path.read_text(encoding='utf-8').splitlines(keepends=True)
        projected=[]; active=False; depth=0; opened=False; found=[]
        for line in lines:
            code=_strip_comments(line)[0]
            if not active:
                name=next((s for s in sections if re.match(r'^'+s+r'\b',code)),None)
                if name is None: continue
                active=True; depth=0; opened=False; found.append(name)
            projected.append(line)
            depth+=code.count('{')-code.count('}')
            opened=opened or '{' in code
            if opened and depth==0: active=False
        if active or set(found)!=sections or len(found)!=len(sections):
            raise ValueError('Active parameter blocks are missing, duplicated or unbalanced')
        directory=root/'parameter_audit';directory.mkdir(exist_ok=True)
        selected_file=directory/(path.stem+'_selected.par')
        selected_file.write_text('Material="Silicon" {\n'+''.join(projected)+'\n}\n',encoding='utf-8')
        result={}
        for block in parse_parameter_ir(selected_file)['blocks']:
            if block['section'] not in sections: continue
            for p in block['parameters']:
                result[block['section']+'.'+p['raw_name']]=p.get('values',p['raw_lexeme'])
        return result
    historical=values(FIXTURE/'source/models.par'); current=values(dump)
    differences=[{'parameter':k,'source':historical.get(k),'current_default':current.get(k)} for k in sorted(historical.keys()|current.keys()) if historical.get(k)!=current.get(k)]
    for difference in differences:
        key=difference['parameter']
        difference['active']=not (key.startswith('DopingDependence.Ar_') or key in ['Bandgap.CNTdiameter','Bandgap.FermiVelocity','eDOSMass.Nc300','hDOSMass.Nv300'])
    result={'schema':'vela.pn2d.parameter_audit.v1','source_sha256':sha(FIXTURE/'source/models.par'),'current_dump_sha256':sha(dump),
            'sections':sorted(sections),'source_parameters':historical,'current_default_parameters':current,
            'differences':differences,'status':'active_parameters_match' if not any(x['active'] for x in differences) else 'active_parameter_difference',
            'note':'A parameter library is not activation: SRH doping dependence, high-field saturation, Fermi and avalanche are not selected. The campaign explicitly loads the source parameter file.'}
    write_json(root/'parameter_audit.json',result)
    return result


def run_cmd(argv, cwd, log, timeout=3600):
    with Path(log).open('w',encoding='utf-8') as out:
        try:
            p=subprocess.run([str(x) for x in argv],cwd=cwd,stdout=out,stderr=subprocess.STDOUT,timeout=timeout)
            return p.returncode
        except subprocess.TimeoutExpired:
            out.write('\nTIMEOUT\n')
            return 124


def selected(root, ids, meshes, allow_modeling=False):
    entries=[e for e in read_json(root/'campaign.json')['entries'] if e['id'] in ids.split(',') and e['mesh'] in meshes.split(',')]
    if not entries: raise ValueError('selection contains no experiments')
    if not allow_modeling and any(e['stage'] in ['B','C'] for e in entries):
        report=read_json(root/'results.json')
        if not report['stage_A_accepted'] or report['contract_sha256']!=sha(root/'contract.json'):
            raise ValueError('Stage A has not passed the frozen contract; B/C simulation is gated')
        pair=read_json(SPEC.parent/'mesh_refinement.json')['additional_pair']
        evidence={(r['id'],r['mesh']):r for r in report.get('results',[])}
        prerequisites=['P0','P1','P2','P3']
        if any(e['stage']=='C' for e in entries):
            if any(report.get('refined_mesh',{}).get(i,{}).get('status')!='pass' or
                   any(evidence.get((i,m),{}).get('status')!='pass' for m in pair) for i in ['D1','D2']):
                raise ValueError('Stage B must finish qualification before stage C simulation')
            prerequisites+=['D1','D2']
        for identifier in prerequisites:
            prerequisite_stage='A' if identifier.startswith('P') else 'B'
            for mesh in pair:
                case=evidence.get((identifier,mesh),{})
                if case.get('status')!='pass':
                    raise ValueError(f'Stage {prerequisite_stage} gate requires every current refined-pair result')
                d=root/identifier/mesh
                required=['manifest.json','pn2d_sde.cmd','pn2d_msh.tdr','inputs/mesh.json','inputs/doping.csv',
                          'vela_forward.json','vela_reverse.json','materials.json','models.par',
                          'forward_sdevice.cmd','reverse_sdevice.cmd','vela_forward_run.json','vela_reverse_run.json','sentaurus_run.json']
                required.extend(name for name in ['vela_forward_runtime.json','vela_reverse_runtime.json']
                                if case.get('evidence_sha256',{}).get(name))
                for name in required:
                    digest=case.get('evidence_sha256',{}).get(name)
                    if not digest or not (d/name).exists() or sha(d/name)!=digest:
                        raise ValueError(f'Stage {prerequisite_stage} evidence changed or is missing: {identifier}/{mesh}/{name}; compare again')
    return entries


def remote_run(root, ids, meshes, target, license_server=None, workers=1, mesh_only=False):
    entries=selected(root,ids,meshes,allow_modeling=mesh_only)
    # Largest meshes first reduces the queue tail for independent native jobs.
    mesh_rank={name:i for i,name in enumerate(read_json(SPEC)['meshes'])}
    entries.sort(key=lambda entry:mesh_rank[entry['mesh']],reverse=True)
    success=True
    if license_server not in (None,'27000@127.0.0.1','27000@tcad'):
        raise ValueError('Only the documented VM license server is supported')
    license_prefix=f'SNPSLMD_LICENSE_FILE={license_server} ' if license_server else ''
    # Upload only campaign-generated text into a unique, task-owned directory.
    remote_record=root/'remote_root.json'
    if remote_record.exists():
        remote=read_json(remote_record)['directory']
    else:
        completed=root/'P0/original/sentaurus_run.json'
        remote=read_json(completed)['remote'].rsplit('/P0/',1)[0] if completed.exists() else f'/home/tcad/sentaurus_runs/vela_oracle_2022/pn2d_variants_{sha(root/"campaign.json")[:12]}'
        if not re.fullmatch(r'/home/tcad/sentaurus_runs/vela_oracle_2022/pn2d_variants_[a-f0-9]{12}',remote):
            raise ValueError('invalid task-owned remote directory')
        write_json(remote_record,{'directory':remote})
    ssh=['C:/Windows/System32/OpenSSH/ssh.exe','-o','BatchMode=yes','-o','ConnectTimeout=15',target]
    scp=['C:/Windows/System32/OpenSSH/scp.exe','-o','BatchMode=yes']
    banner=subprocess.run(ssh+['sdevice -h'],capture_output=True,text=True,timeout=45)
    version_text=banner.stdout+banner.stderr
    (root/'sentaurus_version.txt').write_text(version_text,encoding='utf-8')
    if read_json(SPEC)['sentaurus_version'] not in version_text:
        raise RuntimeError('Sentaurus version mismatch or unavailable')
    archive=root/'upload.tar'
    with tarfile.open(archive,'w') as tar:
        for e in entries:
            d=root/e['directory']
            for name in ['pn2d_sde.cmd','models.par','forward_sdevice.cmd','reverse_sdevice.cmd']:
                tar.add(d/name,arcname=f'{e["directory"]}/{name}')
    subprocess.run(ssh+[f'mkdir -p {remote}'],check=True)
    subprocess.run(scp+[str(archive),f'{target}:{remote}/upload.tar'],check=True)
    subprocess.run(ssh+[f'cd {remote} && tar xf upload.tar'],check=True)
    def execute(e):
        rel=e['directory']; d=root/rel
        print(f'Sentaurus {rel}',flush=True)
        mesh_command=(f'cp ../../P0/{e["mesh"]}/pn2d_msh.tdr pn2d_msh.tdr' if e['id'] in ['P1','P2','P3'] else f'{license_prefix}sde -e -l pn2d_sde.cmd > sde.out 2>&1')
        mesh_record=d/'mesh_run.json'
        if not mesh_only and mesh_record.exists():
            previous=read_json(mesh_record)
            if previous['returncode']==0 and previous['input_sha256']['pn2d_sde.cmd']==e['input_sha256']['pn2d_sde.cmd']:
                mesh_command='test -f pn2d_msh.tdr'
        record=d/('mesh_run.json' if mesh_only else 'sentaurus_run.json')
        write_json(record,{'returncode':None,'fetch_returncode':None,'status':'running','remote':remote+'/'+rel,'version':read_json(SPEC)['sentaurus_version']})
        command=f'cd {remote}/{rel} && {mesh_command}'
        if not mesh_only:
            command+=f' && {license_prefix}sdevice forward_sdevice.cmd > forward.out 2>&1 && {license_prefix}sdevice reverse_sdevice.cmd > reverse.out 2>&1'
        rc=run_cmd(ssh+[command],root,d/('remote_mesh.log' if mesh_only else 'remote.log'),timeout=10800)
        # Fetch even failed partial runs; never turn an exit code into acceptance.
        # Transfer required full states, curves, inputs and logs in one archive.
        # Extra native Plot endpoints remain in the isolated VM run directory.
        states=' '.join(f'{"reverse" if b<0 else "forward"}_{tag(b)}_des.tdr' for b in read_json(SPEC)['saved_biases_V'])
        pack=f'cd {remote}/{rel} && tar --ignore-failed-read -czf transfer.tgz -- *.cmd models.par pn2d_msh.tdr *.out *.log *.plt {states if not mesh_only else ""}'
        packed=subprocess.run(ssh+[pack],capture_output=True,text=True)
        fetch=subprocess.run(scp+[f'{target}:{remote}/{rel}/transfer.tgz',str(d/'transfer.tgz')],capture_output=True,text=True)
        if fetch.returncode==0:
            with tarfile.open(d/'transfer.tgz') as archive: archive.extractall(d,filter='data')
        write_json(record,{'returncode':rc,'fetch_returncode':fetch.returncode or packed.returncode,'remote':remote+'/'+rel,'version':read_json(SPEC)['sentaurus_version'],'license_server':license_server or 'inherited','input_sha256':e['input_sha256'],'scope':'mesh_only' if mesh_only else 'mesh_and_dd'})
        return rc==0 and fetch.returncode==0 and packed.returncode==0
    # Cases own disjoint directories and use the same frozen solver strategy.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        completed=list(pool.map(execute,entries))
    return all(completed)


def postprocess_states(d, branch, source_config):
    """Recover full fields from Vela's own solved states; never run transport here."""
    manifest_path=d/'spatial_artifacts.json'
    (d/'postprocess').mkdir(exist_ok=True)
    manifest=read_json(manifest_path) if manifest_path.exists() else {}
    for bias in read_json(SPEC)['saved_biases_V']:
        if (bias<0)!=(branch=='reverse'): continue
        token=('m' if bias<0 else '')+f'{abs(bias):.6f}'.replace('.','p')
        state=f'{branch}/accepted_bias_{token}.h5'
        cfg=copy.deepcopy(source_config)
        cfg['output_csv']=f'postprocess/{branch}_{tag(bias)}.csv'
        cfg['solver']['method']='frozen_state'
        sweep=cfg['sweep']
        sweep.pop('initialization',None)
        sweep.pop('write_state_file',None)
        sweep.pop('write_state_every_point_prefix',None)
        sweep.update(start=bias,stop=bias,bias_points=[bias],initial_state_file=state,write_vtk=True,
                     vtk_prefix=f'postprocess/{branch}_{tag(bias)}',diagnostics={})
        path=d/f'postprocess_{branch}_{tag(bias)}.json';write_json(path,cfg)
        rc=run_cmd([REPO/'build-release/vela_example_runner.exe','--config',path],REPO,path.with_suffix('.stdout'))
        if rc: raise RuntimeError(f'Vela state postprocessing failed: {path}')
        candidates=list((d/'postprocess').glob(f'{branch}_{tag(bias)}_*V.vtk'))
        if len(candidates)!=1: raise ValueError('postprocessing VTK is missing or ambiguous')
        manifest[f'{bias:g}']={'vtk':candidates[0].relative_to(d).as_posix(),'state_file':state,'state_sha256':sha(d/state),
                              'classification':'postprocessing_only_of_independently_solved_Vela_state','config_sha256':sha(path)}
    write_json(manifest_path,manifest)


def import_run(root, ids, meshes, compact_spatial=False, mesh_only=False):
    success=True
    for e in selected(root,ids,meshes,allow_modeling=mesh_only):
        d=root/e['directory']; out=d/'inputs'; out.mkdir(exist_ok=True)
        rc=run_cmd([REPO/'build-release/sentaurus_import.exe','--tdr',d/'pn2d_msh.tdr','--export-dir',out,'--inventory-json',out/'inventory.json'],REPO,d/'import.log')
        if rc: raise RuntimeError(f'Mesh import failed: {d}')
        mesh,_=load_mesh(out); write_json(out/'mesh.json',mesh)
        # Do not synthesize doping: the exact exported rows are the solver input.
        write_json(d/'input_audit.json',{'mesh_sha256':sha(out/'mesh.json'),'doping_sha256':sha(out/'doping.csv'),
                   'tdr_sha256':sha(d/'pn2d_msh.tdr'),'node_count':len(mesh['nodes']),'triangle_count':len(mesh['triangles']),
                   'contacts':[{k:c[k] for k in ('name','node_ids')} for c in mesh['contacts']], 'policy':'reported', 'status':'exported_requires_comparison'})
        if mesh_only:
            from compare_pn2d_variants import audit_imported_mesh
            write_json(d/'mesh_input_comparison.json',audit_imported_mesh(d,e,root))
            continue
        for branch in read_json(SPEC)['branches']:
            (d/branch).mkdir(exist_ok=True)
            print(f'Vela independent solve {e["directory"]} {branch}',flush=True)
            config_digest=sha(d/f'vela_{branch}.json')
            source_config=read_json(d/f'vela_{branch}.json')
            runtime_path=d/f'vela_{branch}.json'
            if compact_spatial:
                runtime=copy.deepcopy(source_config);runtime['sweep']['write_vtk']=False
                runtime_path=d/f'vela_{branch}_runtime.json';write_json(runtime_path,runtime)
            write_json(d/f'vela_{branch}_run.json',{'returncode':None,'status':'running','config_sha256':config_digest})
            rc=run_cmd([REPO/'build-release/vela_example_runner.exe','--config',runtime_path],REPO,d/f'vela_{branch}.stdout')
            write_json(d/f'vela_{branch}_run.json',{'returncode':rc,'runner_sha256':sha(REPO/'build-release/vela_example_runner.exe'),
                       'config_sha256':config_digest,
                       'runtime_config':runtime_path.name,'runtime_config_sha256':sha(runtime_path),
                       'compact_spatial':compact_spatial,
                       'linear_backend':os.environ.get('VELA_LINEAR_SOLVER','sparselu'),
                       'compiler':'MSYS2 UCRT64 GNU 16.2.0; Release preset',
                       'initialization':'independent poisson_block; Vela continuation only'})
            if rc==0 and compact_spatial: postprocess_states(d,branch,source_config)
            success=success and rc==0
    return success


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['generate','sentaurus','sentaurus-mesh','import-mesh','vela','audit-parameters'])
    p.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--ids',default='P0')
    p.add_argument('--meshes',default='original,M0,M1')
    p.add_argument('--ssh-target',default='sentaurus')
    p.add_argument('--license-server',choices=['27000@127.0.0.1','27000@tcad'])
    p.add_argument('--export-review-inputs',action='store_true')
    p.add_argument('--workers',type=int,choices=[1,2,3,4],default=1,
                   help='Independent native jobs; use at most the available VM CPUs and memory')
    p.add_argument('--parameter-dump',type=Path)
    p.add_argument('--compact-spatial',action='store_true')
    a=p.parse_args(); root=a.root.resolve()
    if a.command=='generate':
        generate(root)
        if a.export_review_inputs: export_review_inputs(root)
    elif a.command=='audit-parameters':
        if a.parameter_dump is None: p.error('--parameter-dump is required')
        audit_parameters(root,a.parameter_dump.resolve())
    elif a.command in ['sentaurus','sentaurus-mesh']:
        return 0 if remote_run(root,a.ids,a.meshes,a.ssh_target,a.license_server,a.workers,a.command=='sentaurus-mesh') else 1
    else: return 0 if import_run(root,a.ids,a.meshes,a.compact_spatial,a.command=='import-mesh') else 1
    return 0


if __name__=='__main__':
    sys.exit(main())
