"""Four-point, read-only formula/flux/reconstruction audit; no nonlinear solves."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np
import run_simplemos_m80b_state_semantics as upstream
from compare_genius_bjt_transport_fields import read_vtk_point_data
from audit_genius_bjt_srh_fixed_state import RATE_COLUMNS

old = upstream.old
m73 = upstream.original.m73
REPO, ROOT = upstream.REPO, upstream.ROOT
LOCAL = REPO / 'build-release/simplemos_fixed_state_formula_audit'
OUT = ROOT / 'fixed_state_formula_audit'
CONTRACT = OUT / 'contract.json'
EVIDENCE = OUT / 'evidence.json'
SCRIPT = Path(__file__).resolve()
Q = 1.602176634e-19
FIELDS = ('eDensity', 'hDensity', 'EffectiveIntrinsicDensity', 'ElectrostaticPotential',
          'eQuasiFermiPotential', 'hQuasiFermiPotential', 'srhRecombination')


def workflows():
    return [w for w in old.read_json(old.m74.PORTABLE_MANIFEST)['workflows'] if w['device'] in ('n19', 'n23')]


def paths(w):
    tag = f"{w['device']}_vd_{old.m74.voltage_tag(w['drain_voltage_V'])}"
    source = REPO / w['phases'][-1]['source_config']
    export = m73.exports()[(w['device'], float(w['drain_voltage_V']), 'gate')]
    mapped = upstream.LOCAL / w['device'] / f"vd_{old.m74.voltage_tag(w['drain_voltage_V'])}" / 'sentaurus_mapped_state.csv'
    return tag, source, export, mapped


def verify_inputs():
    c = old.read_json(CONTRACT)
    for p, digest in c['input_hashes'].items():
        if old.sha256(REPO / p) != digest:
            raise ValueError(f'Frozen input changed: {p}')
    return c


def prepare():
    if CONTRACT.exists():
        raise FileExistsError(CONTRACT)
    upstream.verify()
    files = {upstream.EVIDENCE, old.m74.PORTABLE_MANIFEST,
             upstream.original.upstream.previous.RUNNER,
             ROOT / 'convergence_acceptance_isolation/case_ledger.csv'}
    for w in workflows():
        _, source, export, mapped = paths(w)
        cfg = old.read_json(source)
        files.update((source, source.parent/'state.csv', source.parent/'curve.csv', mapped))
        files.update(Path(cfg[k]) for k in ('mesh_file', 'node_doping_file', 'materials_file'))
        files.update(p for p in export.rglob('*.csv'))
    files.update((upstream.OUT / f'm80_{name}_ledger.csv') for name in ('case', 'component', 'spatial'))
    old.write_json(CONTRACT, {
        'status': 'frozen_before_execution', 'matrix': [{k:w[k] for k in ('device','drain_voltage_V','gate_voltage_V')} for w in workflows()],
        'new_nonlinear_solves': 0, 'new_sentaurus_runs': 0, 'production_changes': False,
        'mapping': 'Raw Si n,p,ni are used in SRH and explicit density-form SG diagnostics. Vela residual/edge/VTK probes use unchanged baseline or M80b coherent mapped state; never label their recomputed densities raw.',
        'srh_padding': 'Existing BJT probe requires all nodes. Non-Si missing rows are synthetic equilibrium placeholders n=p=ni=1, potentials=0, SRH=0 and ALWAYS excluded from summaries. Original Si rows unchanged. Common Si barycentric volume used for field comparison; other volume choices separately reported.',
        'gates': {'coordinate_um':1e-10, 'terminal_replay_relative':1e-8, 'edge_terminal_relative':1e-8,
                  'field_active_fraction':1e-6, 'srh_formula_normalized_L1_screen':0.01},
        'interpretation': 'Node current/mobility fields are reconstructed outputs, not native edge fluxes or residuals. Arithmetic nodal mobility edge substitution is a declared proxy, not a recovered Sentaurus operator. Large local differences are not causal Id attribution. Retain all strict convergence failures.',
        'input_hashes': {old.portable(p):old.sha256(p) for p in sorted(files)}})
    print(json.dumps({'status':'prepared','inputs':len(files),'cases':4}))


def probe(deck, root, name, **kwargs):
    root.mkdir(parents=True, exist_ok=True)
    dest = root/f'{name}.json'
    if dest.with_suffix('.status.json').exists():
        if old.read_json(dest) != dict(deck, **kwargs):
            raise ValueError(f'Deck changed on resume: {dest}')
        return old.read_json(dest.with_suffix('.status.json'))
    return old.run_probe(dict(deck, **kwargs), dest)


def scalar(export, name):
    return m73.scalar(export/'fields'/f'{name}_region0.csv')


def vector_metric(ref, got, fraction=1e-6):
    ref, got = np.asarray(ref), np.asarray(got)
    if ref.shape != got.shape or not np.all(np.isfinite(ref)) or not np.all(np.isfinite(got)):
        raise ValueError('Nonfinite or mismatched vectors')
    rn, gn = np.linalg.norm(ref,axis=1), np.linalg.norm(got,axis=1)
    keep = (rn >= rn.max()*fraction) & (rn > 0)
    # Zero candidate values are failures, not silently dropped from logarithmic statistics.
    logs = np.abs(np.log10(np.maximum(gn[keep], 1e-300)/rn[keep]))
    return {'nodes':int(keep.sum()),'candidate_zero_nodes':int(np.sum(gn[keep]==0)),
            'p95_dex':float(np.percentile(logs,95)),
            'vector_relative_L2':float(np.linalg.norm(got[keep]-ref[keep])/np.linalg.norm(ref[keep]))}


def srh_case(w, root, deck, export, geo):
    ids = sorted(scalar(export,'eDensity'))
    pad = root/'srh_fields'
    for name in FIELDS:
        data = scalar(export,name)
        if set(data) != set(ids): raise ValueError('Si field support mismatch')
        fill = 1. if name in ('eDensity','hDensity','EffectiveIntrinsicDensity') else 0.
        old.write_csv(pad/f'{name}_region0.csv', [{'node_id':i,'component0':data.get(i,fill)} for i in range(geo.count)])
    probe(deck,root,'srh',simulation_type='srh_fixed_state_probe',fixed_state_fields_dir=str(pad),output_csv=str(root/'srh.csv'))
    raw = {int(r['node_id']):r for r in old.read_csv(root/'srh.csv')}
    rows = [raw[i] for i in ids]
    ref = np.array([float(r['sdevice_srh_cm3_s']) for r in rows])
    records=[]
    for volume in ('barycentric_si','signed_si','all_cell'):
        # m2 -> um2; cm^-3 s^-1 * um2 * 1um -> s^-1 via 1e-12.
        area = geo.volumes[volume][ids]*1e12
        norm = np.sum(np.abs(area*ref))
        for col in RATE_COLUMNS:
            got = np.array([float(r[col]) for r in rows])
            if not np.all(np.isfinite(got)): raise ValueError('Nonfinite Si SRH')
            records.append({'device':w['device'],'vd':w['drain_voltage_V'],'volume':volume,'variant':col,
                'reference_integral_A_per_um':float(Q*1e-12*np.dot(area,ref)),
                'candidate_integral_A_per_um':float(Q*1e-12*np.dot(area,got)),
                'absolute_integral_ratio':float(np.sum(np.abs(area*got))/norm),
                'normalized_L1':float(np.sum(np.abs(area*(got-ref)))/norm)})
    old.write_csv(root/'srh_valid_silicon_nodes.csv',rows)
    identities=[]
    for name,a,b in (('ni','production_ni_eff_cm3','sdevice_ni_eff_cm3'),
                     ('electron_density','qf_reconstructed_electron_density_cm3','electron_density_cm3'),
                     ('hole_density','qf_reconstructed_hole_density_cm3','hole_density_cm3')):
        v=np.array([abs(math.log10(float(r[a])/float(r[b]))) for r in rows])
        identities.append({'device':w['device'],'vd':w['drain_voltage_V'],'field':name,'p95_dex':float(np.percentile(v,95)),'max_dex':float(v.max())})
    # Charge formula uses exact exported dopants/densities on their common Si support.
    n,p,nd,na,charge=[scalar(export,k) for k in ('eDensity','hDensity','DonorConcentration','AcceptorConcentration','SpaceCharge')]
    got=np.array([p[i]-n[i]+nd[i]-na[i] for i in ids]); refq=np.array([charge[i] for i in ids])
    # SpaceCharge exported here has concentration units; check against net number density.
    identities.append({'device':w['device'],'vd':w['drain_voltage_V'],'field':'space_charge_number_density_relative_L2',
                       'p95_dex':None,'max_dex':float(np.linalg.norm(got-refq)/np.linalg.norm(refq))})
    return records,identities


def bernoulli(x):
    x=np.asarray(x,dtype=np.longdouble)
    small=np.abs(x)<1e-5
    result=np.empty_like(x)
    result[small]=1-x[small]/2+x[small]**2/12-x[small]**4/720
    result[~small]=x[~small]/np.expm1(x[~small])
    return result


def electron_density_flux(n0,n1,eta,coef):
    # ScharfetterGummel.cpp: left=B(-eta)*n0, right=B(eta)*n1.
    return coef*(bernoulli(-eta)*n0-bernoulli(eta)*n1)


def sum_contacts(edges,mesh,nflux,pflux):
    result={}
    for contact in mesh['contacts']:
        ids=set(map(int,contact['node_ids']))
        signs=np.array([int(int(r['node0']) in ids)-int(int(r['node1']) in ids) for r in edges])
        # Continuity flux has opposite sign to ContactCurrent's conventional electron flux.
        result[contact['name']]=float(-Q*1e-6*np.sum(signs*(nflux-pflux),dtype=np.longdouble))
    return result


def edge_case(w,role,root,export,mesh,status):
    edges=old.read_csv(root/'edges.csv')
    get=lambda key:np.array([float(r[key]) for r in edges],dtype=np.longdouble)
    nf,pf=get('electron_particle_line_flux_per_m_s'),get('hole_particle_line_flux_per_m_s')
    native=sum_contacts(edges,mesh,nf,pf)
    identity=abs(native['drain']/status['current_A_per_um']-1)
    if identity>1e-8: raise ValueError(f'Edge/terminal identity failed {role}: {identity}')
    electron=get('electron_flux'); hp=get('electron_sg_high_precision_reference_flux')
    active=np.abs(electron)>np.max(np.abs(electron))*1e-6
    precision=float(np.linalg.norm(electron[active]-hp[active])/np.linalg.norm(hp[active]))
    result={'device':w['device'],'vd':w['drain_voltage_V'],'role':role,
            'drain_edge_A_per_um':native['drain'],'drain_functional_A_per_um':status['current_A_per_um'],
            'edge_terminal_relative':identity,'KCL_over_Id':abs(sum(native.values())/native['drain']),
            'electron_high_precision_relative_L2':precision}
    old.write_json(root/'contact_fluxes.json',native)
    variants=[]
    if role=='sentaurus':
        n,p,mu_n,mu_p=[scalar(export,name) for name in ('eDensity','hDensity','eMobility','hMobility')]
        ids0=np.array([int(r['node0']) for r in edges]);ids1=np.array([int(r['node1']) for r in edges])
        valid=np.array([(a in n and b in n) for a,b in zip(ids0,ids1)])
        def on_edges(field,side):return np.array([field.get(int(i),0) for i in side],dtype=np.longdouble)
        eta=get('electron_sg_eta');bm,bp=bernoulli(-eta),bernoulli(eta)
        coef=get('electron_mobility_m2_V_s')*upstream.VT*get('couple_m')/get('length_m')
        # Independent density form, long-double arithmetic; cancellation error is separately visible.
        coherent=electron_density_flux(get('electron_density0_m3'),get('electron_density1_m3'),eta,coef)
        raw=electron_density_flux(on_edges(n,ids0)*1e6,on_edges(n,ids1)*1e6,eta,coef)
        ratios=np.ones(len(edges),dtype=np.longdouble)
        good=valid & (get('electron_mobility_m2_V_s')>0)
        ratios[good]=(on_edges(mu_n,ids0)[good]+on_edges(mu_n,ids1)[good])*0.5e-4/get('electron_mobility_m2_V_s')[good]
        # Hole mobility substituted separately; no density substitution inferred for holes.
        hr=np.ones(len(edges),dtype=np.longdouble);hgood=valid & (get('hole_mobility_m2_V_s')>0)
        hr[hgood]=(on_edges(mu_p,ids0)[hgood]+on_edges(mu_p,ids1)[hgood])*0.5e-4/get('hole_mobility_m2_V_s')[hgood]
        for name,nn,pp in (('coherent_density_electron_control',coherent,pf),
                           ('raw_density_electron_only',raw,pf),
                           ('nodal_arithmetic_electron_mobility_proxy',nf*ratios,pf),
                           ('nodal_arithmetic_hole_mobility_proxy',nf,pf*hr)):
            cur=sum_contacts(edges,mesh,nn,pp)
            variants.append({'device':w['device'],'vd':w['drain_voltage_V'],'variant':name,
                'drain_A_per_um':cur['drain'],'delta_drain_A_per_um':cur['drain']-native['drain'],
                'relative_drain_change':cur['drain']/native['drain']-1,'KCL_over_Id':abs(sum(cur.values())/cur['drain'])})
        old.write_csv(root/'edge_substitution.csv',[{'edge_id':int(r['edge_id']),'node0':int(ids0[i]),'node1':int(ids1[i]),
            'original_electron_line_flux':float(nf[i]),'coherent_density_line_flux':float(coherent[i]),
            'raw_density_line_flux':float(raw[i]),'electron_nodal_mobility_ratio':float(ratios[i]),
            'hole_nodal_mobility_ratio':float(hr[i])} for i,r in enumerate(edges)])
    return result,variants


def run():
    verify_inputs()
    if EVIDENCE.exists(): raise FileExistsError(EVIDENCE)
    old.m74.RUNNER=upstream.original.upstream.previous.RUNNER
    srh=[];identities=[];vectors=[];terminals=[];substitutions=[];mobility=[]
    for w in workflows():
        tag,source,export,mapped=paths(w);root=LOCAL/tag
        geo=m73.Geometry(w['device'])
        if m73.coordinate_error(export,geo)>1e-10: raise ValueError('Coordinate identity')
        base=old.base_deck(source,source.parent/'state.csv',w)
        s,i=srh_case(w,root,base,export,geo);srh.extend(s);identities.extend(i)
        mesh=old.read_json(Path(base['mesh_file']))
        for role,state in (('baseline',source.parent/'state.csv'),('sentaurus',mapped)):
            dest=root/role;deck=old.base_deck(source,state,w)
            status=probe(deck,dest,'functional',simulation_type='terminal_current_functional_probe')
            probe(deck,dest,'vtk',simulation_type='write_dd_state_vtk',output_vtk=str(dest/'state.vtk'),output_diagnostics={'cell_first_sg_current_recovery':True})
            probe(deck,dest,'edges',simulation_type='sg_edge_flux_probe',output_csv=str(dest/'edges.csv'))
            probe(deck,dest,'mobility',simulation_type='edge_mobility_probe',output_csv=str(dest/'mobility.csv'))
            er,sub=edge_case(w,role,dest,export,mesh,status);terminals.append(er);substitutions.extend(sub)
            count,sc,vec=read_vtk_point_data(dest/'state.vtk')
            for carrier,prefix in (('electron','e'),('hole','h')):
                rr={int(r['node_id']):(float(r['component0']),float(r['component1'])) for r in old.read_csv(export/'fields'/f'{prefix}CurrentDensity_region0.csv')}
                ids=sorted(rr);ref=np.array([rr[i] for i in ids])
                for method,key in (('direct',f'DualFaceSg{carrier.title()}CurrentDensityVector'),('cell_first',f'CellFirstSg{carrier.title()}CurrentDensityVector')):
                    got=np.array([vec[key][i][:2] for i in ids])
                    vectors.append({'device':w['device'],'vd':w['drain_voltage_V'],'role':role,'carrier':carrier,'method':method,**vector_metric(ref,got)})
            for carrier in ('electron','hole'):
                rr=old.read_csv(dest/'mobility.csv')
                errors=[abs(float(r[f'{carrier}_limiter_reconstruction_error'])) for r in rr]
                mobility.append({'device':w['device'],'vd':w['drain_voltage_V'],'role':role,'carrier':carrier,'max_limiter_reconstruction_error':max(errors)})
        print(f'Completed {tag}: SRH, 2 state replays, 4 current-field comparisons',flush=True)
    for name,rows in (('srh',srh),('identities',identities),('vectors',vectors),('terminals',terminals),('substitutions',substitutions),('mobility',mobility)):
        old.write_csv(OUT/f'{name}.csv',rows)
    # Preserve the already qualified full-DD equation and spatial ledger, restricted to this matrix.
    for name in ('case','component','spatial'):
        rows=[r for r in old.read_csv(upstream.OUT/f'm80_{name}_ledger.csv') if r['device'] in ('n19','n23')]
        old.write_csv(OUT/f'dd_{name}_ledger.csv',rows)
    strict=old.read_csv(ROOT/'convergence_acceptance_isolation/case_ledger.csv')
    old.write_csv(OUT/'strict_qualification.csv',strict)
    summary={'status':'fixed_state_comparisons_completed_causal_correction_not_qualified','cases':4,
             'new_read_only_probes':36,'new_nonlinear_solves':0,'new_sentaurus_runs':0,
             'maximum_edge_terminal_relative':max(r['edge_terminal_relative'] for r in terminals),
             'maximum_electron_high_precision_relative_L2':max(r['electron_high_precision_relative_L2'] for r in terminals),
             'production_changes':False,'m82_released':False,
             'limitations':['Mapped state is M80b coherent Vela representation, not raw n/p.',
                            'Non-Si placeholders in BJT SRH input are excluded from analysis.',
                            'Raw-density subtraction and node-mobility proxy are diagnostics, not self-consistent solutions.',
                            'Strict qualification remains 1/4; small current movement is not a convergence error bound.']}
    old.write_json(OUT/'summary.json',summary)
    verify_inputs()
    print(json.dumps(summary,indent=2))


def seal():
    verify_inputs()
    files=[SCRIPT,CONTRACT,REPO/'tests/regression/test_simplemos_fixed_state_formula_audit.py',
           REPO/'docs/validation/simplemos_fixed_state_formula_audit_2026-09-05.md']
    files += [p for folder in (LOCAL,OUT) for p in folder.rglob('*') if p.is_file() and p!=EVIDENCE]
    old.write_json(EVIDENCE,{'status':'frozen','hashes':{old.portable(p):old.sha256(p) for p in sorted(set(files))}})


def thermal_control():
    """Test raw-density hybrid artifacts with the already inferred source constant."""
    verify_inputs()
    source=upstream.OUT/'m80b_source_identity.csv'
    contract=OUT/'thermal_control_contract.json'
    if not contract.exists():
        old.write_json(contract,{'status':'frozen_before_execution','source':old.portable(source),
            'source_sha256':old.sha256(source),'script_sha256':old.sha256(SCRIPT),
            'purpose':'Only replace Vt in the independent raw-density SG formula by the previously recorded source Vt; geometry and Vela mobility remain fixed. No solver/material change and no Id fit.',
            'matrix':old.read_json(CONTRACT)['matrix']})
    c=old.read_json(contract)
    if old.sha256(source)!=c['source_sha256']: raise ValueError('Thermal identity changed')
    source_rows=old.read_csv(source);rows=[]
    for w in workflows():
        tag,sourcecfg,export,mapped=paths(w)
        root=LOCAL/tag/'sentaurus'
        edges=old.read_csv(root/'edges.csv')
        mesh=old.read_json(Path(old.read_json(sourcecfg)['mesh_file']))
        vt=np.longdouble(next(r['source_inferred_thermal_voltage_V'] for r in source_rows if r['device']==w['device'] and float(r['drain_voltage_V'])==float(w['drain_voltage_V'])))
        get=lambda k:np.array([r[k] for r in edges],dtype=np.longdouble)
        n=scalar(export,'eDensity')
        n0=np.array([n.get(int(r['node0']),0)*1e6 for r in edges],dtype=np.longdouble)
        n1=np.array([n.get(int(r['node1']),0)*1e6 for r in edges],dtype=np.longdouble)
        eta=(get('psi1_V')-get('psi0_V'))/vt
        coef=get('electron_mobility_m2_V_s')*vt*get('couple_m')/get('length_m')
        nf=electron_density_flux(n0,n1,eta,coef)
        current=sum_contacts(edges,mesh,nf,get('hole_particle_line_flux_per_m_s'))
        base=old.read_json(root/'contact_fluxes.json')['drain']
        rows.append({'device':w['device'],'vd':w['drain_voltage_V'],'source_Vt_V':float(vt),
                     'drain_A_per_um':current['drain'],'relative_drain_change':current['drain']/base-1,
                     'KCL_over_Id':abs(sum(current.values())/current['drain'])})
    old.write_csv(OUT/'thermal_control.csv',rows)
    print(json.dumps(rows,indent=2))


def verify():
    verify_inputs()
    for p,digest in old.read_json(EVIDENCE)['hashes'].items():
        if old.sha256(REPO/p)!=digest: raise ValueError(f'Evidence changed: {p}')
    print('Fixed-state evidence verified')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run','thermal-control','seal','verify'))
    {'prepare':prepare,'run':run,'thermal-control':thermal_control,'seal':seal,'verify':verify}[p.parse_args().action]()
