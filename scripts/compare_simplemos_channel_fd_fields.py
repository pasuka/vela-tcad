"""Compare matched node field responses after signed native/Vela DC solves."""
import math
from pathlib import Path
import numpy as np
import validate_simplemos_vela_channel_charge as run


def compare():
    run.audit.verify(run.FREEZE)
    result = run.audit.read(run.OUT / 'result.json')
    if result['qualified_states'] != 5:
        raise ValueError('Five qualified Vela states required before field interpretation')
    m73 = run.native.upstream.m73
    mask = {int(r['node_id']): r['selected']=='True' for r in run.audit.rows(run.native.OUT/'node_mask.csv')}
    volumes, _ = run.native.upstream.native_arrays()
    vf = {}; sf = {}
    for name, _ in run.CASES:
        state = {int(r['node_id']): r for r in run.audit.rows(run.LOCAL/name/'state.csv')}
        vf[name] = {k: {i: float(r[field])*scale for i,r in state.items()} for k,field,scale in
                    (('Potential','psi',1.),('eDensity','electrons_m3',1e-6),('hDensity','holes_m3',1e-6))}
        sf[name] = {k: m73.scalar(run.native.LOCAL/'exports'/name/'fields'/f'{"ElectrostaticPotential" if k=="Potential" else k}_region0.csv') for k in vf[name]}
    detail = []; summaries = []
    for suffix, amplitude in (('1e13',1e13),('5e12',5e12)):
        block = []
        for i in sorted(sf['zero']['Potential']):
            r = {'amplitude_cm_3': amplitude, 'node_id': i, 'channel': mask[i], 'signed_Si_area_um2': volumes[i]}
            for key,label in (('Potential','delta_psi_V'),('eDensity','delta_n_over_n0'),('hDensity','delta_p_over_p0')):
                for code, fields in (('vela',vf),('native',sf)):
                    scale = fields['zero'][key][i] if key!='Potential' else 1.
                    value = (fields['plus_'+suffix][key][i]-fields['minus_'+suffix][key][i])/2/scale
                    assert math.isfinite(value)
                    r[code+'_'+label] = value
            detail.append(r); block.append(r)
        for support in ('channel','outside_channel','all_si'):
            selected = [r for r in block if support=='all_si' or r['channel']==(support=='channel')]
            w = np.array([abs(r['signed_Si_area_um2']) for r in selected])
            for field in ('delta_psi_V','delta_n_over_n0','delta_p_over_p0'):
                v = np.array([r['vela_'+field] for r in selected]); s = np.array([r['native_'+field] for r in selected])
                norm = float(np.dot(w,s*s))
                if norm==0: raise ValueError('Null native field response')
                summaries.append({'amplitude_cm_3':amplitude,'support':support,'field':field,'nodes':len(selected),
                    'native_area_weighted_rms':float(np.sqrt(norm/w.sum())),
                    'vela_area_weighted_rms':float(np.sqrt(np.dot(w,v*v)/w.sum())),
                    'relative_weighted_L2_difference':float(np.sqrt(np.dot(w,(v-s)**2)/norm)),
                    'projection_vela_on_native':float(np.dot(w,v*s)/norm)})
    run.audit.write_csv(run.OUT/'field_response_nodes.csv',detail)
    run.audit.write_csv(run.OUT/'field_response_summary.csv',summaries)
    for r in summaries:
        if r['support']=='channel': print(r)


if __name__=='__main__':
    compare()
