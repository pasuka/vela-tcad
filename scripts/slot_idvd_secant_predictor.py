"""Diagnostic Fermi IdVd secant seed; never a qualified device state.

Physical potentials are extrapolated in volts. Reference/increment QF pairs
are extrapolated separately. Positive density guesses use log extrapolation,
not a Boltzmann occupation law; the unchanged Fermi Newton assembler evaluates
the constitutive densities. No production continuation policy is enabled here.
"""
from pathlib import Path
import hashlib
import numpy as np
import slot_state_archive as states


def extrapolate_fields(older, parent, ratio):
    if not np.isfinite(ratio) or not 0 < ratio <= 3:
        raise ValueError('Secant ratio must lie in (0,3]')
    if older.keys() != parent.keys() or not states.archive.SPLIT <= parent.keys():
        raise ValueError('Matching referenced QF fields required')
    if any(k in parent for k in states.archive.PACKED + states.archive.PACKED_LOW):
        raise ValueError('Packed-coordinate prediction is outside this diagnostic')
    out={k:np.array(v,copy=True) for k,v in parent.items()}
    for k in ('psi', *sorted(states.archive.SPLIT)):
        if older[k].shape != parent[k].shape:raise ValueError('Field shape mismatch')
        out[k]=np.asarray(parent[k],dtype=np.longdouble)+np.longdouble(ratio)*(
            np.asarray(parent[k],dtype=np.longdouble)-np.asarray(older[k],dtype=np.longdouble))
        out[k]=out[k].astype(float)
    for carrier,q in (('electron','phin'),('hole','phip')):
        out[q]=out[carrier+'_qf_reference_V']+out[carrier+'_qf_increment_V']
    for k in ('electrons_m3','holes_m3'):
        if np.any(older[k]<0) or np.any(parent[k]<0):raise ValueError('Negative density')
        mask=(older[k]>0)&(parent[k]>0)
        out[k][mask]=np.exp(np.log(parent[k][mask])+ratio*(np.log(parent[k][mask])-np.log(older[k][mask])))
        # A zero-density field has no logarithmic secant; retain the parent.
    if any(not np.isfinite(v).all() for v in out.values()):raise ValueError('Nonfinite prediction')
    return out


def predict(config, older_path, parent_path, destination, root, target):
    if config['solver'].get('carrier_statistics')!='fermi_dirac' or config['sweep']['contact']!='drain':
        raise ValueError('Only explicit Fermi drain diagnostics supported')
    old,om=states.read(older_path,config,root);parent,pm=states.read(parent_path,config,root)
    for fields,meta in ((old,om),(parent,pm)):
        if meta.get('state_role')!='accepted' or meta.get('qualified') is False:
            raise ValueError('Both historical states must be accepted')
        for key in ('mesh_file','materials_file','node_doping_file'):
            expected=hashlib.sha256((Path(root)/config[key]).read_bytes()).hexdigest()
            if meta['input_file_sha256'][key]!=expected:raise ValueError('Input identity mismatch')
        if 'electron_quantum_potential_V' in fields and np.any(fields['electron_quantum_potential_V']!=0):
            raise ValueError('Quantum predictor not validated')
    ob,pb=om['contact_biases_V'],pm['contact_biases_V']
    if any(ob[k]!=pb[k] for k in pb if k!='drain'):raise ValueError('Non-drain bias changed')
    step=pb['drain']-ob['drain']
    if step<=0:raise ValueError('Non-increasing drain history')
    ratio=(target-pb['drain'])/step
    output=extrapolate_fields(old,parent,ratio)
    states.write_candidate(destination,output,pm,config,[older_path,parent_path],target,'idl_idvd_secant_diagnostic')
    return dict(qualified=False,ratio=ratio,previous_V=ob['drain'],parent_V=pb['drain'],target_V=target,
        density_guess='log_secant_not_Boltzmann_law',
        max_change_V={k:float(np.max(np.abs(output[k]-parent[k]))) for k in ('psi','phin','phip')})
