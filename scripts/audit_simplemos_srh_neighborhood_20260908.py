"""Native/Vela fixed-state source support ledger on frozen hotspot rings."""
import argparse
import math
from pathlib import Path
import numpy as np
import audit_simplemos_linesearch_precision_20260908 as audit
import audit_simplemos_masetti_box_mobility_20260907 as native

t=audit.t; c=t.c; a=t.a; d=t.d; p=t.p; v=t.v
OUT=audit.OUT/'neighborhood'


def prepare():
    a.verify(t.OUT.parent/'final_evidence.json');jobs=[]
    files=[Path(__file__).resolve(),Path(native.__file__),t.OUT.parent/'final_evidence.json',audit.prior.OUT/'validation_evidence.json']
    for case in a.read(t.OUT/'contract.json')['cases']:
        row=next(x for x in a.rows(audit.prior.OUT/'selected_states.csv') if x['key']==case['key'] and x['axis']=='joint' and x['label']=='finite')
        assert row['qualified']=='True';dest=p.REPO/row['path'];geo,mask=v.m.previous.prior.support(case)
        roots=([793,794] if case['device']=='n19' else [794,1056]) if case['vd']==1. else ([1086] if case['device']=='n19' else [1087])
        anchors=[case['mapped_nodes'][tag]['node'] for tag in ('792','1057')]
        ring=set(roots+anchors)
        for edge in a.rows(dest/'acceptance_edges.csv'):
            i,j=int(edge['node0']),int(edge['node1'])
            if i in roots or j in roots:ring.update((i,j))
        ring=sorted(i for i in ring if mask[i])
        jobs.append(dict(case=case,path=str(dest),roots=roots,anchors=anchors,nodes=ring))
        files += [dest/x for x in ('state.csv','all_row.csv','acceptance_edges.csv','independent_acceptance.json')]
        files += list((c.s.prior.old.l.raw_root(case)/'fields').glob('*_region0.csv'))
    a.write(OUT/'contract.json',dict(jobs=jobs,
        source='Independent native element-box Masetti/SG replay in both qf and density forms; signed Si versus current original-volume SRH. Selected state is qualified joint finite arm.',
        roots='Previously reported worst free-Si phin/phip nodes, one edge ring; preserve original two-node anchors.',
        support_diagnostic_gates=dict(minimum_native_source_over_64eps_absolute_flux=100,maximum_qf_density_net_disagreement_over_source=1e-4,maximum_signed_Si_balance=1e-5),
        boundary='Fixed exported-state support audit, not native same-source DC response or permission for finite volume expansion. Retain cancellation-limited and no-volume-change controls. No fitted volumes.'))
    d.matrix.freeze(OUT/'freeze.json',files+[OUT/'contract.json']);print('Frozen neighborhood sizes',[len(j['nodes']) for j in jobs],flush=True)


def run():
    a.verify(OUT/'freeze.json');scope=a.read(OUT/'contract.json');rows=[];nodesout=[];checks=[]
    geometries={dev:native.geometry(dev) for dev in ('n19','n23')}
    for geo,weights,metrics in geometries.values():assert all(x['qualified'] for x in metrics)
    for job in scope['jobs']:
        case=job['case'];dev=case['device'];geo,weights,_=geometries[dev];vol=c.s.prior.old.l.g.geometry(dev)[3]
        dest=Path(job['path']);raw=c.native_fields(case);state=d.ordered(dest/'state.csv',geo.count);terms=d.ordered(dest/'all_row.csv',geo.count)
        edges=a.rows(dest/'acceptance_edges.csv');conversion=next(float(e['electron_particle_line_flux_per_m_s'])/float(e['electron_flux']) for e in edges if float(e['electron_flux'])!=0.)
        signed={i:float(vol['Si'][i]) for i in job['nodes']};original={i:float(geo.volumes['all_cell'][i]) for i in job['nodes']}
        selected_ratios={rec['node']:rec['ratio'] for rec in case['mapped_nodes'].values()}
        currents={mode:{car:{i:[] for i in job['nodes']} for car in ('electron','hole')} for mode in ('qf','density')}
        for e in edges:
            i,j=int(e['node0']),int(e['node1'])
            if not ({i,j}&set(job['nodes'])) or i not in raw['electrons_m3'] or j not in raw['electrons_m3']:continue
            for mode in currents:
                en,hp=native.sg_flux(raw['psi'][i],raw['psi'][j],raw['electrons_m3'][i],raw['electrons_m3'][j],raw['holes_m3'][i],raw['holes_m3'][j],raw['phin'][i],raw['phin'][j],raw['phip'][i],raw['phip'][j],mode)
                w=weights[tuple(sorted((i,j)))]
                for carrier,value in [('electron',float(en)*w['e']),('hole',float(hp)*w['h'])]:
                    if i in currents[mode][carrier]:currents[mode][carrier][i].append(value)
                    if j in currents[mode][carrier]:currents[mode][carrier][j].append(-value)
        for i in job['nodes']:
            curvol=original[i]*selected_ratios.get(i,1.);native_rate=float(raw['SRH'][i]);nativesource=native_rate*signed[i]*1e4
            vela_rate=float(terms[i]['electron_recombination'])*conversion/(curvol*1e6)
            nodesout.append(dict(key=case['key'],device=dev,vd=case['vd'],node=i,role='root' if i in job['roots'] else 'anchor' if i in job['anchors'] else 'ring',x_um=geo.coords[i][0],y_um=geo.coords[i][1],original_volume_m2=original[i],current_volume_m2=curvol,signed_Si_volume_m2=signed[i],Si_over_current_volume=signed[i]/curvol,
                psi_error_V=float(v.m.previous.prior.physical(state[i],'psi')-raw['psi'][i]),phin_error_V=float(v.m.previous.prior.physical(state[i],'phin')-raw['phin'][i]),phip_error_V=float(v.m.previous.prior.physical(state[i],'phip')-raw['phip'][i]),native_SRH_cm3_s=native_rate,Vela_SRH_cm3_s=vela_rate,SRH_rate_relative=vela_rate/native_rate-1 if native_rate else 'zero_native'))
            for carrier in ('electron','hole'):
                sums={mode:math.fsum(currents[mode][carrier][i]) for mode in currents}
                disagreement=abs(sums['qf']-sums['density'])/max(abs(nativesource),1e-300)
                for mode in currents:
                    total=sums[mode];absolute=math.fsum(abs(x) for x in currents[mode][carrier][i]);scale=max(abs(total),abs(nativesource),1e-300)
                    signal=abs(nativesource)/max(64*np.finfo(float).eps*absolute,1e-300)
                    native_balance=abs(total+nativesource)/scale
                    qualified=signal>=100 and disagreement<=1e-4 and native_balance<=1e-5
                    r=terms[i];vela_scale=max(float(r[carrier+'_flux_abs_sum']),abs(float(r[carrier+'_recombination'])))
                    rows.append(dict(key=case['key'],device=dev,vd=case['vd'],node=i,role=nodesout[-1]['role'],carrier=carrier,mode=mode,Si_over_current_volume=signed[i]/curvol,
                        net_native_flux=total,absolute_native_flux=absolute,native_Si_source=nativesource,native_current_volume_source=native_rate*curvol*1e4,
                        signed_Si_balance=native_balance,current_volume_balance=abs(total+native_rate*curvol*1e4)/scale,
                        source_over_64eps_flow=signal,qf_density_disagreement_over_source=disagreement,source_support_resolved=qualified,
                        Vela_row_ratio=abs(float(r[carrier+'_residual']))/vela_scale,
                        Vela_source_change_internal=float(r[carrier+'_recombination'])*(signed[i]/curvol-1.),
                        native_source_response_calibrated=False))
        print(case['key'],'nodes',len(job['nodes']),flush=True)
    a.write_csv(OUT/'nodes.csv',nodesout);a.write_csv(OUT/'rows.csv',rows)
    for job in scope['jobs']:
        selected=[r for r in rows if r['key']==job['case']['key'] and r['mode']=='qf']
        checks.append(dict(key=job['case']['key'],nodes=len(job['nodes']),rows=len(selected),source_resolved=sum(r['source_support_resolved'] for r in selected),unresolved=sum(not r['source_support_resolved'] for r in selected),max_Vela_row_ratio=max(r['Vela_row_ratio'] for r in selected)))
    a.write_csv(OUT/'summary.csv',checks);d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',OUT/'nodes.csv',OUT/'rows.csv',OUT/'summary.csv']);print(checks,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','run'));globals()[parser.parse_args().action]()
