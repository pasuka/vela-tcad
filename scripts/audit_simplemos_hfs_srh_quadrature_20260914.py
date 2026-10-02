"""Audit SRH quadrature with signed carrier terminal sums at sixteen HFS points.
No candidate solve or source-volume modification. Keep native and Vela q distinct.
"""
import math
from pathlib import Path
import numpy as np
import simplemos_hfs_curves_20260914 as e
a,d=e.a,e.d
QVELA=1.602176634e-19
QNATIVE=1.602192e-19

def main():
    for p in (e.O/'sixteen_completion_evidence.json',e.O/'native_evidence.json',e.O/'low_fields/evidence.json',e.q.O/'fields/evidence.json'):a.verify(p)
    outputs=[];cells=[];files=[Path(__file__).resolve(),e.O/'sixteen_completion_evidence.json',e.O/'native_evidence.json',e.O/'low_fields/evidence.json',e.q.O/'fields/evidence.json']
    state_nodes=a.rows(e.O/'low_fields/nodes.csv')+a.rows(e.q.O/'fields/nodes.csv')
    early=a.rows(e.O/'early_low_attempts.csv');high=a.rows(e.q.O/'selected_states.csv')
    for cc in a.read(e.O/'vela_contract.json')['cases']:
        geo=d.matrix.spatial.m73.Geometry(cc['device'])
        for i in (0,10,40,50):
            key=cc['case']+f'_vg_{i:03d}';rows=[r for r in state_nodes if r['key']==key];assert rows
            ids=np.array([int(r['node_id']) for r in rows]);nr=np.array([float(r['srhRecombination_native']) for r in rows]);vr=np.array([float(r['srhRecombination_vela']) for r in rows])
            path=e.L/'native_raw/bundle'/cc['case']/f'point_{i:03d}_native_des.plt';ports=e.n.exporter.pltrows(path)[0];files.append(path)
            sums={car:math.fsum(ports[n+' '+car+'Current'] for n in ('drain','source','substrate','gate')) for car in ('e','h')}
            conditions={car:math.fsum(abs(ports[n+' '+car+'Current']) for n in ('drain','source','substrate','gate'))/max(abs(sums[car]),1e-300) for car in ('e','h')}
            srh={name:QNATIVE*math.fsum(np.asarray(vol)[ids]*nr) for name,vol in geo.volumes.items()}
            vela_all=QVELA*math.fsum(np.asarray(geo.volumes['all_cell'])[ids]*vr)
            vela_signed=QVELA*math.fsum(np.asarray(geo.volumes['signed_si'])[ids]*vr)
            found=next(r for r in (early if i<40 else high) if r['case']==cc['case'] and int(r['index'])==i and r['arm']=='enormal_seed' and r['qualified']=='True')
            root=Path(found['dest']);cfg=a.read(root/'config.json');mesh=a.read(Path(cfg['mesh_file']));edges=a.rows(root/'acceptance_edges.csv');vs=a.read(root/'config.status.json')['contact_currents_A_per_um'];vports={}
            for contact in mesh['contacts']:
                cn=set(contact['node_ids']);raw={car:QVELA*1e-6*math.fsum((int(int(r['node0']) in cn)-int(int(r['node1']) in cn))*float(r[car+'_particle_line_flux_per_m_s']) for r in edges) for car in ('electron','hole')}
                vports[contact['name']]=dict(e=-raw['electron'],h=raw['hole'])
                assert abs((-raw['electron']+raw['hole'])-vs[contact['name']])<=1e-12*max(abs(vs[contact['name']]),1e-30)
            vh=math.fsum(x['h'] for x in vports.values())
            delta_Id=vs['drain']-ports['drain TotalCurrent']
            delta_sources=-(vela_all-sums['h'])
            out=dict(key=key,device=cc['device'],vd=cc['vd'],vg=e.n.GRID[i],native_e_sum_A_per_um=sums['e'],native_h_sum_A_per_um=sums['h'],native_e_sum_condition=conditions['e'],native_h_sum_condition=conditions['h'],native_node_SRH_signed_Si_A_per_um=srh['signed_si'],native_node_SRH_all_cell_A_per_um=srh['all_cell'],native_node_SRH_barycentric_Si_A_per_um=srh['barycentric_si'],native_signed_to_h_relative=abs(srh['signed_si']/sums['h']-1),native_all_to_h_relative=abs(srh['all_cell']/sums['h']-1),vela_actual_SRH_A_per_um=vela_all,vela_h_sum_A_per_um=vh,vela_source_to_h_relative=abs(vela_all/vh-1),vela_fixed_state_signed_Si_SRH_A_per_um=vela_signed,delta_generation_A_per_um=delta_sources,delta_Id_A_per_um=delta_Id,source_difference_over_delta_Id=delta_sources/delta_Id if delta_Id else None,delta_source_total_A_per_um=vs['source']-ports['source TotalCurrent'],delta_body_total_A_per_um=vs['substrate']-ports['substrate TotalCurrent'])
            outputs.append(out)
            for node,nrate,vrate in zip(ids,nr,vr):
                diff=QVELA*vrate*(geo.volumes['signed_si'][node]-geo.volumes['all_cell'][node])
                cells.append(dict(key=key,node_id=int(node),x_um=float(geo.coords[node][0]),y_um=float(geo.coords[node][1]),native_rate_cm3_s=float(nrate),vela_rate_cm3_s=float(vrate),signed_volume_m2=float(geo.volumes['signed_si'][node]),all_cell_volume_m2=float(geo.volumes['all_cell'][node]),fixed_state_source_volume_delta_A_per_um=float(diff)))
            files += [root/'config.json',root/'config.status.json',root/'acceptance_edges.csv']
            print(key,'signed/native hole',out['native_signed_to_h_relative'],'all/native hole',out['native_all_to_h_relative'],flush=True)
    a.write_csv(e.O/'srh_volume_audit.csv',outputs);a.write_csv(e.O/'srh_volume_node_contributions.csv',cells)
    summary=dict(points=len(outputs),max_native_signed_h_relative=max(r['native_signed_to_h_relative'] for r in outputs),max_native_h_sum_condition=max(r['native_h_sum_condition'] for r in outputs),max_vela_source_h_relative=max(r['vela_source_to_h_relative'] for r in outputs),q_native_C=QNATIVE,q_vela_C=QVELA,interpretation='Existing signed Si quadrature of native nodal SRH compared with signed native hole-port sum. Electron sum cancellation is recorded and not used as the primary strong-inversion source witness. Fixed-state volume differences do not predict a self-consistent drain response without a sensitivity calibration.',production_changed=False,self_consistent_candidate_run=False,acceptance_changed=False)
    a.write(e.O/'srh_volume_summary.json',summary)
    d.matrix.freeze(e.O/'srh_volume_evidence.json',files+[e.O/'srh_volume_audit.csv',e.O/'srh_volume_node_contributions.csv',e.O/'srh_volume_summary.json'])
    print(summary,flush=True)

if __name__=='__main__':main()
