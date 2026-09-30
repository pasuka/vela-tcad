"""Portable original 816-target mesh/doping/reference inventory; no solve."""
from pathlib import Path
import copy,hashlib,json,shutil,tarfile
import simplemos_hfs_cloud_20260926 as h
R=Path(__file__).resolve().parents[1]
O=R/'build/hfs_srh_20260926/original_inputs'
M=R/'build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics'
def main():
    assert not O.exists();O.mkdir(parents=True)
    old=h.read(M/'vela/vela_matrix_manifest.json');ref=h.read(M/'references/reference_manifest.json')
    assert ref['status']=='qualified';cases=[];scope=[]
    template=h.read(R/'build/hfs_merged_20260926/payload_v5/m65_n23_vd_1p000000_endpoint/template.json')
    material=R/'build/hfs_merged_20260926/payload_v5/m65_n23_vd_1p000000_endpoint/materials_file.json'
    shutil.copy2(material,O/'materials.json')
    for x in old['inputs']:
        dev=x['device'];p=O/dev;p.mkdir()
        oldcfg=h.read(M/'vela'/dev/'base.json')
        mesh=Path(oldcfg['mesh_file']);doping=Path(oldcfg['node_doping_file'])
        assert h.sha(mesh)==x['mesh_sha256'];assert h.sha(doping)==x['doping_sha256']
        shutil.copy2(mesh,p/'mesh.json');shutil.copy2(doping,p/'doping.csv')
        m=h.read(mesh);si={r['id'] for r in m['regions'] if r['material']=='Si'};ids=set();contacts={i for c in m['contacts'] for i in c['node_ids']}
        for t in m['triangles']:
            if t['region_id'] in si:ids.update(t['node_ids'])
        geo=dict(count=len(m['nodes']),free_si=sorted(ids-contacts),all_si=sorted(ids),contacts=sorted(contacts))
        h.write(p/'geometry.json',geo);scope.append(dict(device=dev,nodes=geo['count'],carrier_rows=2*len(geo['free_si'])))
        cfg=copy.deepcopy(template);cfg.pop('simplemos_m65',None)
        cfg.update(mesh_file=f'{dev}/mesh.json',node_doping_file=f'{dev}/doping.csv',materials_file='materials.json')
        cfg['solver']['region_resolved_interface_assembly']['srh_signed_transport_volume_fraction']=1.
        for c in cfg['contacts']:c['bias']=0.
        h.write(p/'template.json',cfg)
        for rr in ref['artifacts']:
            if rr['device']!=dev:continue
            f=M/'references'/rr['path'];assert h.sha(f)==rr['sha256'];points=h.rows(f)
            assert len(points)==51 and all(abs(float(r['gate_voltage_V'])-i*.05)<1e-12 for i,r in enumerate(points))
            shutil.copy2(f,p/rr['path']);cases.append(dict(case=rr['case'],device=dev,vd=rr['drain_voltage_V'],reference=f'{dev}/{rr["path"]}',points=51))
    assert len(cases)==16
    h.write(O/'contract.json',dict(cases=cases,scope=scope,targets=816,reference='Frozen M8 original-deck native curves; no interpolation. Keep distinct from strengthened HFS 204-point references.',source_manifest_sha256=h.sha(M/'references/reference_manifest.json'),mesh_manifest_sha256=h.sha(M/'vela/vela_matrix_manifest.json'),launch_gate='SRH fixed-state, derivative, two-amplitude, 16-point dual and field checks must pass first.',current_relative_percent=2.,cold_start='Poisson block -> coupled equilibrium -> drain ramp -> gate sweep; only accepted predecessors may propagate.'))
    h.write(O/'hashes.json',{p.relative_to(O).as_posix():h.sha(p) for p in O.rglob('*') if p.is_file()})
    with tarfile.open(O.parent/'original_inputs.tgz','w:gz') as t:t.add(O,arcname='original_inputs')
    print(json.dumps(dict(cases=len(cases),targets=816,scope=scope,archive_bytes=(O.parent/'original_inputs.tgz').stat().st_size)))
if __name__=='__main__':main()
