"""Validate and package the read-only investigation for independent review."""
import csv,hashlib,json,math,shutil,zipfile
from pathlib import Path
import h5py
ROOT=Path(__file__).resolve().parents[1]
B=ROOT/'build/outlier_analysis_20260928'
O=B/'analysis'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rows(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def main():
    cloud=json.loads((B/'cloud_hashes.json').read_text())
    for p,s in cloud.items():assert sha(B/p)==s,p
    points=rows(B/'outliers.csv');nodes=rows(B/'detailed/outlier_node_fields.csv')
    assert len(nodes)==9066
    source={(r['case'],int(r['index'])):float(r['free_Si_SRH_A_per_um']) for r in rows(B/'srh_source_integrals.csv')}
    for r in points:
        case=r['case'];idx=int(r['index']);group=[n for n in nodes if n['case']==case and int(n['index'])==idx]
        mesh=json.loads((B/f"inputs/{r['device']}/mesh.json").read_text())
        assert len(group)==len(mesh['nodes']) and [int(n['node_id']) for n in group]==list(range(len(group)))
        with h5py.File(B/f'detailed/{case}_vg_{idx:03d}.h5') as f:
            metadata=json.loads(f.attrs['metadata_json'])
            assert abs(metadata['bias_V']-float(r['vg']))<1e-10
            for name,data in f['fields'].items():
                assert all(float(n[name])==float(v) for n,v in zip(group,data[()]))
            for n,m in zip(group,mesh['nodes']):assert float(n['x_um'])==m['x'] and float(n['y_um'])==m['y']
            for k,v in metadata['input_file_sha256'].items():
                path={'materials_file':B/'inputs/materials.json','mesh_file':B/f"inputs/{r['device']}/mesh.json",'node_doping_file':B/f"inputs/{r['device']}/doping.csv"}[k]
                assert sha(path)==v
        integral=math.fsum(float(n['SRH_integrated_A_per_um']) for n in group if n['free_si']=='True')
        assert math.isclose(integral,source[case,idx],rel_tol=1e-12)
    # Local source is included only if it is byte-identical to the frozen cloud source.
    frozen=json.loads((B/'frozen/source_hashes.json').read_text());copied=[];unavailable=[]
    wanted=['src/equation/CoupledDDAssembler.cpp','src/equation/SplitDDRuntime.cpp',
        'include/vela/numerics/SplitDDState.h','include/vela/equation/SplitDDOperator.h','include/vela/discretization/ScharfetterGummel.h',
        'include/vela/discretization/StableSGDerivative.h','include/vela/discretization/ElementQfGradient.h',
        'src/physics/MobilityModel.cpp','src/physics/RecombinationModel.cpp',
        'include/vela/core/PhysicalConstants.h','include/vela/core/UnitScaling.h',
        'include/vela/core/UnitScalingSystem.h','scripts/simplemos_srh_cloud_20260926.py']
    for name in wanted:
        p=ROOT/name
        if p.exists() and name in frozen and sha(p)==frozen[name]:
            dest=O/'verified_source'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest);copied.append(name)
        else:unavailable.append(name)
    for name in ['collect_simplemos_outliers_20260928.py','analyze_simplemos_outliers_20260928.py','package_simplemos_outliers_20260928.py','sentaurus_import.py']:
        d=O/'analysis_scripts'/name;d.parent.mkdir(exist_ok=True);shutil.copy2(ROOT/'scripts'/name,d)
    report=(ROOT/'docs/validation/simplemos_outlier_analysis_2026-09-28.md').read_text(encoding='utf-8')
    prompt='''# 给协助定位的大模型

请分析以下 SimpleMOS 六个低 Vd 电流超限点，依据实测数据提出有区分力的根因诊断。
请先读后附报告，再读核心包 CSV；缺少原生同点空间场和日志，请明确写出无法判定的事项。
你的输出需要：①按证据排序的原因候选及反证；②端口电子/空穴分量和 KCL 的独立复算；③说明为什么衬底电子分量差异可以与总 KCL 通过同时出现；④优先两个最小复验及预期可证伪结果；⑤具体需要补导出的物理量/精度/坐标。
请勿凭“reference”标签默认 Sentaurus 正确，也勿凭 Vela 数值收敛默认其物理正确；不能据此直接改算法、调参、豁免 2% 或更换原始基线。
公式中必须区分单位、粒子流/常规电流、带符号面积和 split 高低位。不要把缺失的原生节点值推断为已有测量。

附件用法：此 Markdown 可单独上传；核心 ZIP 增补曲线、端口表、原生 PLT/CMD、配置、摘要与校验；详细 ZIP 再补六点 H5、全节点 CSV、网格/掺杂、残差/边诊断、日志及冻结源码子集。两包均含 SHA256_MANIFEST.json，包外有整体校验清单。
程序/原生 TDR 并非完整仿真运行包。只需分析时不应启动仿真。

---

'''
    (B/'LLM_HANDOFF.md').write_text(prompt+report,encoding='utf-8')
    (O/'report.md').write_text(report,encoding='utf-8')
    check=dict(cloud_files_verified=len(cloud),node_rows_verified=len(nodes),states_verified=6,
        exact_h5_csv_fields=True,input_hashes_verified=True,SRH_node_integrals_match=True,
        verified_frozen_source=copied,requested_source_not_available_or_not_hash_matched=unavailable,
        codespace_final_state='Shutdown',solver_invoked=False)
    (O/'package_validation.json').write_text(json.dumps(check,indent=2)+'\n')
    allfiles=[p for p in B.rglob('*') if p.is_file()]
    core=[]
    keep={'LLM_HANDOFF.md','all_points.csv','case_summary.csv','outliers.csv','neighborhoods.csv','high_vd_controls.csv',
          'vela_terminals.csv','field_statistics.csv','srh_source_integrals.csv','selected_attempts.csv'}
    for p in allfiles:
        rel=p.relative_to(B).as_posix()
        if rel in keep or rel.startswith('analysis/') and not rel.startswith(('analysis/verified_source/','analysis/analysis_scripts/')) or rel.startswith('frozen/') or rel.startswith('sweeps/') and p.name=='config.json':core.append(p)
    manifests={}
    for kind,files in [('core',core),('detailed',allfiles)]:
        dest=B.parent/f'simplemos_outliers_20260928_{kind}.zip'
        manifest={p.relative_to(B).as_posix():sha(p) for p in sorted(files)}
        with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in files:z.write(p,p.relative_to(B).as_posix())
            z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
        with zipfile.ZipFile(dest) as z:
            assert z.testzip() is None
            for p,s in manifest.items():assert hashlib.sha256(z.read(p)).hexdigest()==s
        manifests[dest.name]=dict(bytes=dest.stat().st_size,sha256=sha(dest),files=len(files))
    (B.parent/'simplemos_outliers_20260928_packages.json').write_text(json.dumps(manifests,indent=2)+'\n')
    print(json.dumps(dict(validation=check,packages=manifests),indent=2))

if __name__=='__main__':main()
