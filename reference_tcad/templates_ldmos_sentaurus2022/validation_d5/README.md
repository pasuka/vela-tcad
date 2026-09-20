# 可直接运行的 LDMOS D5 验证算例

本目录把已经验证的中性数据补齐到版本管理，不再依赖本机的
`reference_staging/`、Sentaurus 安装或虚拟机。包含精确网格、掺杂、
AverageBox 传输权重、材料、双栅压合格零漏压种子和原生参考曲线。
原始商业工具二进制输出和参数文件未纳入；来源及逐文件 SHA-256 见
[provenance.json](provenance.json)，运行依赖合同见 [inputs.json](inputs.json)。

## 算例与默认值

- 10,241 节点、19,782 三角形；坐标 m，掺杂输入 cm⁻³，状态密度 m⁻³，
  电势 V、电流 A/μm。参考文件保留原始有符号漏极电流及零压记录，
  精确参考点由既有解析器提取。
- D5 等温经典漂移扩散，300 K，Vg=4/8 V，各 31 个精确漏压点，0–40 V。
  原 Fermi/BGN、掺杂与高场迁移率、SRH/Auger（关闭 Auger 生成）和离散口径
  原样保留；D5 不启用 IALMob、量子、自热或雪崩，不能代替 D4/D0 验收。
- 默认 UMFPACK、主 Newton 跨点符号分析复用、单线程环境。
  可按 UMFPACK、SparseLU、STRUMPACK、MUMPS、SuperLU_MT 的优先次序显式选择。
  指定后端必须在实际构建中可用，失败不会静默切换。
- 从合格的 Vg=4/8 V、Vd=0 V 种子开始。此入口不重跑栅压预偏置，
  也不把初始化时间纳入曲线时间。

## 运行

使用 MSYS2 UCRT64 Release，Python 3；既有 linked 驱动使用 Windows 进程计时 API。
从仓库根目录运行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target vela_example_runner --parallel 2

# 快速输入检查，不启动仿真、不创建输出。
python scripts/run_templates_ldmos_reference.py --preflight --output build/reference_tcad/ldmos_preflight

# 双栅压各前 8 点（默认）；输出目录必须是新的。
python scripts/run_templates_ldmos_reference.py --output build/reference_tcad/ldmos_d5_first8

# 完整 62 点并执行原双栅压联合验收。
python scripts/run_templates_ldmos_reference.py --points 31 --output build/reference_tcad/ldmos_d5_full

# 显式选择另一已编译后端。
python scripts/run_templates_ldmos_reference.py --linear-solver sparselu --output build/reference_tcad/ldmos_d5_sparselu
```

入口沿用 [linked 驱动](../../../scripts/run_templates_ldmos_linked_d5.py)的电压推进、
失败恢复、逐行闭合及参考系门限；不做新的参数优化。
前八点只取得该区间数值资格，完整曲线才检查全部电流及双栅压联合指标。
网格、种子与材料哈希在执行前检查；每个已接受状态、原块门限、逐行门限和
KCL 在结束后独立审计。完整曲线的工程/最终判据来自
[既有评分器](../../../scripts/analyze_templates_ldmos_stage4_d5.py)，未修改。

输入回归：`python -m unittest tests.regression.test_ldmos_reference_fixture`。
输出与性能文件保留在忽略的 `build/` 内。历史完整曲线资格与本次默认值验证
分别记录于 [当前状态](../../../docs/validation/templates_ldmos_current_status.md)。
