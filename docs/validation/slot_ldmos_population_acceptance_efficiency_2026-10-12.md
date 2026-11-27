# SLOT 接受口径与诊断开销修订（2026-10-12）

## 范围和结论

本批只执行用户选定的两项：统一接受口径、减少重复装配及大文件输出。
旧 `build/s12home60` 证据保持不变；未恢复 0→60 V 扫描或计时器，未连接远程。
实现是针对该冻结 Python 驱动的显式选用包，不改变全局 C++ Newton 默认验收。

在相同网格、相同精确偏压和第16次保存态上，重放 d003/d007 的一次更新。
两次场量均与旧第17次候选态逐位一致；原来因零更新 RHS 复读失败而继续迭代的
候选态，现在按原始更新误差分支接受。这里验证的是这两份冻结候选态的处理，
不是从零开始重跑整条路径，也不是不同初值独立复求。

## 接受规则

`scripts/slot_population_acceptance.py` 统一使用：

```text
(rho <= 1e-10 OR raw_physical_update_error < 1)
AND closure AND boundary AND linear_quality AND inverse_mapping AND positivity
```

更新误差取原始线性方向在物理人口坐标中的归一化误差，不取投影/阻尼后的步长。
载流子行比继续报告；闭合0.01、边界、人口正性、逆映射1e-7和线性混合误差保护保留。
合格出口仍做边界物理残差检查、保存态重读及全部接触分载流子电流诊断。

零更新 C++ 求解器无法重建前一步原始更新误差，所以其 RHS 是否达标单独记录，
不再二次否决已满足更新分支的候选态。额外保留保存态一致性检查：
`abs(rho_readback-rho) <= max(1e-15,1e-6*abs(rho))`，使用相同内部残差单位。
这是重读一致性保护，不是新的器件收敛阈值；两次验证中的残差完全相同。
零更新探针的退出码1因此是预期诊断结果，不表示本驱动的更新分支失败。

|冻结候选|物理Vg/Vd (V)|整体残差|原始更新误差|复读RHS|统一资格|
|---|---|---:|---:|---|---|
|d003第17次|5 / 0.016856999999999997|1.7011140352543726e-10|0.5317051528707994|未通过；重读一致|通过|
|d007第17次|5 / 0.014022620833333332|2.55809161842238e-10|0.9221515641846297|未通过；重读一致|通过|

准确偏压以 `build/s12home60/*_17_matrix/config.json` 为准。
首次 `build/s12efficiency` 用手写舍入偏压，发现末位偏压不同，不能用作逐位一致性对照；
该记录保留，正式证据使用 `build/s12efficiency_exact`。
验证编排曾在完整闭合补充探针启动前漏设后端环境，已用独立标签补齐；没有覆盖求解记录。

## 开销修订

1. 每步闭合探针关闭 `candidate_port_diagnostics`：继续计算真实全局连续性闭合，
   不再在同一探针重复做四个接触泛函和整体残差装配。整体残差取已有逐行探针
   `status.json` 的 `block_residuals.combined`；合格出口仍计算全部接触电流。
2. 不再为了取得整体残差而重复读取并归约整份 rows CSV；逐行诊断本身仍保留。
3. `newton_linear_system_probe` 新增 `linear_system_output_profile=solve_only`，
   保留实际 `stabilized.mtx` 和含 residual/RHS/direction 的 `vectors.csv`，
   省去每步 `physical.mtx`、`regularized.mtx`。默认 `full` 不变，非法值报错。
   装配和求解过程完全相同，未削减线性后向误差所需证据。

本批尚未把逐行输运与闭合两项必要装配融合成一次 C++ 内存调用；不能称为已消除
全部重复装配，也不能把单项提速当作整个扫描的提速倍数。

完整/精简闭合探针逐字段比较相同。d003：29.7175→4.8467秒；
d007：29.4136→5.5061秒。这是同态各一次测量，主要反映被省去的重复诊断。
同一d007冻结系统的矩阵输出从71,084,544字节降至26,307,584字节（减少63.0%）；
实际求解矩阵及向量文件SHA-256相同，向量文件同时包含残差、RHS和原始方向。
全量矩阵探针22.7447秒，精简20.5052秒；此处节省主要为输出空间，不能推导
线性求解器本身变快。哈希和完整计时见下列机器证据。

## 实现与复现

新脚本 `scripts/prepare_slot_population_efficiency.py` 从审计过的 s12home60 快照创建
全新目录，绑定输入哈希，遇到快照结构变化拒绝生成；只准备文件，不启动求解。
生成的 `continuation.py` 已移除两处父态链的零更新 RHS 二次门槛。
后续使用仍需独立的用户运行窗口、监督和状态监测配置。

```powershell
D:/msys64/ucrt64/bin/python.exe -B scripts/prepare_slot_population_efficiency.py --source build/s12home60 --destination build/NEW_ISOLATED_PACKAGE --runner build-release/vela_example_runner.exe
```

固定物理：SLOT 原网格、Vg=5V、AverageBox/Fermi/有效带边SG/原电离模型；
PARDISO严格消元，单线程，beta=0.001浓度有限步，电势对数阻尼，载流子正则化和
闭合方程替换关闭。物理电流单位A/µm，电压V，矩阵/残差采用原内部缩放。

关键证据：

- `build/s12efficiency_exact/validation.json`：两候选资格、场量逐位一致、闭合一致、矩阵文件哈希。
- `build/s12efficiency_exact/new_d003_records.json`、`new_d007_records.json`：全部保护与接触电流。
- `build/s12efficiency_exact/runs.json`：每个探针的退出码和实测耗时。
- `build/s12efficiency_exact/final_bindings.json`：实现、冻结runner、配置及父态的哈希绑定。
- `build/s12efficiency_profile_rejection/result.json`：非法输出配置在装配前拒绝。
- `build/s12efficiency_resume.py`：有界单更新及全量/精简诊断验证脚本。

## 验证与限制

- UCRT64 Release `vela_example_runner` 构建通过。
- `slot_population_acceptance`、`slot_linear_quality` CTest两项通过，共13个Python用例。
- `test_slot_candidate_conservation` 4个用例通过；初次默认临时目录受ACL限制，
  改用工作树下临时目录后通过，未改测试断言。
- 非法profile拒绝验证通过；全量/精简数值一致性见机器证据。
- 既有h5py/HDF5版本警告仍存在，未因本任务安装或升级依赖。

该修订解决已选标准被重复否决与诊断开销问题，不证明所有减半都能消除。
此前21.614mV等目标的实际非线性振荡仍需另行处理；未宣称60V路径完成，
未进行Sentaurus定量比较或全场独立复现。
