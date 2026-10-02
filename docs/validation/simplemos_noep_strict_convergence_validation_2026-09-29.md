# SimpleMOS：仅删除 ExtendedPrecision(128) 的四次判别实验

2026-09-29。按用户授权在 Sentaurus 虚拟机执行。生产代码、参考数据、2% 电流标准及冻结数值门槛未修改。

## 结论

四次运行全部正常完成，各输出 51 个 Vg 点。**双精度保留 Digits=15、RhsMin=1e-15 可以跑完，但没有复现 EP 对照的精度。**
默认端口算法的 102 个偏置点中，Vela 相对本次原生参考的最大绝对 Id 误差为 **0.196732302%**，在本轮范围内通过 2% 电流比较标准；这不等于全部数值资格通过。
Direct 在深关断区仍有大幅偏差。全部 248 次非线性求解通过更新误差条件停止，最终 RHS 没有一次小于 1e-15。

因此，不能把提高 Digits、降低 RhsMin 当作提高浮点精度的替代。本轮也未显示比 M60 默认端口参考更稳定、更准确的普遍改善。

## 实验约束与执行核验

- n23/n24 × default/Direct，共四次；Vd=0.05 V，Vg=0–2.5 V、间隔 0.05 V。沿用 EP 完整 deck 的扫描范围，没有扩展器件或 Math 组合。
- 从 [前轮 EP 实验](simplemos_ep_and_three_state_validation_2026-09-28.md) 的冻结输入生成，仅删除 `ExtendedPrecision(128)`；输出名称单独改名。脚本逆向还原名称后逐字符验证无其它差别。
- 保持 `Extrapolate RelErrControl Digits=15 RhsMin=1e-15 ErrRef(Electron)=1e2 ErrRef(Hole)=1e2 Iterations=20 ExitOnFailure CNormPrint`；Direct 两份仍含 `DirectCurrent`。未增加 AND 收敛或其它参数。
- 保持同一 TDR、接触、300 K、Boltzmann、OldSlotboom、PhuMob、Enormal、HighFieldSaturation、SRH(DopingDependence)，以及初始化和偏置推进流程。电流按原生 A/µm，电子+空穴=总电流。
- TDR SHA256：n23 `d4c50bf96ccaf817389255071fa07271de20eb95fa8ef5be74ecfd204bbf5016`；n24 `336f0dac984a1821674c7fb558e59d59c5f347ba1f487cd7d513527a142ad128`。均与 EP manifest 相同。
- Sentaurus T-2022.03-SP2，主机 tcad。四份日志均明确回显 `Use 64 bit (double) normal precision floating point arithmetic.`、15 digits、RHS 下限 1e-15；求解后端回显为 `blocked decomposition`，未显式改变 Method。不用其它模型的 Super 回显推断本 DD 内层后端。
- 08:53:20 至 09:00:34（UTC+8）。两个器件并行，同器件 default、Direct 顺序运行；四次退出码均为 0，无重试、无参数修补。
- 8 个上传输入哈希通过；下载包 SHA256 和 68 个唯一输入/结果文件哈希通过。4×51=204 条电流记录、102 个 default/Direct 配对偏置点；每次五份状态 TDR 随包保留，本轮未新增逐节点场对比。

## 默认算法的电流比较

误差统一定义为 `(Id_Vela / Id_对应原生默认法 − 1) × 100%`，Vela 使用既有冻结结果，没有重新求解。

| 器件 | Vg/V | M60 双精度，Digits=8 / % | 本次双精度，Digits=15 / % | EP128，Digits=15 / % |
|---|---:|---:|---:|---:|
| n23 | 0 | −0.027081451 | +0.196732302 | +0.001257559 |
| n23 | 0.05 | −0.186617395 | +0.181668766 | +0.002273396 |
| n23 | 0.10 | +0.018573936 | +0.004300982 | +0.003008902 |
| n24 | 0 | −0.072726302 | −0.054143826 | +0.001463572 |
| n24 | 0.05 | +0.126094009 | +0.122856436 | +0.002441642 |
| n24 | 0.10 | −0.012937167 | −0.002833557 | +0.003133248 |

本次 n23 最大误差在 Vg=0，为 0.196732302%；n24 最大误差在 Vg=0.05，为 0.122856436%。
n23 在 Vg=0.05 相对 M60 的误差变号，绝对量级仍约 0.18%，不能解释为 Digits 增大后单调逼近真值。

## 观测算法与守恒：Vg=0.05 V

| 量 | n23 | n24 |
|---|---:|---:|
| 本次 default Id / A/µm | 1.04045784792144e-16 | 1.12989659540853e-16 |
| 本次 Direct Id / A/µm | 1.42227355588159e-15 | −7.46547509923595e-15 |
| `(Direct/default−1)×100%` | +1266.968935% | −6707.219749% |
| 衬底电子 `(default−Direct)/abs(Id_default)` | +3.117521929e-3 | +9.330086839e-4 |
| default `ΣI/abs(Id_default)` | −6.044819367e-4 | +3.599943622e-15 |
| Direct `ΣI/abs(Id_Direct)` | +1.987462141 | −2.002302654 |

本次默认法全曲线最大 `abs(ΣI)/abs(Id)` 为 6.044819367e-4（0.060448194%）。n24 在该点的 KCL 接近零，但电流差仍为 0.122856436%，再次说明 KCL 不是电流误差条。
前轮 EP 对照中 default/Direct 与 KCL 均已接近输出舍入量级；本次删除 EP 后，这些观测差没有消失。
上述 default/Direct 是独立运行的配对输出，不宣称原生内部状态位级相同，也不据此唯一分离端口误差和状态误差。

## 为什么 RHS 下限为 1e-15 仍然结束

四份日志共 248 个非线性求解表，包括各自的 Poisson 初始化：全部停止原因为 `Error smaller than 1`，没有一次最终 RHS 小于 1e-15。
排除每次最初的 Poisson 求解后，244 个耦合 DD 表的最终 RHS 范围为 **0.0102–0.0366**，这是求解器内部尺度，不是安培或相对电流误差。

Vd=Vg=0.05 V 的默认算法日志中：

| 器件 | Newton 次数 | 第 2/3/4 步 RHS | 最后一步 step | 最后更新 error |
|---|---:|---|---:|---:|
| n23 | 4 | 0.0153 / 0.0153 / 0.0153 | 9.96e-16 | 0.34173 |
| n24 | 4 | 0.0233 / 0.0233 / 0.0233 | 9.70e-16 | 0.31252 |

RHS 停在平台，更新继续减小，随后满足 error<1；降低 RhsMin 没有使它变成必须通过的条件。当前逻辑允许更新误差条件接受，本次未启用 `RhsAndUpdateConvergence`。
该观察与双精度求值限制一致；结合唯一删除 EP 的对照，支持算术精度对剩余差异有实际影响。它不证明所有历史误差仅由某个端口公式的舍入引起，也不把 Digits 解释为最终解保证获得的有效数字。

## 运算代价（描述性记录）

| 运行 | 本次 double CPU/s | 前轮 EP CPU/s | 本次 wall/s | 前轮 EP wall/s | 本次/EP Newton 步数 |
|---|---:|---:|---:|---:|---:|
| n23 default | 24.46 | 124.48 | 214.96 | 314.66 | 267 / 256 |
| n23 Direct | 22.75 | 115.79 | 213.54 | 306.18 | 267 / 256 |
| n24 default | 24.82 | 119.19 | 215.51 | 310.21 | 267 / 256 |
| n24 Direct | 22.79 | 140.94 | 213.61 | 331.50 | 267 / 256 |

本次 CPU 用时明显较少。但前轮顺序运行，本轮两个器件并行，扫描的自适应内部求解表数也不同（本次 62、EP 65）；wall 含初始化等等待。本表不是受控性能基准，不用于声称通用 EP 减速倍数。

## 证据和复现

- 准备脚本：[prepare_simplemos_noep_20260929.py](../../scripts/prepare_simplemos_noep_20260929.py)。专用输出目录必须不存在，防止覆盖证据。
- 分析脚本：[analyze_simplemos_noep_20260929.py](../../scripts/analyze_simplemos_noep_20260929.py)。支持失败/部分曲线记录；精确偏置匹配、重复点检查、保留全部电流分量。
- 本地证据目录：`build/simplemos_noep_20260929/`；远端：`sentaurus:/tmp/vela_simplemos_noep_20260929`。
- 输入与来源：`manifest.json`、`deck_changes.diff`、`inputs.sha256`；结果：`summary.json`、`run_checks.csv`、`comparison.csv`、`observer_pairs.csv`、`terminal_components.csv`、`nonlinear_blocks.csv`、`newton_iterations.csv`；原始日志和状态位于 `raw/`。
- 完整包 SHA256：`3ffa3d4e79c267a8f3ceed0ae5679065553424304a17752d5b693b8a9d5e4e25`。本地逐文件核验记录为 `verified_download_hashes.json`。
- 分析器先用既有 n23 EP 控制验证 51 点解析、末步 RHS/error、wallclock，随后处理本次四条完整曲线；原始日志另行核查关键停止行和 Newton 行。

本轮没有重新执行 Vela 双初始化，也没有扩展 EP 到其它器件或 Vd。结论仅覆盖上述两个器件、同一模型和网格、Vd=0.05 V 的扫描；不更新已有全工况参考资格。
