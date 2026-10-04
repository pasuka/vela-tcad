# SimpleMOS 统一状态契约生产迁移与十六点复验

日期：2026-09-12。当前完成 **32/32 状态、16/16 双初始化**。32 次首次求解全部通过，无重载。16 个合格状态的 336 项正式方向 Jv 检查全部通过，最大相对差 6.543020944e-7（原门槛 1e-4）。

## 实现与范围

新增显式 `solver.split_dd_state=true`，默认关闭。Newton 线搜索累积 ψ、电子与空穴准费米增量的低位；Poisson 电荷、稳定 SG、实时 PhuMob、SRH、端口与验收边/行探针共同读取同一状态。状态相关残差在 100 位十进制精度下求值，公共输出投影为 double。方向残差差分在投影前相减，避免再次引入残差消减。

CSV 保留三块低位、完整参考场、归一化尺度、模式版本和网格指纹。陈旧或不完整的低位状态被拒绝；不再把低位当作可丢弃的提示。无载流子约束行继续采用原编码参考定义。只读端口残差积分显式接收完整状态；未适配的弧长更新继续拒绝。

当前限定为 300 K、Boltzmann、OldSlotboom、plain PhuMob + element_box_phumob 和掺杂 SRH。保留现代 SI 常数、Si 基础 ni=1.4638805412559193e10 cm⁻³、原 SRH 体积、原接触条件、四次线性迭代修正、更新限幅及所有接受门槛。原解析 double Jacobian 保留，不能把残差的 100 位求值误称为整个 Jacobian 也使用相同精度。

几何仍由独立配置明确选择：Delaunay box 转移、cell_material 介电装配、element_box 硅侧输运及 signed_transport 三项 Poisson 电荷体积。本选项不替代这些几何配置。n19/n23 分别为 1480/1482 节点、2742/2746 个三角形；NWell 为 1e17/2e17 cm⁻³。网格坐标 μm，状态文件浓度 m⁻³，电流 A/μm，宽度 1 μm。

## 冻结的比较与结果

同原来的 32 份独立种子，覆盖 n19/n23 × Vd=0.05/1 V × Vg=0/0.2/0.8/1 V。每个失败目标最多同偏置重载一次，旧失败记录保留。本轮未使用额外重载。

| 检查 | 最坏值 | 原门槛 |
| --- | ---: | ---: |
| 全部 1814 个自由载流子行 | 7.354958334e-7 | 1e-6 |
| KCL/Id | 4.222310596e-14 | 1e-8 |
| 端口残差积分/提取差 | 4.440892099e-16 | 1e-8 |
| 双初始化 Id 相对差 | 4.329869796e-14 | 1e-6 |
| 双初始化 ψ 最大差 | 3.22e-15 V | 1e-6 V |
| 双初始化 φn 最大差 | 9e-15 V | 1e-6 V |
| 双初始化 φp 最大差 | 2.938020871e-8 V | 1e-6 V |
| 双初始化密度相对差 | 1.136477783e-6 | 1e-4 |

最差逐行比值位于 n23、Vd=1 V、Vg=0.8 V 的 Vela 初始化状态，仍较接近原门槛，未称为大裕量。最多 5 次 Newton 更新。先前六个原生初始化失败点均在首次尝试通过，原有 26/32、10/16 的历史失败不删除。

全局连续性检查仍保留原 1e-10 源下限；64 个载流子源积分分量均低于此下限，因此不宣称它们已通过任意相对源精度或任意局部源扰动校准。

## n23 电流对照

误差为 `100 × (Id_Vela/Id_Sentaurus − 1)`。这里只列本轮合格的四个栅压，不当作完整曲线，也不新增电流百分比验收线。

| Vg/V | Vd=0.05 V | Vd=1 V |
| --- | ---: | ---: |
| 0 | -0.018196719% | -0.808135237% |
| 0.2 | +0.003584516% | -0.016321308% |
| 0.8 | +0.005693624% | +0.006529258% |
| 1 | +0.005077499% | +0.005109832% |

高 Vd 深关断点：Vela 为 3.68309865058048e-16 A/μm，Sentaurus 为 3.71310556502558e-16 A/μm。数值迁移解决了统一实现下的收敛资格缺口，没有据此声称消除了这项物理电流差。

## 前置检查、回归与失败保留

- 最终两个真实状态均逐字节重启一致；边通量对节点通量重构相对项尺度误差最大 2.142884917e-16，端口差最大 2.220446049e-16。
- 两态、三块坐标、1e-6/1e-20 V 局部方向的非弱交叉块 Jv 最大相对差 6.340501512e-10，通过 1e-4 门槛。这不代替任意弱交叉源块的完整性证明。
- 两个默认关闭控制：与本轮修改前的程序相比，状态 CSV 字节和数值状态字段均完全一致。
- UCRT64 Release 构建成功。七个相关 Catch2 程序全部通过：267 个用例、6875 个断言，覆盖 Newton、DC sweep、CSV、SG、box、线性修正及新状态契约。本轮未把它们描述为全库 CTest。
- 保留最初 CSV 白名单遗漏导致的测试失败、首个只读端口探针被弧长保护拦截的失败，以及均匀载流子测试发现的新电子漂移诊断符号问题。最终漂移/扩散分解满足均匀密度无扩散的性质；总电流公式未变。`LongDoubleReference` 命名的旧诊断字段在 split 模式下来自同一宽精度算子，不是独立精度证据。
- 最终算例使用 Eigen SparseLU/COLAMD；构建也检测到 HDF5、SPQR 和 UMFPACK，但没有把后两者作为这批 DC 的实际求解后端。

## 顺序与证据

生产迁移与 16/16 双初始化阶段已通过，包括 16 态 336 项正式 Jv 检查。PhuMob 四条 0–1 V、步长 0.02 V 的曲线输入已冻结，Vela 延续路径已启动。Enormal → 高场饱和 → 原始 n17–n24 × 两个 Vd、0–2.5 V 范围保持后续阶段；当前尚未运行这些恢复模型。Sentaurus 运行版本为 T-2022.03-SP2；本轮曲线上传/执行与结果另行记录，未取得完整曲线结果之前不放行模型恢复。

- [生产验证合同](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/validation_contract.json)
- [全部 32 次尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/attempts.csv)
- [16 点比较与双初始化](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/comparison.csv)
- [比较证据](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/comparison_evidence.json)
- [两态预检](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/preflight.json)
- [16 态 Jv 结果](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/jvp_summary.json)
- [默认关闭一致性](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v3/default_identity_evidence.json)
- [完整曲线合同](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/vela_contract.json)
- [验证执行器](../../scripts/validate_simplemos_split_production_20260912.py)
- [相关回归记录](../../build-release/split_prod_0912/regressions/summary.json)
- [第一版探针失败及归档](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/candidate_v1_archive.json)
- [端口诊断分量修正前的单点证据](../../reference_tcad/simplemos_sentaurus2022/phumob_split_production_20260912/v2/pilot_evidence.json)

本轮尚未提交或推送；已有无关本地修改保留。
