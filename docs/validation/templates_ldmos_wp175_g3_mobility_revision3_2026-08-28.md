# Templates/LDMOS WP1.75 G3 mobility 合同 revision 3（2026-08-28）

## 结论

WP1.75/G3 mobility 合同和生成逻辑已修正，并在 exact mesh 上完成无 predictor 的
Sentaurus 精确漏压路径复跑。原先阻断于 `0.031641657941 V` 的路径现已接受全部
15 个点并到达 `0.1 V`，说明上一轮失败的主要触发因素确实包含错误的 mobility
语义与约四个数量级的单位缩放。

但这只关闭了“路径可达性”问题，没有通过物理曲线校核。Vela 深截止电流仍由约
`274 A/um` 的 drift/diffusion 高位相消得到，净电流随偏压非物理振荡；在已有
Sentaurus 对照点上相差 `6--7 dex`。因此阶段 3 继续阻断，不启动 IALMob。

## 合同修正

Sentaurus `G3-no-IALMob` 明确使用：

```text
Mobility(HighFieldSaturation)
```

其中没有 `DopingDependence`。revision 2 却生成 `masetti_field`，还把合同中的
`1417 cm2/(V s)`、`1.07e7 cm/s` 等 TCAD 内部数值再次转换为 SI 后交给
`scaling.mode=unit_scaling`，造成 mobility 约低 `1e4`。

revision 3 作以下单因素修正：

- 新增生产模型 `constant_field`：低场 mobility 直接取材料合同，随后只应用现有
  HighFieldSaturation limiter，不读取掺杂 mobility 参数；
- G3 生成 `constant_field`，G4/零偏置 classical 分支生成 `constant`；
- G3 deck 只写 high-field 驱动力、离散方式、饱和速度和 beta；
- legacy `unit_scaling` 的这些键直接接收 TCAD 内部值 `1.07e7/8.37e6 cm/s`，不再
  乘 `1e-2`；
- 合同显式记录 `doping_dependence_enabled=false`、
  `high_field_saturation_enabled=true`，revision 从 2 升至 3；
- Masetti 查询参数仍保留在合同中，但在 deck 明确启用 `DopingDependence` 前不参与
  G3 计算。

本次没有启用 predictor、没有修改偏压点、Newton、line search、接触、BGN、复合或
离散 profile。

## Exact-mesh 无 predictor 结果

运行链为零偏置 Poisson bootstrap、coupled equilibrium、Save/Load reclose，随后按
Sentaurus 原始路径运行 15 个漏压点。三个前置步骤均通过，主路径结果如下：

| 指标 | revision 2 | revision 3 |
| --- | ---: | ---: |
| 接受点数 | 9，下一行记录失败点 | 15/15 |
| 最后接受偏压 | `0.023155977422 V` | `0.1 V` |
| `0.031641657941 V` | 100 次后 `max_iterations` | 37 次，`reltol` |
| 最大迭代数 | 100（失败） | 76（`0.097722431245 V`） |
| 主路径墙钟 | 未作为本轮重测基准 | `782.53 s` |

路径可达性硬门通过，但电流没有通过：

| Vd | Sentaurus drain TotalCurrent | Vela revision 3 | 幅值误差 |
| ---: | ---: | ---: | ---: |
| `0.023155977422 V` | `2.0832e-15 A` | `4.5484e-9 A/um` | `6.34 dex` |
| `0.031641657941 V` | `2.3352e-15 A` | `7.1084e-8 A/um` | `7.48 dex` |
| `0.1 V` | 本次 capture 未覆盖 | `3.2205e-12 A/um` | 不评分 |

Sentaurus 原生电流以二维默认深度输出为 A；当前表格沿用既有 `1 um` 对比约定，端口
定义仍为诊断等级。即便不对前两点作硬评分，Vela 曲线自身在相邻偏压间跨越多个数量级
且改变符号，已经足以判定状态/通量闭合不合格。

## 数值观察

revision 3 在 `0.1 V` 的平均电子 mobility 为 `1409.69`、最小值为 `3.77`
（TCAD 内部 `cm2/(V s)`），证明 constant-field 与 high-field limiter 已生效。
同点电子 drift/diffusion 分量分别为 `+274.273966858849` 和
`-274.273966858846 A/um`，净电子电流仅 `3.2205e-12 A/um`。约 14 位有效数字的
相消使 residual 达标并不等同于端口电流可信。

上一轮固定 Sentaurus 状态回放已经显示，正确 HFS-only mobility 能把主载流边 SG
中位误差从 `3.884 dex` 降至 `0.119 dex`。本轮自洽曲线仍失败，因而剩余主问题位于
状态闭合/Poisson 合同与高位相消精度，而不是继续修改 mobility 公式。

## 验证与后续门

已通过：

- `test_mobility [mobility]`：26 个用例、129 项断言；
- 新增 `constant_field` 掺杂不敏感与 HFS 公式单测；
- Templates/LDMOS phase 0/1 与 phase 2/3 共 36 项 Python 回归；
- revision 3 合同通过 JSON schema 严格校验；
- exact-mesh 无 predictor 全路径 15/15 收敛。

下一步保持单因素纪律：

1. 在 `Vd=0.023155977422 V` 重新保存 revision 3 checkpoint，并对比 Sentaurus/Vela
   `psi/phin/phip/n/p`；
2. 对固定 Sentaurus 状态拆解 Vela Poisson 的介电、净掺杂、载流子、BGN/Fermi、
   接触和介质边界贡献；
3. 单独审计仍沿用 legacy `unit_scaling` 的非 mobility 字段，尤其 BGN/SRH 浓度和
   Auger 系数，禁止与 Poisson 修复合并；
4. Poisson 和端口高位相消关闭前，不启用 secant predictor，不进入 IALMob。

## 产物

- 版本化合同：`reference_tcad/templates_ldmos_sentaurus2022/contracts/physics_contract.json`
- 生成器：`scripts/prepare_templates_ldmos_phase23.py`
- mobility 实现：`src/physics/MobilityModel.cpp`
- ignored 运行目录：
  `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`
  `stage1_v4/phase23_t2022_contract_v3/`
