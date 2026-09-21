# Templates/LDMOS G3 continuity、SG 与漏极接触分支审计（2026-08-28）

## 结论

`Vd=0.0231559774221138 V` 的六个数量级漏电差已经定位到 **Vela 相对残差提前
收敛留下的漏极第一圈电子准费米势尾差**，不是 SG 公式、端口积分、浮点相消或
IALMob 缺失。

原 revision-4 扫描状态虽然以 `reltol` 通过，但电子 continuity L2 仍为
`1.055080877e-4`；其中节点 702（`x=-9.976355 um, y=10.076172 um`）一项为
`1.055080874e-4`，几乎占满整个范数。该节点与漏极接触节点 5546 相连的 edge 2268
有 `8.21284e-12 V` 电子 QF drop，并贡献 `4.57735e-9 A/um`，主导原端口电流
`4.65890e-9 A/um`。

保持全部 revision-4 物理、网格、离散和偏压不变，只把同点 reclose 的
`reltol/abstol/stall_residual_floor` 收紧到 `1e-10/1e-14/1e-12` 后，最后有限状态将：

- 电子 continuity L2 降到 `1.608768487e-9`，改善 `65583` 倍；
- 漏极第一圈电子 residual L2 降到 `7.1764e-12`；
- 漏极电流降到 `2.195335497e-15 A/um`；
- 相对 Sentaurus PLT 的 `2.083166960e-15 A/um` 仅差 `5.3845%`，即
  `0.02278 dex`。

该 strict reclose 在第 35 次迭代因 `line_search_non_decrease` 停止，未满足形式上的
`1e-10` 相对门槛，因此它是根因控制而不是可提交的正式曲线点。阶段 3 仍保持阻断，
下一开发项应是 block-aware 绝对收敛/重闭合和 line-search 数值地板处理；暂不启动
IALMob，也不启用 predictor。

## 审计合同

- 网格：Templates/LDMOS exact imported topology，10241 节点、30022 条 SG 边。
- 物理：revision 4 的 Fermi statistics、OldSlotboom BGN、SRH/Auger 和
  `constant_field` HFS-only mobility，审计期间不改变。
- 偏压：`Vg=0 V, Vd=0.0231559774221138 V`。
- Vela 状态：revision-4 无 predictor 精确偏压序列最后接受 checkpoint。
- Sentaurus 状态：同偏压 TDR 导出的 `psi/QF/n/p`，只作固定状态算子回放。
- Sentaurus 电流真值：`g3_capture_min_transition_states.plt` 中同偏压、最高 time 的
  drain `TotalCurrent`，不从 TDR 的节点 QF 差反推。
- SG 守恒：生产 `sg_edge_flux_probe`；节点散度按 node0 加、node1 减重建。
- continuity：生产 `newton_carrier_term_probe`。
- 端口：生产 `ContactCurrent` contact cut，默认 `1 um` 深度。
- strict control：无 predictor、无物理变化、从原 checkpoint 同偏压 reclose。

## 三条生产恒等式

| 状态 | SG cut vs ContactCurrent | ContactCurrent vs 曲线 | SG 节点散度 vs carrier flux 最大绝对差 | stable vs long-double |
| --- | ---: | ---: | ---: | ---: |
| Vela 原接受状态 | `8.34e-15` | `0` | `9.02e-23` | `8.94e-15` |
| strict reclose 最后状态 | `4.67e-15` | `0` | `2.48e-24` | `4.33e-15` |
| Sentaurus 导出状态 | `5.32e-4` | `2.52e-16` | `1.65e-24` | `5.43e-15` |

Vela 自洽状态的 SG 边算子、节点散度、端口 cut 与曲线电流全部闭合。Sentaurus 导出
状态上 `5.32e-4` 的 SG/contact 差不是 SG 分支错误：Newton SG 算子从 `psi/QF`
重构载流子密度，`ContactCurrent` 在 frozen replay 中使用随 TDR 一起提供的原始
`n/p`。漏极边两种密度的中位/最大差分别为 `2.796e-4/4.934e-4 dex`，足以解释该
诊断差。Vela 自洽状态的对应差为零。

## 漂移/扩散高位相消

原 Vela 状态的漏极电子端口拆分为：

| 项 | 电流（A/um） |
| --- | ---: |
| stable SG | `4.658902479e-9` |
| long-double SG | `4.658902479e-9` |
| drift | `+274.537715651` |
| diffusion | `-274.537715646` |

终端拆分条件数为 `1.17855e11`。strict reclose 最后状态的条件数进一步达到
`2.50099e17`。单边上更明显：Sentaurus 导出状态的 edge 15948 stable SG 为
`-2.6075e-13 A/um`，而 drift+diffusion 的 double 拆分为约
`-3.21e-2 A/um`。因此 drift/diffusion 只能用于展示物理大项，不能用于低漏电评分；
生产 stable SG 与 long-double 的闭合才是数值判据。

## 漏极接触第一圈定位

原接受状态的第一圈电子 residual L2 为 `1.055080875e-4`，第二圈仅
`3.0118e-9`。主要记录为：

| 节点/边 | 位置或连接 | 指标 |
| --- | --- | ---: |
| node 702 | `(-9.976355, 10.076172) um`，drain ring 1 | electron residual `1.055080874e-4` |
| edge 2268 | contact 5546 -> interior 702 | stable electron `4.577354057e-9 A/um` |
| edge 2268 | 同上 | electron QF drop `8.212844e-12 V` |
| edge 15922 | contact 5547 -> interior 701 | stable electron `3.043646128e-11 A/um` |

strict reclose 最后状态中，漏极主导边的 QF drop 已降到约 `0.7--1.4e-17 V`，端口
电流与 Sentaurus PLT 接近。此时最大 electron residual 转移到
`(-9.967628, 2.578125) um`，不再位于漏极第一圈，证明原异常属于接触邻域的未闭合
状态分支。

## Sentaurus TDR 固定状态电流的使用限制

Sentaurus 导出状态回放得到约 `-2.00047e-13 A/um`，比 PLT 真值大约两个数量级。
漏极边上的绝对 QF 值约为 `0.023155977422... V`，真正决定 fA 电流的边差已落在
`1e-17 V` 量级；普通绝对 double CSV 无法稳定保存这种差。Vela 的自洽状态使用
`qf_reference + qf_increment` 坐标保存小差，因此 strict reclose 可得到 fA 电流。

据此，后续低漏电真值顺序必须是：Sentaurus PLT 端口电流 > Sentaurus 内部残差/
NewtonPlot > TDR 绝对 QF 固定状态回放。TDR 回放仍可用于空间趋势、密度和大于其精度
地板的 SG 边，不得作为 fA 端口硬门。

## strict reclose 的剩余阻断

strict control 从初始 combined residual `1.055081067e-4` 降到
`4.352647783e-8`，其中 `psi/phin/phip` 分别为
`4.349600515e-8 / 1.608768487e-9 / 2.523340795e-10`。第 35 次迭代的 Newton step
norm 为 `178.75`，13 次 line search 均未降低 merit，最后报告
`line_search_non_decrease`；接触多数载流子 QF drop 诊断约 `2.1275 V`。

这说明简单地把相对门槛从 `1e-7` 收到 `1e-10` 不是生产修复。原门槛会在大初始
残差路径上过早通过，严格门槛又会越过可用的低电流状态并撞上 Poisson/line-search
地板。需要显式区分：

1. continuity block 的绝对/端口分辨率门；
2. Poisson 数值地板；
3. 最佳有限迭代状态保存与正式接受规则；
4. contact-majority QF branch guard。

## 决策和后续顺序

1. **SG 与端口积分关闭**：不修改 SG 公式，不修改 contact cut，不用漂移/扩散拆分
   反推物理参数。
2. **当前根因归 WP1.5 求解器资格**：新增 block-aware absolute convergence/reclose
   设计，使大初始残差不能掩盖 leakage-dominant continuity 尾差。
3. 对 strict reclose 的迭代 0--35 做 residual/current 同步轨迹，确定电流已稳定而
   Poisson 到达数值地板的安全接受窗口；不得直接接受任意 rejected state。
4. 增加 exact-mesh 回归：同点 reclose 后 drain current 相对 Sentaurus PLT 误差
   `<=0.15 dex`，且漏极第一圈 electron residual 不再主导全局范数。
5. 上述资格门通过后再恢复 G3 曲线；**IALMob 继续保持未授权**。

## 可复现交付物

- 审计器：`scripts/audit_templates_ldmos_g3_continuity_sg_contact.py`
- 回归测试：`tests/regression/test_templates_ldmos_phase23.py`
- ignored 主汇总：
  `reference_staging/ldmos_g3_align_20260827/continuity_sg_contact_v4/summary.json`
- strict reclose nonlinear trace：同目录 `vela_strict_reclose/newton_iterations.csv`
- strict 最后状态回放：同目录 `vela_strict_reclose/final_state_replay/`

审计器单元测试与既有 phase-2/3 工具测试共 20 项全部通过。
