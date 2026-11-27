# SLOT IDL：Bennett–Wilson 模型映射与接入验证

日期：2026-10-05。工作树：`codex/slot-ldmos-bvds-next`。
本批补齐显式 `bennett_wilson` 选项与本算例配置，保留旧批次证据及原验收门槛。
模型/导入回归和真实网格零偏压接入检查通过；正式非零漏压 IDL 对比仍未完成。

## 依据与公式

[Sentaurus Device User Guide O-2018.06（镜像）](https://studylib.net/doc/27467177/sdevice-ug---%E5%89%AF%E6%9C%AC)
第 12 章，式 159–165、表 36–37给出带隙、亲和势及 Bennett–Wilson 定义。
[Synopsys 培训材料 17.5.2（镜像）](https://ghzphy.github.io/Sentaurus_Training/sd/sd_17.html#5.2)
独立列出同一公式与总掺杂口径：

\[
N_{\rm tot}=N_D+N_A,\qquad
\Delta E_g^0=\begin{cases}
E_{\rm ref}\,[\ln(N_{\rm tot}/N_{\rm ref})]^2,&N_{\rm tot}\ge N_{\rm ref},\\
0,&N_{\rm tot}<N_{\rm ref}.
\end{cases}
\]

对数是自然对数。不是净掺杂绝对值，也不取自由电子/空穴密度的最大值。
本批无 Fermi 修正；对称带边移动下，`ni_eff = ni_base exp(DeltaEg/(2 Vt))`。
正窄化项降低有效带隙。基础带隙参考修正 `dEg0` 与上述窄化项必须分别处理，不能重复加入。

原始文献为 H. S. Bennett、C. L. Wilson，*Statistical comparisons of data on
band-gap narrowing in heavily doped silicon: Electrical and optical measurements*，
J. Appl. Phys. 55(10), 3582–3587 (1984)，DOI `10.1063/1.332950`。
[NBS/NIST 1984 文献目录](https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nbsspecialpublication305supp16.pdf)
记录了该论文及摘要。本次没有取得论文全文；实现依据是可核对的 Sentaurus 公式与实际参数文件，
不把这条经验式宣称为普适第一性原理模型。

## 本算例的显式映射

原生证据：`build/s05idl/native_resolved/console.log:335` 明确打印
`Bennett/Wilson ... (no Fermi)`；同目录 `idl.par` 给出下表。

| 原生定义 | Vela 配置或处理 | 本算例值 |
|---|---|---|
| BennettWilson | `solver.bandgap_narrowing.model` | `bennett_wilson` |
| Bennett.Ebgn | `coefficient_eV` | 0.00684 eV |
| Bennett.Nref | `reference_doping_m3` | unit_scaling 输入 3.162e18 cm⁻³ |
| SI 下同一参考浓度 | 同一字段，legacy_si 输入 | 3.162e24 m⁻³ |
| Bandgap.dEg0(Bennett) | 基础材料 ni 的参考修正，非额外窄化 | 0 eV |
| Bandgap.Bgn2Chi | 现有对称带边处理 | 0.5 |
| no Fermi | Boltzmann，`fermi_statistics_correction=false` | 保持 |

注意字段名的历史 `_m3` 后缀不能代替实际 `scaling.mode` 单位约定。
字符串形式 `"bennett_wilson"` 的工厂会自动换算默认参考浓度。

| 总掺杂 / cm⁻³ | 窄化 / eV |
|---:|---:|
| 1e17 | 0 |
| 3.162e18 | 0 |
| 1e19 | 0.00906762876516 |
| 1e20 | 0.0816003608035 |

以上为公式校验点，不是新求解的器件电流数据。

## 实现范围与回归

- `BandgapNarrowing.h/.cpp` 增加独立模型和工厂分支，载流子注入不改变窄化值。
  检查有限正参考浓度、有限非负系数及非负总掺杂；拒绝非零 smoothing/offset，避免混入 Slotboom 或重复参考修正。
- Newton、Gummel 的既有 JSON 解析使用同一工厂；装配沿用总杂质浓度和 variable-ni 输运路径。
- 文本导入器识别显式 BennettWilson。混合 BGN 选择、BennettWilson+Fermi 自动导入拒绝执行。
  未写模型的普通 `.cmd` 仍须核对实际原生默认值；不把所有旧项目静默改成 BennettWilson。
  通用 `.par` 材料参考量/温度律的自动发射仍是独立限制，本批没有宣称其全部完成。
- `none`、Slotboom/OldSlotboom 默认值及 PN2D 模板不改。BV 旧配置及冻结求解器不改。

验证使用 UCRT64 C++20 Release；实际编译定义包含 HDF5 state、SPQR、UMFPACK。
编译 `vela_example_runner` 及 4 个相关测试目标成功。

| 检查 | 结果 |
|---|---:|
| test_recombination 全部 | 43 cases / 234 assertions 通过 |
| test_gummel_high_doping 全部 | 6 / 160 通过 |
| test_sg_flux `[bgn]` | 12 / 199 通过 |
| test_newton_solver `[bgn]` | 5 / 47 通过 |
| test_sentaurus_import_tools.py | 54 tests 通过 |
| test_sentaurus_parameter_map.py | 51 tests 通过 |

新增测试覆盖阈值、文档数值、注入独立性、两种单位下两类求解器的字符串/对象解析、
高补偿掺杂节点的 ni、Gummel 平衡态质量作用关系及非法参数。
未重跑全 CTest；上述结果不引用上一批测试作为本批通过证据。

## 真实网格接入检查与边界

新目录 `build/s05bennett/`；可追溯脚本 `verify_mapping.py`，输入/runner/参考资料
哈希在 `manifest.json`，源码哈希和完整 IV 行在 `analysis.json`。
对旧 `s05idl` 所有文件前后 SHA-256 核对一致。

与旧零偏压 deck 比较：除注释和输出位置外，唯一配置变化为 BGN 从 none 改成上述显式对象。
沿用 SLOT 原网格（13,881 节点、27,275 三角形）、材料、掺杂、外部 AverageBox 输运、
300 K、四个 Ohmic 接触、无雪崩；`unit_scaling`，电流报告 A/µm。
没有载入旧 none-BGN 态，采用冷 Poisson 初始化。

限时 120 s，实际 20.83 s，退出码 0：8 次 Newton 更新、停止原因 `abstol`，
载流子超标行数 0，最大行比 `2.4368851e-8 < 0.001`。
全局闭合 `satisfied=1`，此零偏压点走绝对判据，电子/空穴相对闭合 qualified 标志均为 0；
不可据此声称高电流相对闭合已验证。零偏压电流约 `-1.01e-51 A/µm`，只表示数值平衡。

本次未启动门压长扫描，正式非零漏压 IDL 点仍为 0；未改任何验收门槛。
下一步应从这份 Bennett-Wilson 配置的新平衡态开始，先验证低漏压建立，再进行门压路径，
之后才可与原生 IDL 曲线作定量比较。旧 none 批次失败仍不能用于判断门槛过严。
