# SimpleMOS M65 no-BGN 本征浓度归因

## 结论

M65 分类为 `base_intrinsic_density_convention_material_but_not_dominant`。Sentaurus 显式 no-BGN 的硅区 `EffectiveIntrinsicDensity` 为 `10750038488.84424 cm^-3`，Vela 原基准为 `14638914958.76762 cm^-3`；M65 只在复制材料文件中把 Si `ni` 匹配到前者，未修改生产材料、BGN、HFS、SG、接触、网格或收敛参数。

8 个冻结 NWell 配对的误差增幅由中位 `0.111220 dex` 变为 `0.037029 dex`，中位闭合率 `68.05%`，其中 `0/8` 个配对达到 80% 闭合。

状态账本覆盖 16 个固定诊断点。Sentaurus `ni` 最大相对空间展宽为 `0.000e+00`；匹配后的电子浓度差由 `ni + (psi-phin)` Boltzmann 恒等式闭合，p95 残差 `1.827e-06 dex`。这使电势/准费米势垒账本可用于解释匹配后剩余项，而不是把其误记为基础 `ni` 差异。

机器报告：`reference_tcad/simplemos_sentaurus2022/nobgn_intrinsic_density_attribution/m65_nobgn_intrinsic_density_attribution_report.json`。
