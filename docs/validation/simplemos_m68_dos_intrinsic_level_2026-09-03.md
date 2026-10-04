# SimpleMOS M68 DOS / 本征能级分配审计

M68 分类为 `dos_ratio_inactive_in_frozen_boltzmann_path`。在 `ni` 不变时，把 Si `Nc` 放大 4 倍、`Nv` 缩小到 1/4（若该比值生效，对应本征能级参考移动 `-35.838 mV`），n23 固定状态电流仍逐位一致，差值 `0.000e+00 A/um`、`0.000e+00 dex`。

源码路径同时证明 Boltzmann 浓度、平衡态和 SRH 平衡积不使用 `Nc/Nv`。因此 M65 剩余项不能通过当前 classical noBGN 路径中的 DOS 比值调参闭合；若 Sentaurus 的本征能级参考分配确有影响，它应表现为电势/边界参考差异，而不是 Vela 现有 `Nc/Nv` 参数响应。机器报告：`reference_tcad/simplemos_sentaurus2022/dos_intrinsic_level/m68_dos_intrinsic_level_report.json`。
