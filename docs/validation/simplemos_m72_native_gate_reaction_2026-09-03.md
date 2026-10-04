# SimpleMOS M72 原生栅反力与完整界面响应

M72 分类为 `native_charge_qualified_full_profile_material_but_not_dominant`。Vela 原生 Poisson Dirichlet 反力不再使用直接 P1 梯度代理；32 个冻结状态观察器回放保持状态逐点不变。

原生微分栅电荷相对 Sentaurus `ContactCharge` 的中位相对误差为 `1.853%`，最大为 `2.239%`，资格结论为 `True`。

完整 Si/SiO2 界面分布的 NWell 配对代理对 M65 剩余电流配对增长的中位闭合率为 `32.89%`，相关系数为 `0.2739`。界面逐点栅控响应误差的中位数 / P95 / 最大值分别为 `0.00107281` / `0.00266732` / `0.00522976` V。

本任务没有新增 Sentaurus 或 Vela 自洽求解，没有改动 HFS、SG、接触提取、准费米打包、BGN、Nc/Nv、网格、收敛参数或生产默认值。原生栅电荷用于观察器资格；完整界面势响应用于静电归因，但两者都不是局部电流连续性分解。
