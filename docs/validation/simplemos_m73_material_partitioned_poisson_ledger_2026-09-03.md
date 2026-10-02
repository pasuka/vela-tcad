# SimpleMOS M73 材料分区 Poisson 账本

M73 分类为 `material_partition_material_but_not_dominant`。任务复用 M65/M69 的 32 个 drain/gate 冻结状态，没有新增 Sentaurus 或 Vela 自洽求解。

五种离线 Poisson 观察器中，对剩余电流配对增长解释最好的材料分区候选为 `region_local_barycentric_si`：同号 `8/8`，中位闭合率 `73.57%`。对 M65 静电势垒配对代理解释最好的是 `legacy_signed_si`：同号 `8/8`，中位闭合率 `87.40%`。生产 `legacy_all_cell` 冻结修正本身对电流目标的中位闭合率为 `64.46%`。

Vela 生产状态回放的最大 Poisson 线性修正为 `4.913e-09` V，分项线性叠加最大误差为 `1.033e-12` V。上述结果只用于固定状态归因；未修改 HFS、SG、接触提取、准费米打包、BGN、Nc/Nv、网格、收敛参数或生产默认值。
