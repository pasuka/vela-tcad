# SLOT Sentaurus阻尼及初始化顺序核查

2026-10-26通过vela VM只读核对原Workbench项目的BVDS_des.cmd、IdVd_des.cmd、IdVd1_des.cmd、IdVd2_des.cmd，并与已运行的低偏压诊断脚本和日志比较。未修改脚本、未启动求解。具体VM路径保留于用户本地运行记录。

## 原Workbench设置

| 项目 | IdVd_des.cmd | BVDS_des.cmd |
| --- | --- | --- |
| Math.NotDamped | 25（第93行） | 100（第88行） |
| Math.Iterations | 25（第94行） | 25（第89行） |
| 初始Poisson.Iterations | 100（第117行） | 1000（第106行） |
| 初始Poisson.LineSearchDamping | 1e-8（第117行） | 1e-8（第106行） |
| 初始耦合及扫描Coupled.LineSearchDamping | 未显式设置，默认1 | 未显式设置，默认1 |
| 初始耦合方程 | Poisson Electron Hole，可选Temperature | Poisson Electron Hole |
| 后续流程 | 耦合升栅压→Save/Load→耦合扫漏压 | 耦合扫漏极目标电压 |

IdVd_des.cmd第117—118行和BVDS_des.cmd第106—107行都明确先求Poisson，再求电势与载流子连续性方程耦合系统。后续偏压点继续求解耦合系统，不是在每个偏压点重新执行最前面的独立Poisson初始化。

原项目另有IdVd1_des.cmd及IdVd2_des.cmd，两份的初始Poisson设置同为Iterations=100、LineSearchDamping=1e-8，但后续Coupled只列Poisson Electron及可选Temperature，没有Hole。不能用这两份单载流子变体替代主IdVd_des.cmd的双载流子流程。

## 参数语义

VM保存的Sentaurus Device User Guide U-2022.12第206页及参数表确认：

- LineSearchDamping默认1，线搜索默认关闭；指定小于1时启用线搜索，此值是允许的最小阻尼系数，实际每步因子自动决定。1e-8不表示每步固定只更新一亿分之一。
- NotDamped是Bank–Rose阻尼启动前允许的迭代数量；超过该值才启动该阻尼。它既不是固定步长，也不是残差或更新量阈值。
- Math设置为Coupled提供默认值；具体Coupled括号中的设置可覆盖该阶段默认值。

因此原脚本普通耦合阶段最多25步，IdVd.NotDamped=25、BVDS.NotDamped=100，通常在达到Bank–Rose启动条件前就已经结束该次耦合尝试。初始Poisson的独立局部控制与后续扫描控制不能混为一谈。

## 已运行诊断脚本的差别

| 诊断步骤 | 实际Solve设置 |
| --- | --- |
| BVDS独立Poisson采集 | Coupled(Iterations=100 LineSearchDamping=1e-4){Poisson} |
| BVDS冷态耦合采集 | Load(Poisson保存态)，再Coupled(Iterations=25){Poisson Electron Hole} |
| 补偿掺杂IdVd冷态采集 | Coupled(Iterations=1000 LineSearchDamping=1e-8){Poisson}，再Coupled(Iterations=25){Poisson Electron Hole} |
| IdVd升栅及5 mV采集 | Load前一阶段保存态，再执行Quasistationary中的耦合求解 |

诊断日志的全局控制回显分别为BVDS：LineSearchDamping=1、undamped iterations=100；IdVd：LineSearchDamping=1、undamped iterations=25。局部Poisson的参数覆盖写在Solve中，不能只读全局回显便认为Poisson未启用线搜索。

这些拆分脚本保留先Poisson、再耦合的状态依赖，但Poisson局部迭代上限及BVDS线搜索下限与原Workbench不完全相同。它们是独立诊断配置，不应被称为逐参数完全复现原脚本。
