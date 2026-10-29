# SLOT IdVd 显式外推驱动集成短验（2026-10-06）

## 范围

此前三个相邻点与一个同目标不同外推初值复求已经通过；纯热启动复求失败。
本批接续完成显式驱动集成，使用实际 CLI 推进两个新漏压目标，不重复上一批任务，
不启动至 60 V 的长扫描。本批两个新目标均一次通过原资格，各 4 次 Newton 更新，
约 326 秒完成，原资格及文件绑定核对通过。

工作树：`D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next`。
Vg=5 V、300 K、物理参考原点 0 V，电压单位 V、电流 A/µm、状态密度 m⁻³。
沿用 13,881 节点 SLOT 网格、Fermi-Dirac/Bennett-Wilson、Si PhuMob/Lombardi/GradQF HFS、
复合/雪崩、Poly 材料分流和 element-box 配置。配置来自上一批 `p03_secant/config.json`，
仅变更偏压目标、初值和输出位置；非输出物理配置逐项核对。
冻结求解器为 Release/HDF5/UMFPACK，SHA-256：
`d137894b21332b59f20afedba08e99daf334d6f0a931a517a54b38401b920c35`。

原 abstol=1e-9、reltol=1e-8、eps_row=0.001、全局闭合 0.01、QF 检查不变。
门压初始化与漏压外推分别控制，没有将 Fermi 多晶硅门压初始化套入漏压路径。

## 实现

- `--drain-initialization secant` 显式开启；默认仍为 `warm`。
  采用已验证的分离参考量/增量 QF 与对数密度猜测，完整 Fermi Newton 完成资格验收。
- 历史只包含合格漏压态，外推倍率限定 (0,3]；倍率越界直接使用热启动。
  外推求解失败后只增加一次同目标热启动，仍失败则进入原最多三层二分。
  未合格态不能更新父态或外推历史。
- `--drain-import-summary` 从已完成批次最后两份合格漏压记录导入，核对实际 IV、
  配置/状态哈希、物理输入、接触偏压和原资格；导入数量单列。这是保存证据复核，不是独立复求。
- `--drain-stop` 限定终点为原生精确样本。短批 `goal_reached` 仅表示请求区间完成，
  `stop_reason=requested_drain_segment_qualified`，不能解释为 60 V 全路径完成。
- 策略、冻结工具、求解器及导入证据哈希写入新目录身份记录；旧生产账本不追加。

## 本批命令与预算

从上一批 Vd=0.172651318567887 / 0.210912960216741 V 两个合格态导入，
新目标为 0.256418806017777 / 0.310540425290476 V。
本批总预算 1800 秒，普通子求解上限 600 秒，串行运行。

```powershell
D:/msys64/ucrt64/bin/python.exe -B scripts/run_slot_ldmos_idvd.py --config build/s06idvd_neighbors/p03_secant/config.json --runner build/s06idvd_neighbors/runner.exe --native build/s05idvd/native_drain_path.csv --output build/s06idvd_integration --budget 1800 --child-cap 600 --drain-initialization secant --drain-import-summary build/s06idvd_neighbors/summary.json --drain-stop 0.310540425290476
```

执行文件与状态位于 `build/s06idvd_integration/`；该命令要求新输出目录，不能覆盖本次证据重跑。

## 验证

### 数值结果

| Vd / V | SG Id / A·µm⁻¹ | 全部尝试成本 / s | Newton 更新 | 最大行比 | 电子 relative 闭合比 |
|---:|---:|---:|---:|---:|---:|
| 0.256418806017777 | 3.9853440131107896e-5 | 162.54 | 4 | 1.8882639e-7 | 2.4663930e-4 |
| 0.310540425290476 | 4.7858992174342794e-5 | 163.30 | 4 | 5.5233248e-7 | 3.2648970e-4 |

两个点的有效载流子超标行数及 QF 越界均为 0，以 `abstol` 收敛；
空穴 relative 闭合比分别为 6.16980e-15 / 9.65544e-15。
外推倍率均为 1.189333；预测成本分别 0.312 / 0.333 秒。
全部尝试共 325.83 秒，批次总墙钟 325.90 秒，未触发重试或缩步。
平均 162.92 秒/点不是成对热启动对照所得的加速倍数，也不能外推为高漏压成本。

同偏压原生电流分别为 3.97702776378211e-5 / 4.77593541876567e-5 A/µm，
本批 Vela 分别高 0.209107% / 0.208625%，取原生精确样本，不插值。
这延续低漏压区的差异水平，尚非全 IdVd 曲线精度结论。

`audit.json` 保存的资格与绑定核对结果为通过；两份导入态单列，本批新增 2 个合格态。
前批及原生产账本文件哈希保持不变。当前代码与预测器、材料、网格、冻结求解器、
初值和终态均核对绑定；没有用预测态的资格标志代替完整器件证明。

### 回归

以下 9 项回归通过：

```powershell
D:/msys64/ucrt64/bin/python.exe -B -m unittest tests.regression.test_run_slot_ldmos_idvd tests.regression.test_slot_idvd_secant_predictor tests.regression.test_slot_idvd_configuration
```

包含原门压/Poisson 隔离、失败父态隔离、外推失败→热启动失败→二分的历史与记账、
导入状态哈希不符拒绝，以及预测器和材料配置检查。控制流程测试使用模拟求解器，
不作为数值证据。首次测试因沙箱临时目录写权限失败，改用工作树下 `build/s06idvd_test_tmp`
作为 TEMP/TMP 后全部通过；没有调整求解门槛。未修改 C++，未重新编译。

`audit.py` 用于本批结束后的原资格、IV、配置、状态和资源绑定核对，
不能替代独立数值复求或网格/模型验证。本批不新增同目标独立复求；
上一批成功的复求范围仅为 Vd=0.210913 V。

## 结论与后续边界

显式外推已接入实际驱动并通过两个低漏压相邻目标的短验。可作为下一段有界漏压续算的
初值策略，继续周期性抽样同目标复求，分别记录新增点、失败与回退成本。
本批仅模拟回归覆盖失败后的热启动/二分控制链；数值短验没有触发该分支，
不能据此宣布实际失败恢复可靠，也不覆盖高漏压或雪崩区。
本批已正常结束，没有自动扩展预算、启动长扫描或恢复定时监测。
