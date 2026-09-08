# 浅沟道原生固定电荷双幅度验证：输入已冻结，远程执行待授权

> 后续状态：用户已明确授权，五次 DC 已完成且双幅度校准通过。见同目录 `simplemos_native_channel_charge_fd_execution_2026-09-05.md`；以下保留授权前的准备记录。

2026-09-05。本轮已完成输入制作、空间支撑核对、归档逐文件校验和冻结输入验证；**尚未上传，也未运行任何新的 Sentaurus DC 计算**。不能据此宣称浅沟道 Green 权重已通过原生有限差分校准。

## 冻结实验

工作点为 n23，Vd=0.05 V、Vg=0.9 V，沿用 M81 网格、状态、物理模型及 DC 重闭合设置。共五次 DC：零扰动，以及 ±1e13、±5e12 cm^-3 固定体电荷；不含新 AC 或偏压扫描。

依据本地 Sentaurus Device User Guide T-2022.03 第 548 页公式 523，使用 `SpatialShape=Uniform`、`SpaceMid`、`SpaceSig` 定义窗口。中心为 (0.030045938170817486, 0, 0) um，半宽为 (0.025000000001, 0.125000000001, 1) um。窗口选中原浅沟道集合的全部 186 个节点，与冻结权重积分支撑完全一致，且不含 Poisson Dirichlet 节点。边界所加 1e-12 um 裕量不改变节点成员。

运行后将逐个 Si 节点核验 `ΔSpaceCharge - Δp + Δn - ΔND + ΔNA` 是否等于有符号固定电荷乘以冻结掩码。此处导出 SpaceCharge 使用 cm^-3 数密度，需同时检验窗口内部和外部，不能只凭指令文本判断实际支撑正确。

## 判据与预期

- 五次计算退出码均为零，运行版本为 T-2022.03-SP2，偏压与冻结值偏差不超过 1e-12 V。
- 零扰动电流相对原 M81 电流变化不超过 1e-5 dex。
- 正负中心差分与原生 Green 预测同号、相对差不超过 5%；双幅度归一化导数变化不超过 5%；偶次非线性比例不超过 1%；信号/零控制漂移不小于 100。
- 所有 Si 节点必须通过实际电荷掩码校验，容限及浮点舍入处理已写入 contract.json。

在 +1e13 cm^-3 下，原生 Green 预测为 8.207584787557664e-12 A/um，严格合格 Vela 状态的预测为 9.285954416133963e-12 A/um。这些是待校准预测，当前没有本轮原生有限差分结果。

通过后也仅能资格化本工作点的浅沟道体电荷响应，不能外推至界面片电荷、连续性源项、n19 或高 Vd；M82/M83 仍未放行。

## 载荷和目的地

本地归档：`build-release/simplemos_native_channel_charge_fd/input.tgz`，4,210,148 字节，SHA256：

`f8473838e3f0d3f605b7965f9b626800c126f7b4a6ea76ca4aec319ab4f0132d`

拟上传到用户此前使用的 Sentaurus 虚拟机：`sentaurus:/tmp/vela_simplemos_channel_fd_20260905/input.tgz`。只在该新目录解包并执行 `bash run.sh`，依次运行五次 `sdevice m81_des.cmd`，随后取回输出进行电流与字段校验。

归档含 21 个文件：五组各自的 `input_fps.tdr`、`state_des.tdr`、`state_circuit_des.sav`、`m81_des.cmd`，以及一个 `run.sh`。全部 15 份网格/状态文件已证实与此前 M81 n23 输入逐字节相同；五份 deck 仅固定电荷窗口和幅度所在行变化；所有归档字节与本地冻结文件一致。

## 当前阻塞和复现入口

自动审批审查两次拒绝创建远程目录、上传和执行的组合命令；第二次已提交上述逐字节审计作为补充依据。拒绝理由是：已有授权覆盖 M81 两工况，但未明确覆盖这批五工况载荷及新远程目录，审查将网格/状态视为可能敏感数据。因此本轮未再尝试其他传输路径，等待用户对上述具体载荷、目的地和五次 DC 执行确认。

输入合同及掩码：`reference_tcad/simplemos_sentaurus2022/native_channel_charge_fd/{contract.json,node_mask.csv}`。

载荷逐文件审计：同目录 `upload_payload_audit.json`。其中传输和运行状态为本次冻结准备时的快照。

执行入口：`scripts/validate_simplemos_native_channel_charge.py`；`verify` 已通过，取回 `results.tgz` 后依次运行 `unpack`、`export`、`analyze`。新的运行结果尚不存在。
