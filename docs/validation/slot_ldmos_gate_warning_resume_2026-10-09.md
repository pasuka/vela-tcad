# SLOT 门压 warning-only 有界续算（2026-10-09）

本批从已合格诊断态 Vg=0.4962890625 V、Vd=0 V 导入，先完整 Newton 复验；导入与复验不计新增。

用户授权窗口：09:07:03–12:07:03 北京时间，准备、求解与复验共三小时。只建立门压至 5 V。

0.1 V 仅触发日志 warning；自适应初始/最大跨度 0.25 V，成功增步上限 1.5，失败按实际跨度减半；原完整物理及资格不变。

停止原因：gate_5V_qualified；达到 5 V：True；本批新增门压点：19。

最后新增合格点：5.0 V；SG Id：1.5200005573942599e-47 A/µm。

记录/哈希审计：True；这不是独立数值复求，也不是 Sentaurus 定量对比通过。

配置：SLOT 原网格、AverageBox 配套边/体积、Fermi 和接触电离回退；abstol=1e-9、reltol=1e-8、eps_row=.001、闭合=.01、QF 原门槛。Vd=0，参考原点0 V。

证据：build/s09gate_resume 下 preflight.json、summary.json、audit.json、gate_steps.csv、各尝试 config/iv/rows/progress/exit；求解器哈希及输入绑定见 preflight。
