# SLOT T470p 新冷态漏压自适应续算

授权窗口：2026-10-12北京时间12:11:51至15:11:51，包含准备、续算和配对检查。Vg固定5 V，Vd尝试从0到60 V。

结束原因：`paired_reference_not_confirmed`；60 V目标达成：False。主路径新增1点，尝试1次。

最后合格Vd=0.006 V，SG Id=9.650052572455e-07 A/µm，rho=1.855905e-10。

## 固定配置与验收

继承新冷态5 V/0 V终态及全部物理。UMFPACK严格消元，浓度有限步，电势对数阻尼；载流子正则化和闭合方程替换关闭。
RHS≤1e-10或原始归一化更新误差<1，并要求原闭合、边界与正性保护；行比仅报告。零更新RHS检查单列，不等同于OR资格。
6 mV起步，热启动，成功增长最多1.5倍且跨度不超过1 V，失败减半。每目标最多40更新、软1080秒/硬1200秒，均受总截止约束。
人口逆映射相对误差1e-7和线性分量后向误差1e-6为保护；不改变方程。PARDISO配对是同父态不同后端，并非不同初值或独立冷启动。

## 配对检查

```json
[
  {
    "milestone_V": 0.006,
    "result": {
      "name": "ref000",
      "Vg_V": 5.0,
      "Vd_V": 0.006,
      "parent_sha": "d04ba8e3f14c977c93116c300e4d053f04d70baadb0d750da2597e1b8c4c5f37",
      "parent_state": "D:\\code-repo\\vela-bench\\s12drain_20261012\\inputs\\parent.h5",
      "span_V": 0.006,
      "backend": "pardiso_mkl_eliminated",
      "role": "same_parent_backend_check",
      "qualified": false,
      "seed_boundary_residual": 0.0,
      "seed_sha": "e748d23c410561012275f25d6c47b225572f555bfa1588cf02961e77b1e97850",
      "error": "AssertionError(('linear_backward_guard', 0.8598548430395853))",
      "traceback": "Traceback (most recent call last):\n  File \"D:\\code-repo\\vela-bench\\s12drain_20261012\\run.py\", line 39, in attempt\n    z=p.arm(name,seed)\n  File \"D:\\code-repo\\vela-bench\\s12drain_20261012\\population.py\", line 35, in arm\n    c.put(O/(name+'_linear_quality.json'),dict(component_backward=berr,eliminated_dofs=eliminated));assert np.isfinite(berr) and berr<=1e-6, ('linear_backward_guard',berr)\n                                                                                                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nAssertionError: ('linear_backward_guard', 0.8598548430395853)\n",
      "wall_s": 43.78383180000128
    }
  }
]
```

## 证据与限制

证据下载及逐文件哈希核对完成。记录审计不代替独立数值复求；电流复现不代表低密度QF全场复现，路径达成不等同Sentaurus定量验证。
本地build/s12drain/results及completion_audit.json；远端D:/code-repo/vela-bench/s12drain_20261012。
