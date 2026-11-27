# SLOT T470p：全零偏压新冷态与门压建立

用户授权三小时窗口：2026-10-12 北京时间10:46:46–13:46:46，含准备、求解与复核。
Vd固定0 V，从全零偏压冷平衡建立Vg至5 V；无旧状态导入，不包含漏压扫描。

## 结果

- 控制器结束：True；5 V目标通过：True；原因：`new_cold_gate_5V_qualified`。
- 合格冷态 1，合格门压点 28。
- 最新合格Vg=5 V、Vd=0 V，SG Id=-2.525984272734e-32 A/µm，整体残差=4.556825e-11。

## 配置与资格

冻结最新s12材料/物理配置：AverageBox、定义积分Fermi、effective_band_sg、SRH/Auger/电离及接触回退。
UMFPACK严格消元主路径；载流子正则化和闭合替换关闭，闭合最终验收保留，行比仅报告。
门压5 mV起步，成功按完整更新次数自适应增长、上限0.25 V，失败减半，最多8次连续失败。
使用隔离多晶硅等势平移预测→Poisson→完整器件检查；Poisson态和预测态不计合格。
必要时最多8步浓度有限步修正，保留1e-7浓度恢复保护；本批还要求零更新重读通过RHS检查，故比仅RHS或更新误差更严格。
零漏压闭合可经绝对误差分支通过；并不等于每条载流子行达到1e-3。
最终不同门压路径/PARDISO复核若未完成或未通过，不称为独立复现成功。

## 复核与旧基准

```json
{
  "mid": {
    "qualified": true,
    "Vg_V": 4.9083251953125,
    "Vd_V": 0.0,
    "state": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\independent_mid_zero\\audited.h5",
    "sha": "1051baef3514ec13557e6eea91770e7d7d1b688b4fa507308295a5fe2a79ad6e",
    "updates": 0,
    "origin": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\independent_mid_poisson\\state.h5",
    "rho": 6.080882800339252e-11,
    "closure": {
      "electron": {
        "contact_flux": 2.8353664969072735e-27,
        "integrated_source": -7.680739487541805e-33,
        "mismatch": 2.835374177646761e-27,
        "qualified": false,
        "ratio": 2.8353741776467614e-16
      },
      "enabled": true,
      "enforced": true,
      "hole": {
        "contact_flux": 8.029179820452637e-28,
        "integrated_source": -7.680739487541805e-33,
        "mismatch": 8.029256627847513e-28,
        "qualified": false,
        "ratio": 8.029256627847513e-17
      },
      "satisfied": true
    },
    "rows": {
      "violation_count": 10647,
      "max_ratio": 1.0000923011153147,
      "satisfied": false,
      "enforced": false
    },
    "Id_A_per_um": -2.525984272734077e-32,
    "boundary_residual": 0.0,
    "backend": "pardiso_mkl_eliminated",
    "parent_sha": "4071697076088e85c76ac24a3f75c2fab3fce53f8a19bc9a6c8c54f011f10d35",
    "parent_state": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\g026_zero\\audited.h5",
    "span_V": 0.0916748046875,
    "wall_s": 80.85243479999917
  },
  "final": {
    "qualified": true,
    "Vg_V": 5.0,
    "Vd_V": 0.0,
    "state": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\independent_final_zero\\audited.h5",
    "sha": "ae3d02eae4def948492a2b04298ef5707232539189eaa3904dc079a2f6d92613",
    "updates": 0,
    "origin": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\independent_final_poisson\\state.h5",
    "rho": 5.011244587578822e-11,
    "closure": {
      "electron": {
        "contact_flux": 2.8353690539123782e-27,
        "integrated_source": -7.709466184382585e-33,
        "mismatch": 2.8353767633785625e-27,
        "qualified": false,
        "ratio": 2.8353767633785627e-16
      },
      "enabled": true,
      "enforced": true,
      "hole": {
        "contact_flux": 8.029179884607833e-28,
        "integrated_source": -7.709466184382585e-33,
        "mismatch": 8.029256979269677e-28,
        "qualified": false,
        "ratio": 8.029256979269677e-17
      },
      "satisfied": true
    },
    "rows": {
      "violation_count": 10647,
      "max_ratio": 1.0000926267130545,
      "satisfied": false,
      "enforced": false
    },
    "Id_A_per_um": -2.525984272734077e-32,
    "boundary_residual": 0.0,
    "backend": "pardiso_mkl_eliminated",
    "parent_sha": "1051baef3514ec13557e6eea91770e7d7d1b688b4fa507308295a5fe2a79ad6e",
    "parent_state": "D:\\code-repo\\vela-bench\\s12coldgate_20261012\\independent_mid_zero\\audited.h5",
    "span_V": 0.0916748046875,
    "wall_s": 80.79369100000076
  },
  "fields": {
    "psi": {
      "max_abs": 7.105427357601002e-15,
      "node": 13452
    },
    "phin": {
      "max_abs": 0.0,
      "node": 0
    },
    "phip": {
      "max_abs": 0.0,
      "node": 0
    },
    "electrons_m3": {
      "max_abs": 1133871366144.0,
      "node": 5359
    },
    "holes_m3": {
      "max_abs": 1889785610240.0,
      "node": 12551
    }
  }
}
```

旧5 V基准为s09gate_resume/g020，同为Vd=0。两者存在配置差异，字段差异不直接等于旧态误差。

```json
{
  "psi": {
    "max_abs": 5.653233536295765e-06,
    "node": 10237,
    "new": 0.607249127446743,
    "old": 0.6072434742132067
  },
  "phin": {
    "max_abs": 2.7712637859103416e-20,
    "node": 5408,
    "new": 2.7712637859103416e-20,
    "old": -4.867137158745196e-44
  },
  "phip": {
    "max_abs": 2.721383502503687e-20,
    "node": 5408,
    "new": -2.721383502503687e-20,
    "old": 7.48126699378585e-48
  },
  "electrons_m3": {
    "max_abs": 1.8810422284580626e+22,
    "node": 10237,
    "new": 4.2886638181140266e+26,
    "old": 4.288475713891181e+26
  },
  "holes_m3": {
    "max_abs": 1.4794801689080627e+17,
    "node": 13291,
    "new": 2.17481614309588e+26,
    "old": 2.1748161416164e+26
  }
}
```

## 证据

本地证据：`build/s12coldgate/package/`、`results/`、`completion_audit.json`、`archive.json`。
远端目录：`D:/code-repo/vela-bench/s12coldgate_20261012`。下载归档及逐文件哈希核对通过。
记录/哈希审计不是独立数值复求。实际字段单位为V、m⁻³，端口电流A/µm。
h5py/HDF5版本不一致警告保留；本次未升级依赖或修改全局默认。
