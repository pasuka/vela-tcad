# SLOT LDMOS · 2026-09-21 日报素材包

- [日报文字与配图说明](report.md)
- [器件结构与验证范围 PNG](figures/01_device_structure.png) / [SVG](figures/01_device_structure.svg)
- [10 V 电势与验证进展 PNG](figures/02_potential_and_progress.png) / [SVG](figures/02_potential_and_progress.svg)
- [雪崩边界差异定位 PNG](figures/03_avalanche_boundary.png) / [SVG](figures/03_avalanche_boundary.svg)
- [指标、输入来源与 SHA-256](data/summary.json)

图像由真实网格、已合格 10 V 保存态和同原生状态源项审计生成。
`data/figure_data.npz` 仅保存重绘所需的几何/场量快照，不是求解重启包。
SVG 的文字/坐标标注为矢量轮廓，密集网格场量为 220 dpi 嵌入图层。
PNG 尺寸为 3520×1980，可直接插入日报。

在仓库根目录使用已安装 NumPy/Matplotlib 的 Python 重绘：

```powershell
$env:MPLCONFIGDIR = Join-Path (Get-Location) 'build/mpl-slot-daily'
D:/msys64/ucrt64/bin/python.exe -X utf8 scripts/plot_slot_ldmos_daily_report.py
```

Windows 自动使用微软雅黑。其他环境需安装兼容中文字体并调整字体列表。
只有显式添加 `--prepare` 才会重新读取本机 `build/` 下的实验输入；通常
重绘不需要该选项。脚本只读数据、绘图，不运行仿真或连接虚拟机。

代码状态：`076ab59c`。完整 Release 回归 888/888 通过；此素材包只是
当日验证快照。高压及 BVDS 未完成，仿真继续保持暂停。
