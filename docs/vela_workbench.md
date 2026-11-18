# Windows 仿真工作台原型

工作台通过 `sd-sim.exe` 运行只读算例库。25 个可运行输入对应 23 项 CTest，
包括 PN2D 零偏、正向 I–V、反向 BV 和 16 个变体，以及三个栅极算例、肖特基、
BJT M1 和 BVMethods。栅极的三个输入分别显示在列表中。
界面默认使用英文，工具栏「中文」和「English」可随时切换。
用 Ctrl 或 Shift 多选算例，点击「Run Selected」开始串行批跑。
「Stop Queue」取消等待中的算例，当前求解结束后可以关闭窗口。
单跑时进度条以动画表示求解仍在进行。批跑时进度条按已处理的算例计数，
状态栏从运行日志显示当前阶段和最后接受的偏压。
点击算例可查看曲线、三角网格、末态场量、日志及只读的 Sentaurus 原始脚本。

## 源码启动

从当前 worktree 根目录执行。开发机需要 MSYS2 UCRT64 Python、
Tkinter、Matplotlib 及准备参考算例所用的回归依赖。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target sd-sim vela_example_runner sentaurus_import
python scripts/vela_workbench.py --stage-cases build/vela-workbench-cases
python scripts/vela_workbench.py --cases build/vela-workbench-cases --launcher build-release/sd-sim.exe --runs build/workbench-runs
```

已有 `sentaurus_u2022_12` 常数配置的 Vela runtime 时，可以将
`--launcher` 指向该运行包的 `sd-sim.exe`，直接试用原型。

`--stage-cases` 使用与 CTest 相同的输入准备逻辑。工作台额外传入
`--plot-final-only`，省去每个接受点的常规云图输出。命令要求新的目录，
可用不同路径重新准备算例。省略 `--runs` 时，运行记录存放于
`%LOCALAPPDATA%/VelaTCAD/runs`。每次运行创建独立目录并校验输入哈希。
脚本页展示仓库中的 Sentaurus 源文件或参考输入，供阅读和复制。
工作台求解使用回归准备的 GRD/DAT 输入。已有算例库可继续运行，
重新准备算例库后即可在界面查看原始脚本。
「达到目标」来自退出码、阶段计划与 `run_manifest.json`，参考回归仍通过
CTest 执行。

分阶段的算例按扫描端子绘制电流曲线。云图读取末态 ASCII Tri3 VTU。
BVMethods 达到漏极电流停止阈值即为正常结束，工作台单独显示该状态。

## Windows ZIP

在 MSYS2 UCRT64 Release 构建后安装 [SDevice runtime](sdevice_packaging.md)。
设置 `$mklRoot` 为本机 oneMKL 2026.1 安装目录，然后执行：

```powershell
cmake --build --preset windows-ucrt64-release --target sd-sim vela vela_example_runner sentaurus_import
python scripts/package_vela_workbench.py --build-dir build-release --mkl "$mklRoot" --output build/sds-release
```

解压整个 `SDS-windows-x64.zip` 后，双击其中的 `SDS.exe`。目录包带有
Python、Tk、绘图库、`sd-sim`、只读算例库和 oneMKL 顺序运行时。
启动器为 BJT M1 和 BVMethods 设置随包 PARDISO 所用的库路径。内部发行的
测试记录应注明使用的 `sd-sim` 版本和验证过的 Windows 环境。
