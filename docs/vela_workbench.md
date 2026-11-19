# Windows 仿真工作台原型

工作台通过 `sd-sim.exe` 运行只读算例库。25 个可运行输入对应 23 项 CTest，
包括 PN2D 零偏、正向 I–V、反向 BV 和 16 个变体，以及三个栅极算例、肖特基、
BJT M1 和 BVMethods。栅极的三个输入分别显示在列表中。
界面默认使用英文，可通过「设置」菜单或工具栏的语言图标切换。应用会保存显示语言。
「文件」「仿真」「设置」菜单下方的图标工具栏提供相同操作，悬停可查看说明。
用 Ctrl 或 Shift 多选算例，从「仿真」菜单或工具栏选择「运行选中」开始串行批跑。
「Stop Queue」取消等待中的算例，当前求解结束后可以关闭窗口。
单跑时进度条以动画表示求解仍在进行。批跑时进度条按已处理的算例计数，
状态栏从运行日志显示当前阶段和最后接受的偏压。
点击算例可查看曲线、三角网格、末态场量、日志及只读的 Sentaurus 原始脚本。
所有算例的曲线默认以端子内压为横轴，曲线页可切换到外压。旧结果缺少内压时，
页面会提示可选用外压。网格页支持顺时针或逆时针旋转 90° 和复位，同一次运行的
末态云图跟随旋转，坐标轴以 x′、y′ 标记旋转后的方向。
查看已有结果时，先选择对应算例，再点击「打开已有运行」并选择运行目录。
工作台核对算例输入哈希，随后读取已保存的结果。首次打开有结果的图形页时，
工作台才加载绘图库，Matplotlib 缓存存放于 `%LOCALAPPDATA%/VelaTCAD/matplotlib`。
也可用 `SDS.exe --open-run <算例 ID> "<运行目录的绝对路径>"` 直接打开。
内部性能验收可将 `VELA_WORKBENCH_TIMINGS` 设为可写的 JSON Lines 文件路径。

## 用户配置

首次启动会创建 `%LOCALAPPDATA%/VelaTCAD/sds-settings.json`。从「设置」菜单选择
「打开配置文件」可编辑此文件。修改后保存并重启工作台。也可用
`SDS.exe --config "<配置文件的绝对路径>"` 读取另一份配置。

```json
{
  "language": "en",
  "font_size": 10,
  "log_level": "debug"
}
```

- `language` 支持 `en` 和 `zh-CN`。
- `font_size` 支持 9 至 16 磅的整数，用于控件、日志、脚本和绘图文字。
- `log_level` 支持 `debug`、`info`、`warning` 和 `error`，控制每次新运行的
  `*_des.log`。`debug` 保留原有完整输出，`info` 省去常规求解器明细，
  `warning` 和 `error` 进一步限制消息等级。各级均保留运行头和总结。
  `sd-sim.log` 继续记录控制台进度和错误，单独的阶段诊断日志仍可在运行目录查看。

配置值无效时，界面会提示并使用默认值，原文件会保留供修改。

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

## Ubuntu 26.04 源码启动

在 Ubuntu 原生环境中构建 Linux 版求解器。图形界面需要 Python 的 Tkinter、
Pillow ImageTk 和 Matplotlib。Ubuntu 系统软件包分别为 `python3-tk`、
`python3-pil.imagetk` 和 `python3-matplotlib`。中文绘图需要系统有支持
中日韩文字的字体，例如 `fonts-noto-cjk`。

```bash
cmake -S . -B build-linux26 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DVELA_PHYSICAL_CONSTANTS_PROFILE=sentaurus_u2022_12
cmake --build build-linux26 --target sd-sim vela_example_runner sentaurus_import --parallel 2
python3 scripts/vela_workbench.py --stage-cases build-linux26/sds-cases
python3 scripts/vela_workbench.py --cases build-linux26/sds-cases \
  --launcher build-linux26/sd-sim --runs build-linux26/sds-runs
```

在桌面会话中运行最后一条命令。自动化也可以使用 `xvfb-run` 测试 Tk 窗口。
Linux 命令行运行包可按 [SDevice 运行包指南](sdevice_packaging.md) 生成。
Windows ZIP 使用 PE 可执行文件和 DLL。Linux 图形发行包需要单独打包并验证
ELF 程序、共享库及 Tk 图形依赖。

## Windows ZIP

求解器使用 MSYS2 UCRT64 Release 构建。GUI 和 SDevice 前端用独立的
Windows x64 CPython 3.12 构建，需在该解释器中安装 Nuitka 和 Matplotlib。
Nuitka 首次构建时会下载 Zig 工具链。设置 `$nuitkaPython` 为此
解释器路径，设置 `$mklRoot` 为本机 oneMKL 2026.1 安装目录。

```powershell
cmake --preset windows-ucrt64-release "-DVELA_SDS_NUITKA_PYTHON=$nuitkaPython" "-DVELA_SDS_MKL_ROOT=$mklRoot"
cmake --build --preset windows-ucrt64-release --target sds_package
```

配置时提供这两个路径即可启用 `sds_package`。该目标按需执行，并先构建
求解器和导入工具。产物为 `build-release/sds-package/SDS-windows-x64.zip`
和同目录下的 `SDS/SDS.exe`。再次执行时，通过检查的新包会替换该目录中的旧包。

解压整个 `SDS-windows-x64.zip` 后，双击其中的 `SDS.exe`。目录包带有
已编译的 GUI 和 SDevice 前端、Tk、绘图库、`sd-sim`、只读算例库和
oneMKL 顺序运行时。打包程序会检查最终 ZIP 及其内嵌 ZIP 格式归档中的
Python 源码、字节码和 zipapp 文件。Sentaurus 原始输入仍作为算例资料提供。
启动器为 BJT M1 和 BVMethods 设置随包 PARDISO 所用的库路径。内部发行的
测试记录应注明使用的 `sd-sim` 版本和验证过的 Windows 环境。
