"""Chinese translations for the workbench's English interface."""

CASE_NAMES = {"pn2d_ci_0v": "PN2D Zero Bias", "pn2d_ci_iv": "PN2D Forward I–V"}

GUIDE_EN = """SDS User Guide

Run a case
Select a built-in case and choose Run Current. Use Ctrl or Shift to select multiple cases, then choose Run Selected. Stop Queue cancels the cases waiting to run. The current case can finish and save its results.

Open and inspect results
Select the case, then choose File → Open Saved Run and select its run folder. Use the Curve toolbar to choose inner or outer voltage, set either axis to Linear or Log, and zoom around the center with + or −. Log scales omit non-positive values.

In the Mesh and Final Field toolbars, use ↶, ↷ or 0° to rotate or reset both views. In Final Field, use the toolbar dropdown to select a field and turn on Show Mesh to overlay the triangle edges.

Find run files and settings
Use the Log tab for the main run log and the Sentaurus Script tab for the original input. File → Open Run Folder opens the full set of artifacts. Change the language or open the configuration file from Settings.
"""

GUIDE_ZH = """SDS 使用文档

运行算例
选中内置算例后点击「运行当前」。按 Ctrl 或 Shift 可选中多个算例，再点击「运行选中」。「停止队列」会取消等待中的算例，当前算例继续运行并保存结果。

打开和查看结果
选中对应算例，从「文件 → 打开已有运行」选择运行目录。在曲线工具栏可切换内压和外压，单独设置 X、Y 轴的线性或对数坐标，并用 + 或 − 围绕图形中心缩放。对数坐标会省略非正数。

「网格」和「末态云图」页的绘图工具栏中，↶、↷ 和 0° 可同步旋转或复位两页视图。在末态云图页可在工具栏选择物理量，也可勾选「显示网格」叠加三角单元边线。

查找文件和设置
在「日志」页可查看主日志，在「Sentaurus 脚本」页可查看原始输入。从「文件 → 打开运行目录」可查看完整输出。从「设置」菜单可切换语言或打开配置文件。
"""

ZH = {
    "Semiconductor Device Simulator": "半导体器件仿真器",
    "File": "文件",
    "Exit": "退出",
    "Simulation": "仿真",
    "Help": "帮助",
    "Documentation": "文档",
    "About SDS": "关于软件",
    "Close": "关闭",
    "Desktop workbench for semiconductor device simulations with sd-sim.": "基于 sd-sim 的半导体器件仿真工作台。",
    "Solver version {} ({})": "求解器版本 {}（{}）",
    "Run Current": "运行当前",
    "Run Selected": "运行选中",
    "Stop Queue": "停止队列",
    "Open Run Folder": "打开运行目录",
    "Could not open the run folder: {}": "无法打开运行目录：{}",
    "Open Saved Run": "打开已有运行",
    "Settings": "设置",
    "Open Settings File": "打开配置文件",
    "Switch to English": "切换到 English",
    "Using defaults. Check the settings file: {}": "正在使用默认设置。请检查配置文件：{}",
    "Could not save the language: {}": "无法保存显示语言：{}",
    "Could not open the settings file: {}": "无法打开配置文件：{}",
    "Could not open run for {}: {}": "无法打开 {} 的运行结果：{}",
    "Built-in cases": "内置算例",
    "Case": "算例",
    "Status": "状态",
    "PN2D Zero Bias": "PN2D 零偏",
    "PN2D Forward I–V": "PN2D 正向 I–V",
    "PN2D Reverse BV": "PN2D 反向击穿",
    "Gate: Zero Barrier": "栅极·零势垒",
    "Gate: Positive Barrier": "栅极·正势垒",
    "Gate: Work Function": "栅极·功函数",
    "Schottky Forward I–V": "肖特基正向 I–V",
    "BVMethods Full Breakdown": "BVMethods 完整击穿",
    "Curve": "二维曲线",
    "Mesh": "网格",
    "Final Field": "末态云图",
    "Log": "日志",
    "Sentaurus Script": "Sentaurus 脚本",
    "Voltage": "电压坐标",
    "Inner voltage": "内压",
    "Outer voltage": "外压",
    "Linear": "线性",
    "Log10": "对数",
    "Zoom In": "放大",
    "Zoom Out": "缩小",
    "Log scales omit non-positive values.": "对数坐标已省略零值和负值。",
    "No points can be displayed with the selected log scales.": "当前对数坐标下没有可显示的点，请切回线性坐标。",
    "Inner voltage unavailable. Select outer voltage.": "该运行缺少内压数据，请切换到外压。",
    "Inner voltage unavailable for some points.": "部分数据点缺少内压，可切换到外压查看。",
    "Rotate Left 90°": "逆时针转 90°",
    "Rotate Right 90°": "顺时针转 90°",
    "Reset Rotation": "复位方向",
    "Show Mesh": "显示网格",
    "Original Sentaurus input · read only": "原始 Sentaurus 输入 · 只读",
    "Original script unavailable. Prepare a new case catalog.":
        "当前算例库缺少原始脚本，请重新准备算例库。",
    "Field": "物理量",
    "Potential": "电势",
    "Electric field magnitude": "电场幅值",
    "Electron density": "电子浓度",
    "Hole density": "空穴浓度",
    "Net doping": "净掺杂",
    "Showing the end of the log. Open the run folder for the full log.":
        "仅显示日志末尾。完整日志见运行目录。",
    "No curve available": "当前运行尚无曲线",
    "No mesh available": "当前运行尚无网格",
    "No final field available": "当前运行没有可显示的末态场量",
    "Could not read curve: {}": "无法读取曲线：{}",
    "Could not read mesh or field: {}": "无法读取网格或场量：{}",
    "Could not read original script: {}": "无法读取原始脚本：{}",
    "Anode": "阳极",
    "Drain": "漏极",
    "Swept terminal": "扫描端子",
    "{} inner voltage (V)": "{}内电压 (V)",
    "{} outer voltage (V)": "{}外电压 (V)",
    "Anode total current (A/µm)": "阳极总电流 (A/µm)",
    "{} outer voltage {:.6g} V": "{} 外电压 {:.6g} V",
    "Swept terminal current (A/µm)": "扫描端子电流 (A/µm)",
    "Stage {} · {}": "第 {} 阶段 · {}",
    "Last accepted point": "最后接受点",
    "Final state": "末态",
    "Ready": "待运行",
    "Queued": "排队中",
    "Running": "运行中",
    "Target reached": "达到目标",
    "Current stop reached": "达到电流停止阈值",
    "Some stage targets were not reached": "求解未达到全部阶段目标",
    "Solver exited (code {})": "求解未完成（退出码 {}）",
    "Run failed: {}": "运行失败：{}",
    "Canceled in queue": "已取消排队",
    "Use Ctrl or Shift to select multiple cases": "可用 Ctrl 或 Shift 选择多个算例",
    "{}/{} cases processed": "已处理 {}/{} 个算例",
    "Running {}": "正在运行 {}",
    "Starting solver": "正在启动求解器",
    "Stage {} · last accepted bias {:.6g} V": "第 {} 阶段 · 最后接受偏压 {:.6g} V",
    "Batch finished": "批次结束",
    "Queue stopped": "队列已停止",
    "Finishing the current case before stopping the queue": "当前算例结束后停止队列",
    "Close the window after the current run finishes": "当前批次结束后可关闭窗口",
}
