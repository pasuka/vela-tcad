# SLOT 无限幅 Bank–Rose 对照

2026-10-11。按用户要求关闭原merit回溯线搜索，启用已有Bank–Rose Algorithm Global实现；仍关闭两种QF分量限幅。

本批完成2组同父态10 mV试算，新策略收敛点0个，总墙钟297.2秒。

## 固定条件与参数

同一5 mV合格父态，Vg=5 V、Vd目标10 mV，300 K、原点0 V；13881节点/27275三角形；冻结oneMKL PARDISO binary64后端、1线程。
保留AverageBox/Fermi/effective_band_sg、全目标物理、原正则化与闭合替换、用户选择的整体残差或原始更新误差停止条件，以及全局闭合和边界保护；行比继续报告。
`line_search=false`，`bank_rose.enabled=true`，`not_damped=0`，`damping_factor=1`，`stable_merit_comparison=false`；两种QF限幅=0，因此不触发截断配套Poisson重修正。
NotDamped=0使Bank–Rose从第1步启用；25会让前25步先走完整步，不能检验当前首步阻尼问题。
Bank–Rose仍是Newton方向加全局阻尼，并非取消所有试探：alpha=1/(1+K*rho)，拒绝时K乘10（K=0先取seed=1），接受后K乘0.1，omega=1e-4，每步最多32次试探。
两组最小alpha均为1e-14（显式低于原默认1e-8），为大K初值提供可执行范围；这不是改变收敛门槛。

A组initial_K=0；B组initial_K=1762853859.37，由同态已保存原始方向与rho0=1.14682716398计算：
alpha0=0.05/max(|原始QF方向[V]|)，K0=(1/alpha0-1)/rho0。
B组预设首步alpha=4.94635995858e-10，使首步最大QF位移约0.05 V；之后按Bank–Rose递推K，不逐分量裁剪，也不逐步强制该0.05 V上限。
每组最多40更新、软1080秒/硬1200秒，串行总上限2500秒。未修改生产默认或C++。

## 结果

| 组 | 更新数 | 新策略收敛 | 返回码 | 耗时 s | 最后rho | 最后原始E |
|---|---:|---|---:|---:|---:|---:|
| A_K0 | 0 | False | 1 | 32.8 | — | — |
| B_initial_damped | 6 | False | 1 | 264.2 | 1.1468259097603868 | 6391.326708371905 |

### A_K0

求解器原因字段：`None`；stderr：`vela_example_runner failed: SG source requires finite positive populations`。


### B_initial_damped

求解器原因字段：`None`；stderr：`vela_example_runner failed: SG source requires finite positive populations`。

末个接受试探alpha=4.94636048e-07，电子/空穴闭合比=0.999999987/0.0426269696，坏行=22/1987，最大行比=1.71039021。
实际有限更新误差为0.00315906543，与原始方向E分开保存，前者不参与成功判定。

## 核对与解释

本次A组在首次试探退出；B组接受6步后在后续试探同样触发SG人口检查。B组rho由1.146827163975降至1.146825909760，相对改善约1.094e-6（0.0001094%）；电子/空穴闭合比约1/0.04263，未通过原门槛。没有把极小有限更新误差当作收敛。
结果说明较大K能绕过第一个完整步立即出域，但本次参数下没有实质性建立输运平衡；不能据两组失败否定Bank–Rose算法本身。下一步宜先为可识别的无效人口试探增加有界拒绝/增K重试，并检查实际人口可行性；保留所有最终验收，不采用捕获任意异常后静默继续的做法。
两组首次物理矩阵、正则化矩阵、闭合系统和raw方向均与前一无限幅10 mV系统逐字节相同；所有已接受步的raw与processed向量逐分量相同，未恢复隐藏QF限幅。
当前Bank–Rose对有限候选调用残差函数时，没有将SG人口有效域异常转换为“拒绝、增大K、重试”。因此K=0首次完整步的源项异常会直接传播退出，不能说已经穷尽Bank–Rose的所有阻尼尝试。
较大K初值用于区分“首个完整步立即出域”与“极小统一步能否推进”。极小实际更新不等于接近解，必须同时查看原始E、残差和闭合。此实现依据公开Algorithm Global，不宣称完全复刻Sentaurus内部算法。
停止逻辑10个边界/非有限/闭合用例、零更新和物理单位RMS检查通过；输入、二进制、运行库、执行脚本哈希不变；配置审计通过。未重编译或改动求解器；既有h5py/HDF5版本警告存在但文件核对完成。
证据目录：`build/s11row_bankrose/`，包含contract、summary、audit及两组config/stderr/exit、逐步Bank–Rose K与alpha历史、矩阵和状态。字段`independent_accepted`为复用控制器的容器名，本批仅表示不同K诊断臂，不是独立路径复现或生产续算。
