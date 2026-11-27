# SLOT float80 类型替换试验与可选混合精度后端

日期：2026-10-11。范围：用户要求尝试以 float80 替换默认 float64。
本批完成隔离类型替换编译检查、可用的扩展精度线性后端及同系统对照；
**未完成全链路 float80，未改变默认后端、物理模型或验收门槛，新增合格点 0。**

## 1. 可用实现与全局替换障碍

新增显式环境选项 `VELA_LINEAR_SOLVER=sparselu_long_double`：

1. 将已经装配的 binary64 矩阵和 RHS 精确提升为 `long double`；
2. 以 Eigen SparseLU/COLAMD 做扩展精度分解和回代；
3. 将方向转回现有 binary64 公共接口。

不复用该后端的分解缓存；平台的 `long double` 不比 `double` 更宽时明确报错。
UMFPACK/SparseLU 原默认选择顺序不变。状态、物理残差和 Jacobian 装配仍是 double。
尤其生产闭合行的 Woodbury 小矩阵及修正组合仍用 double：本选项不等于整个
稳定化线性系统的所有运算都升精度。

本机 MSYS2 UCRT64 GCC 16.2.0 实测：

|量|double|long double|
|---|---:|---:|
|有效二进制尾数位数|53|64|
|epsilon|2.220446049250313e-16|1.0842021724855044e-19|
|对象存储字节|8|16|

因此是通常称为 x87 80-bit 的扩展精度，带填充占 16 字节，不是 128-bit 精度。
数值数组存储加倍；整个稀疏求解器内存还有索引、工作区与填充，不能直接称总内存恰好两倍。

在 `build/s11float80/overlay` 同时将 Real、向量、稀疏矩阵和点坐标改为
long double，仅用于 `-fsyntax-only` 编译，没有覆盖生产头文件。
实际失败包括 UMFPACK 的 double 指针接口、CHOLMOD 缺失 long-double 映射；
Newton 编译也在 `ElementLombardi.h:38` 的 `std::max(0.,distance)` 出现
double/Real 模板类型冲突。日志见 `global_type_trial.log` 与
`global_newton_trial.log`。这只是代表性检查，不是全库迁移缺口的穷尽清单。
全链路还需逐项迁移物理核、标量常量、AD/导数、状态和 HDF5/JSON 接口，禁止
靠隐式降回 double 让编译通过后宣称已完成升精度。

## 2. 保存的同一矩阵直接求解

输入来自 `build/s10backend_newton/sparselu/data/` 的 stabilized.mtx/vectors.csv；
41,643 阶、739,295 个存储项，内部无量纲变量，SLOT Vg=5 V/Vd=5 mV 的
旧冻结配置。先按 double 解析文件，再提升，保证两臂输入相同；不是把十进制
文件分别按两个精度解析。矩阵 SHA-256：
`f8947dd2a7d063ade4fbba28e3ac5b4ffe40dea9ae3f0a8ab76c2f5c155bb94c`。

独立 C++ 直接分解完整稳定化矩阵，没有采用生产 Woodbury 路径。
其 Boost 113-bit 残差评估与 60 位十进制独立复算一致：

|方向保存方式|max abs(Ax-b)|相对 double 的改善|
|---|---:|---:|
|double SparseLU|1.213175473516464e-12|1|
|long double SparseLU，保留扩展方向|1.267023647058956e-15|约958倍|
|同一扩展方向转回 double|4.519674807822645e-14|约27倍|

以上是**内部线性方程残差**，不是器件电流、非线性残差或载流子行比。
分量后向误差定义为 `abs(Ax-b)_i / (abs(A)*abs(x)+abs(b))_i`，最大值两臂均为1。
最大值所在24156行只有对角1.00000001、RHS=0；方向分别约
-2.3303e-15/+1.6960e-17，只要非零该比例就是1，不能由此推断物理弱行失败。
该指标P99从0.2020降至0.005496；仍需结合绝对量和实际物理资格。

两臂方向 L2 相对差约9.342e-4（统计转为double后的方向），也说明提高线性精度
不只改变最后一位。分解耗时首次观察103.62/105.60秒，曾与编译等操作有重叠，
不作为受控性能基准，也不能推广为所有工况的float80开销。

独立复算发现此环境 NumPy 的 longdouble 字符串读入会先丢到double；分析脚本
因此用十进制及64位二进制有效位舍入恢复 C++ max_digits10 输出，再复算，
没有依赖该读入路径来证明float80收益。

## 3. 接入现有 SLOT 生产方向路径

两臂使用同一份本批冻结runner、同一 `v5_qf.h5`，仅后端不同。
配置继承 `build/s10initial05/v5_qf_matrix/config.json`，AverageBox、Fermi、
既有迁移率/电离配置及全部门槛不变，参考原点0 V；采用本批当前编译常数。
本节不能与上一节旧常数runner的绝对方向直接归因为精度差异。

运行 `newton_linear_system_probe`，每臂180秒上限，两臂400秒总兜底；
实际24.79/27.44秒，均退出0。未形成或接受Newton试探。
物理/正则化/稳定化矩阵各自两臂SHA相同，RHS和残差逐值相同。
原始方向范数为699607.14445/699587.06358，方向L2相对差7.039746e-4，
最大内单位差190.32868（索引38216，即空穴节点10454）。
这些差异还包含double Woodbury修正的响应，不能当作原始矩阵精确解误差。

两臂摘要的 `converged=true` 仅指该线性探针成功；`qualified=false`、
`accepted_updates=0`。本批没有进行完整非线性复求，不能说器件收敛已改善或已失败。

## 4. 验证与下一步

Release构建通过：vela_example_runner、test_linear_solver。
线性测试通过14例/348条断言，覆盖新后端非对称跨尺度系统、变化RHS、变化系数、
奇异失败及平台精度检查；实际 SLOT 两臂也成功。
重新链接后的完整Newton测试通过176例/4649条断言（默认后端）。
显式设置long-double后端后的Bank–Rose测试通过3例/29条断言。

证据在 `build/s11float80/`：源程序compare.cpp、run.py、analyze.py，
线性结果linear_comparison.json、direction_analysis.json、verification.json，
运行身份contract.json、生产探针summary.json及各臂矩阵/向量。

后续优先保持配置和更新策略不变，进行少量重线性化的double/long-double后端
对照，比较实际行比、闭合、限幅及耗时；如改善只停留在线性残差，不将其作为
改变生产默认的依据。若需完整float80，应单独实施全链路标量迁移与回归，而非
继续把本混合精度后端称为全局替换。

资料：[Eigen自定义标量](https://libeigen.gitlab.io/eigen/docs-master/TopicCustomizing_CustomScalar.html)、
[GCC扩展浮点类型](https://gcc.gnu.org/onlinedocs/gcc/Floating-Types.html)。
