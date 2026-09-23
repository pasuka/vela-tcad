# T470p：公开二进制状态格式与序列化候选

日期：2026-09-22。分支：`codex/templates-ldmos-phase-a`。

## 任务范围与状态

按用户授权比较 HDF5＋HighFive、FlatBuffers，以及未压缩 NPY；VDS1 为
连续数组参照。实现、本地/远端检查及 T470p 三组配对全部完成，12 曲线/96 精确点通过；
另有 48 组跨格式精确点的全部保存字段逐位相同。远端仿真进程正常退出。
生产默认仍为 CSV，新增构建选项 `VELA_ENABLE_STATE_FORMATS` 默认关闭。
本轮不更改 Newton、物理模型、收敛门限、预测规则或保存频率。

本轮只验证 DD 三方程状态，不包含 D0 晶格温度状态。小测试中的量子场只是
序列化字段覆盖，不能解释为新增量子模型的曲线资格。

## 文件合同与实现

| 候选 | 后缀 | 保存方式 | C++ / Python |
|---|---|---|---|
| VDS1 | `.vds` | 定长小端头部＋按字段连续的 float64 数组 | 现有轻量编解码器 |
| HDF5 | `.h5` | 命名一维 float64 dataset；连续、未压缩 | HighFive / h5py |
| FlatBuffers | `.vfb` | 有版本的 DDState table、字段名、连续 float64 vector | flatc 生成代码 / FlatBuffers runtime |
| NPY | `.npy` | 一维结构化数组，每个节点为一条记录，字段为 `<f8` | 有限格式读写器 / NumPy |

共同必需字段为 `psi, phin, phip, electrons_m3, holes_m3`；可选量子字段为
`electron_quantum_potential_V, electron_quantum_potential_like_V`，后者不能单独存在。
四个 QF 分裂字段必须同时出现：电子/空穴各自的参考值和增量。它们独立保存，
不能由已舍入的 `phin/phip` 反推。节点 ID 隐含为网格的原顺序 `0..N-1`。
电势为 V，密度为 m⁻³；内部浓度单位在公共入口转换。

HDF5 根属性为 `schema=vela.ddstate/1`、`units=V;m^-3`、`node_count`；
`field_names` 固定字段次序，各字段 dataset 长度必须等于网格节点数。
FlatBuffers 文件标识为 `VDSF`，根表版本 1、单位及字段顺序均验证；C++ 使用
生成的 verifier 检查结构，再检查节点数和物理字段。其公开合同见
[dd_state.fbs](../../schemas/dd_state.fbs)。

NPY 接受 1.0/2.0 头部、一维 C-order 的命名 `<f8` 记录，字段顺序按本合同，
不接受对象、pickle、任意 dtype 或任意形状；本项目写出 NPY 1.0。
NPY 本身不携带通用单位、网格或模型元数据，本合同规定 SI 单位，运行计划继续
冻结网格、配置及来源散列。因此它是公开格式中的受限状态 schema，不是任意
NumPy 文件都可以作为重启文件。

所有候选保留既有正常 binary64 数值与 QF 分裂精度；沿用 CSV/VDS1 的显式
subnormal→0 写出政策，不宣称对所有浮点比特无损。读取拒绝非有限字段、错误
尺寸和不一致的 QF 坐标。文件格式校验不能替代外层网格与物理配置来源审计。

实现入口：

- [C++ 格式实现](../../src/io/DDStateFormats.cpp) 与
  [共享状态入口](../../src/io/DDSolutionBinary.cpp)。
- [Python 格式实现](../../scripts/dd_state_public.py) 与
  [D5 适配器](../../scripts/ldmos_binary_adapter.py)。
- [跨语言及失败恢复测试](../../tests/regression/test_public_state_worker.py)。

第一版公共格式通过内存中的 VDS1 规范缓冲区复用单位/字段检查；公开文件内保存
实际 typed arrays，并非 VDS1 blob。该桥接包含额外复制和遍历，当前性能结果评价
这一版完整实现，不能解释成 HDF5、FlatBuffers 或 NPY 的理论吞吐上限；也不能称
为端到端零拷贝。

## 构建、环境和检查

在已有 UCRT64 Release preset 上显式启用：

```powershell
cmake --preset windows-ucrt64-release -DVELA_ENABLE_STATE_FORMATS=ON
cmake --build --preset windows-ucrt64-release --target vela_example_runner test_dc_sweep --parallel 4
ctest --test-dir build-release -R 'public_state|binary_state|dc_worker_protocol' --output-on-failure
./build-release/test_dc_sweep.exe '[dc_sweep]'
```

需实际检测到 HDF5、HighFive 头文件和 flatc/FlatBuffers 头文件。C++ binding 在
构建目录生成；Python binding 入库，并由 `flatc --python -o scripts/generated
schemas/dd_state.fbs` 再生。开启该可选功能后的跨语言测试另需 NumPy、h5py 和
FlatBuffers Python；默认关闭构建不新增这三项 Python 要求。

本机及 T470p 均通过：4 个 CTest 协议/回归条目；104 个 DC 扫描测试、3817 条断言。
另外 32 项驱动/缓存/状态 Python 回归通过；在独立构建目录将新选项关闭后，
共享状态入口对象编译通过，未向关闭功能的路径引入 HighFive/FlatBuffers 符号。
其中两种单位制、CSV 对照及 QF/量子字段覆盖的 restart 测试为 2 例、114 条断言。
跨语言测试覆盖 9 种输入/输出格式组合、C++ 输出再作为输入、节点数不符、截断、
NaN、错误单位、float32 错误 dtype，并确认错误请求不污染后续正常请求。

本机与 T470p 的 HighFive 均为 3.3.0，并检测 NumPy 2.5.3、h5py 3.16.0、flatc/FlatBuffers Python
25.12.19。h5py 报告编译 HDF5 2.1.1、运行时 2.2.0 的版本警告；保留日志，
未屏蔽。实际互读测试通过只说明已测试路径可用，不扩展为其他 HDF5 功能兼容性证明。

## T470p 配对协议

独立目录 `D:/code-repo/vela-bench/ldmos_public_formats_20260922`，工作站证据目录
`build/t470p_public_formats_20260922`。源码清单冻结 1098 个文件，远端解包后逐项
SHA256 校验；远端同一 Release 程序用于所有控制与候选。

保持既有 D5 的 10241 节点网格、物性、参考曲线、零压种子、0.2 V 最大推进及现有预测规则。
使用 UMFPACK 1 线程、BLAS 1 线程，保留线性分析和装配结构跨点复用；Fermi
端点缓存关闭。每种候选按 `VDS1/Vg4 → 候选/Vg4 → 候选/Vg8 → VDS1/Vg8`
串行运行，共三组、12 条短曲线、96 精确点；每曲线为从 Vd=0 开始的前 8 个
参考点，终点 9.333333333 V。

Python 依赖导入在配对计时开始前完成。端到端墙钟包括父驱动、子求解器、状态
读写、状态比较、冷读及来源审计；原二进制/CSV 缓存策略相同。逐曲线核验原
块/逐行/电流闭合门限、精确点状态、推进轨迹、所有记录的计数及来源散列。
不将一次配对解释成重复稳定性，不沿用旧程序的墙钟作为本轮基线。

再现时在冻结包与基线路径上分别使用：

```text
python scripts/run_templates_ldmos_d5_cost_controls.py --package <frozen-package> --runner <same-release-exe> --baseline <qualified-baseline> --output <new-directory> --points 8 --control binary --candidate hdf5
```

另两组将 `hdf5` 改成 `flatbuffers`、`npy`。完整命令和实际进度写入远端
`status.json`；源文件与运行库散列在 `source_manifest.json` 及各组计划中。

## 配对结果

三组均正常结束，原门限、来源审计及全部尝试轨迹通过，所有已记录性能计数相同。
各格式 Vg4/Vg8 均为 367/349 次 Newton 更新、68/68 次点服务；无外部减步回退。
这不表示没有内部恢复，内部恢复也保持原计数。本轮是从零压种子开始的漏压扫描，
不包含栅压预偏置。

| 候选 | Vg/V | 同轮 VDS1/s | 候选/s | 候选耗时变化 |
|---|---:|---:|---:|---:|
| HDF5＋HighFive | 4 | 141.7455 | 141.7624 | +0.012% |
| HDF5＋HighFive | 8 | 140.7249 | 141.5477 | +0.585% |
| FlatBuffers | 4 | 142.1680 | 143.4034 | +0.869% |
| FlatBuffers | 8 | 141.3457 | 140.9726 | −0.264% |
| NPY | 4 | 142.0308 | 142.4606 | +0.303% |
| NPY | 8 | 141.8551 | 143.3338 | +1.042% |

双栅压合计相对各自配对 VDS1 的耗时增加分别为 **0.297%、0.304%、0.672%**。
未观察到可认领的端到端加速，单轮也不能建立稳定的格式速度排名。同程序三组
VDS1 控制本身在 Vg4 为 141.75–142.17 s、Vg8 为 140.72–141.86 s，提示应保留
对细小差异的判断余地；这些控制不是三种候选各自的三轮重复。

对旧冻结基线，精确点 ψ/fn/fp 和 n/p 的最大差均为 0。补充独立冷读审计将候选与
同组 VDS1 的全部字段重新编码为共同规范字节：48 组精确点全部相同，包括量子
辅助字段和 QF 的参考值、增量。原程序输出文件格式和文件 SHA256 不同是预期行为。

### 状态文件体积及开销

每条曲线保持 128 个状态文件（68 个求解输出＋60 个预测种子），未减少落盘频率。
同格式两栅压的文件总字节数相同。下表只统计状态，不含日志、图和程序：

| 格式 | 每曲线字节数 | 相对 VDS1 |
|---|---:|---:|
| VDS1 | 104,870,912 | 基准 |
| HDF5 | 106,199,040 | +1.266% |
| FlatBuffers | 104,908,800 | +0.036% |
| NPY | 104,916,992 | +0.044% |

`dc.record_point` 包含状态写入等记录工作。HDF5 对照为
5.7102→5.7657 / 5.7154→5.8272 s；FlatBuffers 为
5.6775→5.5011 / 5.7127→5.5737 s；NPY 为
5.6651→5.6348 / 5.7645→5.7071 s。数值求解约 107–109 s，仍占主体。
FlatBuffers 的记录阶段节省没有转化成可靠总收益；父驱动的格式解码、行字典和审计
仍有成本。当前数据不足以把少量时间差全部归因于文件格式。

所有曲线记录的实际 BLAS 线程均为 1。子求解器峰值工作集为 241.3–249.4 MB
（十进制 MB），该口径不含 Python 父驱动，不能作为整个批次总内存。

## 结论与后续边界

- **HDF5＋HighFive 保留为正式存储/恢复的优先候选**：本轮同样正确，性能接近
  轻量 VDS1，而数据集名称、单位和版本可直接检查。尚不切换生产默认。
- **NPY 保留为 Python 数组交换与性能参照**：格式简单，但节点/模型来源仍须
  外部合同保证；本轮无额外加速证据。
- **FlatBuffers 保留已验证接口，不以加速理由继续扩大接入**：本轮结构化
  序列化通过，但行字典转换和内存桥接仍在，未兑现端到端零拷贝收益。

本轮完成候选筛选，不是 0–40 V 完整曲线、跨参考系长轨迹、重复稳定性、D0 温度
状态或任意平台资格。若将 HDF5 纳入推荐配置，应先补完整曲线与重复计时；原
CSV 默认、原物理模型及门限保留。本轮不继续扩大格式候选或改造内存传输。

## 证据

源码归档 SHA256：
`af52c0b89fd89347de2a74feff636ab63d10218dbef32ffef005833d6f3b1244`。
全部三组实际 Release runner SHA256：
`85fefd25b7b731134b0bb6f634cec49fc0673d170e1a0dfe29653304d48d3ed9`。

下载归档 `build/t470p_public_formats_20260922/public_formats_evidence.tar.gz`
SHA256 为 `853b2db26c579e5e2599fa7b8678762c8f821ee083089e0c509b1185112f8182`，
归档及其中逐文件散列已在本机核验。解包目录 `completed/` 包含各组 `summary.json`、
配对分析、完整 ledger、状态比较、计划、审计、测试/构建日志、`size_summary.json`
及 `all_field_audit.json`。完整逐点状态留在同名远端目录。
当前运行/构建/测试源文件与冻结源码清单一致；后续变动仅为本报告及文档索引。
