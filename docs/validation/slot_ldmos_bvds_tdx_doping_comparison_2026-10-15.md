# SLOT BVDS：tdx GRD/DAT源端掺杂交叉核对

## 结果

已在VM使用Sentaurus Data Explorer U-2022.12，将原生输入`grid.tdr`及
原生BVDS零偏压结果`cold_des.tdr`分别转换为DF-ISE `.grd/.dat`。
两份GRD逐字节相同。按区域和坐标核对后，输入网格中的Nd/Na与旧Vela输入一致，
冷态输出中的Nd/Na与当前Vela对齐输入一致，差别仅为文本转换/浮点表示舍入。

这独立确认了[历史配置比较](slot_ldmos_bvds_history_comparison_2026-10-15.md)的
关键区别：原始TDR保存了大量相互补偿的施主、受主分量；本次原生冷态输出
则为净掺杂正负拆分。不是节点编号错配，也不是Vela旧导入凭空产生了大Nd/Na。
本次没有运行sdevice，没有改变原TDR、物理或任何仿真任务。

## 1. 转换命令和证据完整性

在独立输出目录执行以下命令（输入路径为示意）：

```bash
tdx -dd /path/to/grid.tdr input
tdx -dd /path/to/cold_des.tdr cold
```

两个转换命令均正常返回，生成`input.grd/input.dat`及`cold.grd/cold.dat`。
首次继承的许可证设置无法取得TDX许可；改用用户此前指定的许可证服务后成功，
仅改变本次命令环境，没有修改全局配置。

**转换日志同时报告`Conversion is incomplete!`**：不支持的接触级字段
如`ContactExternalVoltage`、`ContactCurrentFlux`、`ContactCharge`被忽略。
本次所用三个掺杂字段均为vertex数据，在Si和PolySilicon上的数量完整并通过
逐值交叉核验。因此这批文件适合掺杂比较，但不能代替原TDR进行完整接触物理复原。

| 文件 | SHA256 |
|---|---|
| 原输入TDR | `e75b9a626c771e3ca834ab92c146ed8eba892c3d2e7334b51d14319393ecaafe` |
| 原冷态TDR | `a41e3576b2cf179e1b630f47fb21992d4e02a47255290039662e31001bdcfcc4` |
| 两份GRD | `6f62b9e35d73f69830e337d246c98ae273fe12f39ef67046bf0efc6e30df687f` |
| 输入DAT | `e8155f0c5fd80f49b0c8b30b4bbc06c9edcf7e8dd76b92dffd53f8e691435d90` |
| 冷态DAT | `49a98999387d165e563c18670707bfb8c7b29c4ac4d874bdbeac24a255cfe069` |

原TDR转换前后SHA256一致；下载的四个GRD/DAT文件与VM哈希一致。
详细日志、脚本及CSV保存在忽略目录`build/s15tdx_doping/`。

## 2. 节点对应方法

- GRD为13881顶点；包含界面/接触单元的元素数28048，不能与仅统计体三角形的
  27275直接比较。
- 由GRD区域单元及边连接恢复各区域顶点，逐个解析DAT的`validity`，保留同名字段的
  不同区域实例；不能用仅按字段名索引的解析器覆盖另一材料的数据。
- 按二维坐标匹配Vela节点，确认一一对应且顺序也相同；最大坐标差
  9.93e-16 µm，坐标变换为恒等。
- Si为12410节点，PolySilicon为350节点。源接触71节点，源接触及一圈邻接的Si
  节点共144个，逐节点表为`source_neighborhood.csv`。

## 3. 源端关键节点结果

下表Nd/Na均为cm⁻³。3598位于源接触，其他所列节点为附近Si节点。

| 节点 | 原输入DAT Nd | 原输入DAT Na | 原生冷态DAT Nd | 原生冷态DAT Na |
|---|---:|---:|---:|---:|
| 3588 | 9.87344746e20 | 9.83875767e20 | 3.46897864e18 | 0 |
| 3597 | 9.21265242e20 | 9.21584133e20 | 0 | 3.18890344e17 |
| 3598 | 1.18097860e21 | 1.18033350e21 | 6.45099408e17 | 0 |
| 3599 | 9.04882732e20 | 9.06310693e20 | 0 | 1.42796090e18 |
| 3605 | 8.55735202e20 | 8.60490374e20 | 0 | 4.75517197e18 |

3598坐标为(0.004951569545, 4.426645508) µm；3599为
(0.008206610189, 4.428935547) µm。二者正是此前接触边保护涉及的位置。

## 4. 全区域核验

| 指标 | Si | PolySilicon |
|---|---:|---:|
| 输入中Nd、Na同时为正的节点 | 12410/12410 | 350/350 |
| 冷态中Nd、Na同时为正的节点 | 0/12410 | 0/350 |
| 输入/冷态净掺杂字段最大差 | 0 | 0 |
| 输入DAT与旧Vela分量最大相对差 | 6.65e-16以内 | 6.61e-16以内 |
| 冷态DAT与当前Vela分量最大相对差 | 4.99e-16以内 | 4.38e-16以内 |

在所有有效人口节点，冷态导出值均满足
`Nd=max(DopingConcentration,0)`、`Na=max(-DopingConcentration,0)`，
对导出的文本数值检验最大差为0。高浓度文本舍入的绝对差可达约7.9e5 cm⁻³，
但相对差小于7e-16，不能把它误判为掺杂模型变化。

## 5. 解释边界

tdx做的是已有文件格式转换；两份TDR在转换前就不同，不能把Nd/Na变化归因于
tdx执行了净掺杂重建。本次直接证明的是文件内容与Vela输入的对应关系。
原生究竟依据哪些字段优先级重建实际模型输入，还需结合该版本的输入规则及运行
设置解释，不能仅凭一次输出倒推所有内部实现。

此前[同态BGN/接触核对](slot_ldmos_bvds_zero_state_comparison_2026-10-14.md)
已显示，采用原生有效分量后BGN和接触势明显更一致；本批GRD/DAT证据与之相符。
若目标是匹配本次Sentaurus参考，当前有效Nd/Na有直接依据；若器件本意是保留真实
补偿杂质，则应另行明确原生掺杂输入，而不能仅在Vela侧回退旧数组以追求收敛。

## 6. U-2022.12官方规则与运行日志确认：差异源头已定位

继续查阅VM实际安装版本的《Sentaurus Device User Guide》，第2章
“Specifying Doping Species”，印刷页60–62，并核对该版本`datexcodes.txt`、
实际BVDS `run.cmd`及原生运行日志。此前第5节保留的字段解释问题得到以下直接证据。
本节仅查询和整理证据，没有启动新的sdevice仿真。

### 6.1 Grid指定输入文件，但不是把所有同名汇总数组直接作为模型输入

手册规则可概括为：

1. SDevice从`File { Grid=... }`读取净掺杂、总掺杂及独立掺杂物种。
2. 独立物种的识别由`datexcodes.txt`的`doping=donor/acceptor`声明决定；
   对同一物种，存在活性浓度时优先使用活性浓度，否则使用化学浓度。
3. 通常优先读取已有`DopingConcentration`；不存在时由识别出的独立物种重建。
   `ComputeDopingConcentration`用于强制由物种重算，不能使未识别的汇总数组
   自动变成独立物种。
4. 若没有独立物种，使用净/总掺杂重建施主、受主：
   `Nd=(Ntot+Nnet)/2`、`Na=(Ntot-Nnet)/2`。
   本次缺少总掺杂，实际运行及导出证实采用净掺杂正负拆分。

注意：PDF文本提取可能丢失数学绝对值符号。本报告不依赖该符号的文本识别，
缺总量时的非负拆分直接由本次运行日志及12760节点导出数据共同验证。

### 6.2 安装字段定义和实际脚本

安装定义中`DonorConcentration`和`AcceptorConcentration`有单位、插值及显示等
属性，但没有`doping=donor/acceptor`声明。与之相对，
`NDopantConcentration/NDopantActiveConcentration`和
`PDopantConcentration/PDopantActiveConcentration`属于已识别的掺杂物种系列。

现场没有home或运行目录的`datexcodes.txt`覆盖，查询环境`STROOT_LIB`未设置。
更关键的是，原运行日志自己列出的识别物种与上述安装定义一致。

脚本实际使用`Grid="grid.tdr"`，没有启用`ComputeDopingConcentration`，
也没有`Add2TotalDoping`陷阱设置。`Plot`中的
`Doping DonorConcentration AcceptorConcentration`只控制结果输出，
不是强制按输入文件中同名Nd/Na数组初始化的指令。

### 6.3 本次原生日志直接确认的分支

原运行日志第2506行指出未找到独立掺杂物种；随后明确记录：

- 施主/受主浓度：`computed based on net doping`；
- 净掺杂：从掺杂文件取得；
- 总掺杂：文件中未找到。

因此，因果链已明确：输入保存了Nd/Na汇总字段，但未提供本次识别的独立物种，
也未提供`TotalConcentration` → SDevice按已有净掺杂初始化Nd/Na →
冷态保存这些重建数组 → 得到本报告前述差异。
它发生在掺杂初始化/解释阶段，不是Poisson或漂移扩散迭代使固定杂质发生了
扩散、迁移或相互湮灭，也不是tdx转换造成。

对3598节点，原生读入`Nnet=6.45099407914e17 cm^-3`，
于是使用`Nd=Nnet, Na=0`；旧输入同时约1.18e21的补偿分量没有按独立物种参与初始化。

### 6.4 对后续物理选择的意义

“当前Vela有效分量匹配本次原生运行”有充分证据，但不等于已经证明它符合
器件工艺预期。若真实意图是保留补偿杂质，需在独立原生配置中明确输入：

- 合适的已识别施主/受主物种及其活性浓度；或
- 在无独立物种路径提供一致的`DopingConcentration`和
  `TotalConcentration=Nd+Na`，再检查实际重建结果。

选择真实元素或通用物种须与材料、活化和电离模型相符，不能随意把汇总Nd改名为
某一种元素。优先采用独立输入配置，避免修改系统全局datexcodes。
先做初始化输出核验及同态BGN/迁移率检查，再决定生产比较口径；本次未实施改变。

手册摘录、字段定义、脚本、运行日志片段及源文件SHA256存于
`build/s15tdx_doping/manual_audit/`。官方手册正文及机器路径不纳入版本控制。
