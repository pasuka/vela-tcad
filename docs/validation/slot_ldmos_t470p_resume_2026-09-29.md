# SLOT T470p 暂停后续算（2026-09-29）

用户授权检查 T470p 连接，并在可用时继续昨晚暂停的仿真。未修改物理模型、
求解器源码或资格门槛。本轮仅在 T470p 运行一项 Release 仿真。

## 连接与启动

- 主机 `DESKTOP-V2AJOBM` / `desktop-v2ajobm.local`，当前 IPv4 为
  `192.168.0.3`。原 `192.168.0.2` 已成为本机地址。
- 新地址的 ED25519 公钥与本机此前保存的 `192.168.0.4` 记录一致，使用
  `HostKeyAlias=192.168.0.4` 和严格主机公钥验证；未覆盖其他主机的记录。
- 远程工作目录与本机相同：
  `D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next`。
- 初次 SSH 内 Start-Process 未能保持后台进程，尚未生成运行目录或开始求解。
  改用无周期触发的任务 `VelaSlotResume20260929`，admin S4U / Limited，
  隐藏窗口；计划任务最大运行时间 5 小时，控制器预算仍为 4 小时。
- 08:22 后已确认 SSH 断开再连接时任务仍为 Running，Python 控制器 PID
  8364、Release 求解器 PID 6968 均有 CPU 活动。

## 恢复证据

源为 `build/s28homebv/run/checkpoint.json`，用户暂停记录为
`build/s28homebv/user_pause.json`。该源扫描未完成，最后的 running 尝试
属于主动中断，未作为合格点导入。

最后三个接受状态在本机及 T470p 分别通过现有 `SolverBackend.validate`
与 `validate_confirmation`：验证文件哈希、完整器件 Newton 证据、全部有效
载流子行、全局闭合、负载线及确认态绑定。未执行昨晚所有 163 个新增点的
全程事后审计，不以这三个父态的复验代替全程审计。

最后合格父态：

| 量 | 数值 |
|---|---:|
| 物理内压 | 38.45175602808871 V |
| 物理外压 | 38.56477228340672 V |
| 漏极电流 | 1.1301625531801068e-13 A/µm |
| 参考原点 | 36 V |
| 恢复外压步长 | 0.0018906678529976304 V |

迁移包 `build/s29t470bv/transfer.tar.gz` 为 45,022,130 字节，含 81 个文件。
安装前逐文件核对包内哈希及远程已有文件；已有文件仅在完全相同时复用，
未覆盖远程不同内容。远程 `remote_preflight.json` 记录完整预检通过。

## 运行入口与范围

- 新运行目录：`build/s29t470bv/run`。
- 启动脚本：`build/s29t470bv/run.ps1`。
- 父态：`build/s29t470bv/seed_window.json`；物理配置 `base.json`。
- 状态读取：远程 `build/s29t470bv/status.py`，允许 FILE_SHARE_DELETE。
- 目标：外压 200000 V，或首个完整合格态达到 `abs(Id)>=1e-7 A/µm`。
- 沿用 limited_rational、自适应增长 1.4、外压上限 20000 V、接受态内压
  增量上限 0.1 V、试探内压增量上限 0.05 V、电流增长比例上限 0.25。
- 候选预算 600 秒、子求解超时 300 秒、最多 500 次尝试、总预算 14400 秒。
- OMP / OPENBLAS / MKL 均为单线程。模型配置及冻结 Release runner 与昨晚一致。

08:22 状态：3 个导入父态、0 个新增接受点，第 1 次耦合外电阻尝试运行中，
首个预测模型为 limited_endpoint。新的段只继承最近三个合格端口状态；不声称
继承昨晚全部预测器评分历史。本次沿用已保存的下一步步长，没有重置为 0.025 V。
未达到击穿电流判据，尚未完成本轮续算。原 20 分钟 Codex 监测保持暂停。

既有 h5py/HDF5 版本警告仍存在；实际预检已通过，不能把该 stderr 警告单独
解释为启动失败。结束时应检查 `exit.json`、最终检查点及独立审计后再判定完成。
