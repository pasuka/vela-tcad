# Sentaurus VMware SSH Workflow

本文记录在 Windows 家用电脑上通过 VMware Workstation Pro 运行
CentOS/Sentaurus 2018 虚拟机，并用 SSH 自动传输、运行、取回仿真结果的
操作流程。

当前已验证的连接信息：

- Windows 仓库路径：
  `D:\code-repo\vela-tcad`
- 本地 Sentaurus PN2D 脚本和结果目录：
  `D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source`
- 虚拟机 SSH 别名：`sentaurus`
- 虚拟机用户：`tcad`
- 当前虚拟机 IP：`192.168.119.128`
- 当前虚拟机登录后目录：`/home/tcad`

## 目标

使用 VMware Host-only 网络让宿主机可以 SSH 登录虚拟机，同时不允许虚拟机
通过默认网关访问互联网。之后通过 `ssh` 和 `scp` 完成：

1. 从 Windows 上传 Sentaurus 脚本到 CentOS。
2. 在 CentOS 中创建仿真工作目录。
3. 远程运行 `sde -e -l <sde-script>` 和 `sdevice <sdevice-script>`。
4. 将 `.tdr`、`.plt`、`.log`、`.grd`、`.dat` 等结果拷回 Windows。

## VMware 网络设置

在虚拟机关机状态下打开 `Virtual Machine Settings -> Network Adapter`：

- 勾选 `Connected`
- 勾选 `Connect at power on`
- 选择 `Host-only: A private network shared with the host`
- 不选择 `NAT`
- 不选择 `Bridged`

在 VMware Virtual Network Editor 中确认对应 Host-only 网络：

- 类型为 Host-only，通常是 `VMnet1`
- 启用 `Connect a host virtual adapter to this network`
- 不给该网络配置 NAT

在 CentOS 中检查路由：

```bash
ip route
```

当前已验证输出类似：

```text
192.168.119.0/24 dev eth0  proto kernel  scope link  src 192.168.119.128
```

关键点是没有 `default via ...`。这表示虚拟机可以访问同一 Host-only 网段的
Windows 宿主机，但没有默认网关访问互联网。

可用以下命令验证：

```bash
ping 192.168.119.1
ping 8.8.8.8
```

期望：能 ping 通宿主机的 VMnet1 地址；不能 ping 通公网地址。

如果 DHCP 导致虚拟机 IP 改变，需要更新 Windows 的 SSH config 中
`HostName` 字段。

## CentOS SSH 服务

这台 Sentaurus 虚拟机使用较老的 CentOS/OpenSSH，不能使用 `systemctl`。
使用 SysV 服务命令：

```bash
su -
service sshd start
chkconfig sshd on
service sshd status
```

已验证状态：

```text
openssh-daemon (pid ...) is running...
```

如果 SSH 端口没有监听，可检查：

```bash
netstat -tlnp | grep :22
```

## Windows SSH 配置

旧 CentOS 的 SSH server 只提供旧主机密钥算法，例如 `ssh-rsa`。Windows 新版
OpenSSH 需要显式为这个 Host 放开兼容算法。

Windows 配置文件：

```text
C:\Users\qzw\.ssh\config
```

推荐内容：

```sshconfig
Host sentaurus
    HostName 192.168.119.128
    User tcad
    IdentityFile C:/Users/qzw/.ssh/id_rsa_sentaurus
    IdentitiesOnly yes
    HostKeyAlgorithms +ssh-rsa
    PubkeyAcceptedAlgorithms +ssh-rsa
```

注意：该文件应保存为 ASCII 或 UTF-8 无 BOM。若 OpenSSH 报
`no argument after keyword "\377\376"`，说明文件被保存成了 UTF-16。
可用 PowerShell 重新写入：

```powershell
@'
Host sentaurus
    HostName 192.168.119.128
    User tcad
    IdentityFile C:/Users/qzw/.ssh/id_rsa_sentaurus
    IdentitiesOnly yes
    HostKeyAlgorithms +ssh-rsa
    PubkeyAcceptedAlgorithms +ssh-rsa
'@ | Set-Content -Path "$env:USERPROFILE\.ssh\config" -Encoding ascii
```

如果当前 Windows OpenSSH 不认识 `PubkeyAcceptedAlgorithms`，可改用旧名称：

```sshconfig
    PubkeyAcceptedKeyTypes +ssh-rsa
```

## SSH 密钥免密登录

旧 CentOS 通常不支持 Ed25519，建议使用 RSA key：

```powershell
ssh-keygen -t rsa -b 4096 -f $env:USERPROFILE\.ssh\id_rsa_sentaurus
```

在 CentOS 的 `tcad` 用户家目录中准备权限：

```bash
mkdir -p /home/tcad/.ssh
touch /home/tcad/.ssh/authorized_keys
chown -R tcad:tcad /home/tcad/.ssh
chmod go-w /home/tcad
chmod 700 /home/tcad/.ssh
chmod 600 /home/tcad/.ssh/authorized_keys
```

从 Windows 追加公钥：

```powershell
type $env:USERPROFILE\.ssh\id_rsa_sentaurus.pub | ssh sentaurus "cat >> ~/.ssh/authorized_keys"
```

如果出现：

```text
bash: /home/tcad/.ssh/authorized_keys: Permission denied
```

说明 `.ssh` 或 `authorized_keys` 可能由 root 拥有。回到 CentOS root 修复：

```bash
su -
mkdir -p /home/tcad/.ssh
touch /home/tcad/.ssh/authorized_keys
chown -R tcad:tcad /home/tcad/.ssh
chmod go-w /home/tcad
chmod 700 /home/tcad/.ssh
chmod 600 /home/tcad/.ssh/authorized_keys
```

## 已验证的连通性测试

Windows PowerShell：

```powershell
ssh sentaurus "hostname; pwd; whoami"
```

已验证输出：

```text
sentaurus
/home/tcad
tcad
```

验证写入虚拟机：

```powershell
"hello from windows" | ssh sentaurus "cat > ~/ssh_test.txt && cat ~/ssh_test.txt"
```

已验证输出：

```text
hello from windows
```

验证从虚拟机取回文件：

```powershell
scp sentaurus:~/ssh_test.txt $env:TEMP\ssh_test.txt
type $env:TEMP\ssh_test.txt
```

已验证输出：

```text
hello from windows
```

## Sentaurus 命令检查

先确认远程 shell 能找到 Sentaurus 命令：

```powershell
ssh sentaurus "command -v sde; command -v sdevice"
```

如果没有输出，需要先在 CentOS 中加载 Sentaurus 2018 环境。常见做法是在
`~/.bashrc` 或专门的启动脚本中加入站点安装提供的环境初始化命令，例如：

```bash
source /path/to/sentaurus/setup.sh
```

具体路径取决于本机 Sentaurus 安装位置。确认后再重新运行：

```powershell
ssh sentaurus "command -v sde; command -v sdevice"
```

## 上传 PN2D 脚本

Windows 当前源目录：

```text
D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source
```

建议在虚拟机中使用独立运行目录，避免覆盖手工实验：

```powershell
ssh sentaurus "mkdir -p ~/sentaurus_runs/pn2d_sentaurus2018/source"
scp -r "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\*" sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/
```

检查远程文件：

```powershell
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && ls -1 pn2d_sde.cmd pn2d_*_sdevice.cmd models.par"
```

## 运行 Sentaurus

常用命令模式：

```bash
sde -e -l <sde-script>
sdevice <sdevice-script>
```

PN2D 当前脚本示例：

```powershell
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && sde -e -l pn2d_sde.cmd"
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && sdevice pn2d_0v_sdevice.cmd"
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && sdevice pn2d_iv_sdevice.cmd"
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && sdevice pn2d_bv_sdevice.cmd"
```

如果希望保留运行日志，可用：

```powershell
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && sdevice pn2d_iv_sdevice.cmd > run_pn2d_iv.out 2>&1"
```

查看日志尾部：

```powershell
ssh sentaurus "cd ~/sentaurus_runs/pn2d_sentaurus2018/source && tail -n 80 run_pn2d_iv.out"
```

## 取回结果

将远程 `source` 目录中的结果同步回当前 Windows 源目录：

```powershell
scp sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/*.tdr "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\"
scp sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/*.plt "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\"
scp sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/*.log "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\"
scp sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/*.grd "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\"
scp sentaurus:~/sentaurus_runs/pn2d_sentaurus2018/source/*.dat "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source\"
```

如果某类文件不存在，`scp` 会对该通配符报错；这不一定表示仿真失败。
以 Sentaurus 日志和期望输出文件为准。

## 一次性手动流程模板

从 Windows PowerShell 执行：

```powershell
$local = "D:\code-repo\vela-tcad\reference_tcad\pn2d_sentaurus2018\source"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$remote = "~/sentaurus_runs/pn2d_sentaurus2018/run_$stamp/source"

ssh sentaurus "mkdir -p $remote"
scp -r "$local\*" "sentaurus:$remote/"

ssh sentaurus "cd $remote && sde -e -l pn2d_sde.cmd"
ssh sentaurus "cd $remote && sdevice pn2d_iv_sdevice.cmd"
ssh sentaurus "cd $remote && sdevice pn2d_bv_sdevice.cmd"

scp "sentaurus:$remote/*.tdr" "$local\"
scp "sentaurus:$remote/*.plt" "$local\"
scp "sentaurus:$remote/*.log" "$local\"
scp "sentaurus:$remote/*.grd" "$local\"
scp "sentaurus:$remote/*.dat" "$local\"
```

上面模板默认使用带时间戳的新目录，便于保留每次运行结果。确认某次远程结果
已经取回且不再需要后，再手动清理对应的 `~/sentaurus_runs/...` 子目录。

## 常见问题

### SSH 报 no matching host key type found

旧 CentOS 只提供 `ssh-rsa` 或 `ssh-dss` 主机密钥。确认
`C:\Users\qzw\.ssh\config` 中有：

```sshconfig
HostKeyAlgorithms +ssh-rsa
```

### ssh sentaurus 仍要求密码

用详细日志确认 Windows 是否送出私钥：

```powershell
ssh -vvv sentaurus 2>&1 | Select-String "Offering public key|Server accepts key|Authentications that can continue|Permission denied"
```

若看到 `Offering public key` 但没有 `Server accepts key`，检查 CentOS：

```bash
ls -ld /home/tcad /home/tcad/.ssh /home/tcad/.ssh/authorized_keys
grep '^ssh-rsa ' /home/tcad/.ssh/authorized_keys
```

权限应为：

```text
/home/tcad                      不应有 group/other 写权限
/home/tcad/.ssh                 drwx------
/home/tcad/.ssh/authorized_keys -rw-------
```

### service sshd status 提示 Permission denied

普通用户 `tcad` 无权管理 SSH 服务。使用 root：

```bash
su -
service sshd status
service sshd start
```

### systemctl 不存在

这是旧 CentOS 的正常现象。使用：

```bash
service sshd status
chkconfig sshd on
```

### 虚拟机突然连不上

先在 VMware 控制台中检查：

```bash
ip route
ifconfig eth0
service sshd status
```

如果 IP 从 `192.168.119.128` 变了，更新 Windows：

```text
C:\Users\qzw\.ssh\config
```

中的 `HostName`。
