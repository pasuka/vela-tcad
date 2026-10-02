"""Immediate batch dependency coordinator: finalize only after both existing arms finish.

Does not launch or retry solvers. A missing/failed driver stops the coordinator.
"""
import argparse
import ctypes
from ctypes import wintypes
import json
import hashlib
import os
import time
from pathlib import Path
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT/'build-release/hfs_curves_20260914/resume_20260915'
OUT = ROOT/'reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914'


def driver_state(markers, alive):
    """Pure gate used by the process monitor and its independent control tests."""
    if any(m['status'] not in ('running', 'complete') for m in markers):
        return 'driver_failed_or_interrupted'
    if all(m['status'] == 'complete' for m in markers):
        return 'ready'
    if any(m['status'] == 'running' and not alive(m['pid']) for m in markers):
        return 'driver_missing'
    return 'waiting'


def process_alive(pid):
    # Query only. Explicit Win32 signatures avoid truncating 64-bit handles.
    kernel = ctypes.WinDLL('kernel32', use_last_error=True)
    kernel.OpenProcess.argtypes = (wintypes.DWORD, wintypes.BOOL, wintypes.DWORD)
    kernel.OpenProcess.restype = wintypes.HANDLE
    kernel.GetExitCodeProcess.argtypes = (wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD))
    kernel.GetExitCodeProcess.restype = wintypes.BOOL
    kernel.CloseHandle.argtypes = (wintypes.HANDLE,)
    kernel.CloseHandle.restype = wintypes.BOOL
    handle = kernel.OpenProcess(0x1000, False, pid)
    if not handle:
        if ctypes.get_last_error() == 87:  # ERROR_INVALID_PARAMETER: no such process.
            return False
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        code = wintypes.DWORD()
        if not kernel.GetExitCodeProcess(handle, ctypes.byref(code)):
            raise ctypes.WinError(ctypes.get_last_error())
        return code.value == 259  # STILL_ACTIVE
    finally:
        kernel.CloseHandle(handle)


def self_test():
    def m(status, pid=1): return dict(status=status, pid=pid)
    assert driver_state([m('running'), m('running')], lambda _: True) == 'waiting'
    assert driver_state([m('complete'), m('running')], lambda _: True) == 'waiting'
    assert driver_state([m('complete'), m('complete')], lambda _: False) == 'ready'
    assert driver_state([m('complete'), m('running')], lambda _: False) == 'driver_missing'
    assert driver_state([m('complete'), m('failed_or_interrupted')], lambda _: True) == 'driver_failed_or_interrupted'
    assert driver_state([m('running', 1), m('running', 2)], lambda pid: pid == 1) == 'driver_missing'
    assert process_alive(os.getpid())
    print('7 coordinator controls passed; no solvers or evidence modified.', flush=True)


def write_report(e):
    a = e.a
    final = a.read(OUT/'final_summary.json')
    closure = a.read(OUT/'resume_20260915/completion_summary.json')
    lines = ['# HFS 四条完整基线曲线验证（2026-09-15）', '',
        f"完成 408/408 Vela 状态、204/204 全场双初始化、204/204 原生点及密度/SRH 重建。数值失败 {final['failed_attempts']} 次、重载 {final['reloads']} 次。旧 64 态与四个行政中断现场逐字节保持不变。", '',
        '## 电流对比', '', '误差定义：100×(Id_Vela/Id_Sentaurus−1)，电流单位 A/µm。', '',
        '| 器件 | Vd (V) | 最小误差 (%) | 最大误差 (%) | 最大绝对误差对应 Vg (V) |',
        '|---|---:|---:|---:|---:|']
    for c in final['curves']:
        lines.append(f"| {c['device']} | {c['vd']:g} | {c['error_min_percent']:.9g} | {c['error_max_percent']:.9g} | {c['worst_vg']:g} |")
    lines += ['', f"全部电流误差是否小于等于 1%：**{'是' if closure['current_agreement_within_one_percent'] else '否'}**。数值及初始化资格与电流精度分别判定。", '',
        '## 物理场对比', '',
        '以下为每条曲线全部栅压、全部硅节点的最大绝对势差；密度为最大绝对 log10 比值，单位 dex。空间加权 RMS 和逐点结果见字段 CSV。这些为描述指标，不新增物理场验收阈值。', '',
        '| 器件 | Vd (V) | 电势 (µV) | 电子准费米势 (mV) | 空穴准费米势 (mV) | 电子密度 (dex) | 空穴密度 (dex) |',
        '|---|---:|---:|---:|---:|---:|---:|']
    fields = a.rows(OUT/'fields/fields.csv')
    names = ['ElectrostaticPotential', 'eQuasiFermiPotential', 'hQuasiFermiPotential', 'eDensity', 'hDensity']
    for device in ('n19', 'n23'):
        for vd in (.05, 1.):
            values = [max(float(r['max_abs']) for r in fields if r['device']==device and float(r['vd'])==vd and r['field']==name) for name in names]
            values = [values[0]*1e6, values[1]*1e3, values[2]*1e3, *values[3:]]
            lines.append(f'| {device} | {vd:g} | '+' | '.join(f'{x:.9g}' for x in values)+' |')
    lines += ['', '## 配置、资格与限制', '',
        '300 K，原 n19/n23 网格，Vd=0.05/1 V，Vg=0–1 V、步长 0.02 V；Boltzmann、OldSlotboom、PhuMob、Enormal、顶点 HFS 后 box 平均、掺杂相关 SRH。沿用显式几何和 split 状态组合，Eigen SparseLU/COLAMD、四次线性修正。未修改 C++、常数、默认模型或接受门槛。', '',
        f"最大逐行比值 {final['max_row_ratio']:.9g}（门槛 1e-6），最大端口差 {final['max_port_relative']:.9g}（门槛 1e-8）；最大双初始化 Id 相对差 {final['max_dual_Id_relative']:.9g}（门槛 1e-6）。", '',
        f"最大双初始化电势/电子准费米势/空穴准费米势差分别为 {final['max_dual_psi_V']:.9g}/{final['max_dual_phin_V']:.9g}/{final['max_dual_phip_V']:.9g} V（门槛 1e-6 V），密度相对差 {final['max_dual_density_relative']:.9g}（门槛 1e-4）。", '',
        '共同 barycentric Si 体积的 SRH 对比只描述局部复合率及状态差，不等于两求解器各自源积分或端口响应。当前基线仍使用 all_cell SRH 体积；原生 signed Si 体积证据、微小空穴迁移率截断差和大步长跨分支 Jv 失败记录继续保留。未进行 SRH 体积修复，也未放行原始 n17–n24、0–2.5 V 全矩阵。', '',
        '下一阶段按已确认顺序，开展独立 SRH 体积扰动与自洽 A/B；不得把固定状态源账本当作漏极电流响应校准。', '',
        '## 证据', '',
        '- 数据目录：`reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914`。完整电流表 `comparison.csv`，字段 `fields/fields.csv`，SRH `fields/srh.csv`，逐节点场 `fields/nodes.csv`。',
        '- 曲线图 `curves.png` / `curves.pdf`；完整程序和源码归档见 `qualified_source_archive.json`。',
        '- 续算闭合：`resume_20260915/completion_summary.json`、`completion_evidence.json`；行政中断证据仍在 `pause_20260914_2213`。',
        '- [续算过程](simplemos_hfs_curve_resume_2026-09-15.md)、[SRH 积分审计](simplemos_hfs_srh_quadrature_audit_2026-09-14.md)。', '']
    report = ROOT/'docs/validation/simplemos_hfs_fullcurve_validation_2026-09-15.md'
    with report.open('x', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))
    e.d.matrix.freeze(OUT/'resume_20260915/report_evidence.json', [report,
        OUT/'resume_20260915/completion_evidence.json', OUT/'fields/evidence.json', Path(__file__).resolve()])
    status = ROOT/'docs/validation/simplemos_branch_status.md'
    content = status.read_text(encoding='utf-8')
    heading = '\n\n**9 月 15 日 HFS 完整基线已完成**：408/408 状态、204/204 全场双初始化与重建通过；电流精度是否满足 1% 见[完整报告](simplemos_hfs_fullcurve_validation_2026-09-15.md)。SRH 修复与原始全矩阵尚未执行。下方续算描述保留历史状态。\n'
    first, rest = content.split('\n', 1)
    status.write_text(first+heading+rest, encoding='utf-8')


def main():
    manifest = json.loads((OUT/'resume_20260915/postprocessing_evidence.json').read_text(encoding='utf-8'))
    for name, digest in manifest['input_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    marker = LOCAL/'completion_active.json'
    def mark(status, **extra):
        marker.write_text(json.dumps(dict(pid=os.getpid(), status=status,
            recorded_at=datetime.now().astimezone().isoformat(), **extra), indent=2), encoding='utf-8')
    # Exclusive start, like the two solvers' drivers. No duplicate coordinators.
    with marker.open('x', encoding='utf-8') as stream:
        json.dump(dict(pid=os.getpid(), status='starting'), stream)
    try:
        while True:
            markers = []
            try:
                markers = [json.loads((LOCAL/(arm+'_active.json')).read_text(encoding='utf-8')) for arm in ('continuation', 'native')]
            except json.JSONDecodeError:
                time.sleep(1)
                continue
            state = driver_state(markers, process_alive)
            mark(state, drivers=markers)
            if state == 'ready':
                break
            if state != 'waiting':
                print(state+'; no postprocessing launched.', flush=True)
                return
            time.sleep(15)
        import simplemos_hfs_curves_20260914 as e
        import summarize_simplemos_hfs_curve_fields_20260914 as fields
        import finish_simplemos_hfs_curves_20260914 as finish
        import seal_simplemos_hfs_resume_20260915 as seal
        for phase, action in [('analysis', lambda: (e.configure(), e.v.analyze())),
                              ('fields', fields.main), ('archive_and_plots', finish.main),
                              ('seal', seal.main), ('report', lambda: write_report(e))]:
            mark(phase)
            print('Starting '+phase, flush=True)
            action()
        mark('complete', report='docs/validation/simplemos_hfs_fullcurve_validation_2026-09-15.md')
        (OUT/'execution_status.json').write_text(marker.read_text(encoding='utf-8'), encoding='utf-8')
    except BaseException as error:
        mark('failed_or_interrupted', error=repr(error))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['self-test', 'run'])
    self_test() if parser.parse_args().action == 'self-test' else main()
