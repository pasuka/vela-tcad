#!/usr/bin/env python3
"""Build the canonical portable-report input for SimpleMOS M9."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m9_hfs_diagnostics_evidence.json"
)
DEFAULT_OUTPUT = REPO / "docs/validation/reports/simplemos_m9/artifact.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows: list[dict[str, Any]] = []
        for raw in csv.DictReader(stream):
            row: dict[str, Any] = dict(raw)
            for key, value in raw.items():
                if value in (None, ""):
                    continue
                try:
                    row[key] = float(value)
                except ValueError:
                    pass
            rows.append(row)
        return rows


def source_path(evidence: dict[str, Any], key: str) -> Path:
    return REPO / evidence["artifacts"][key]["path"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    evidence = read_json(args.evidence.resolve())
    comparison = read_json(source_path(evidence, "comparison_report"))
    parameter = read_json(source_path(evidence, "parameter_audit"))
    edge = read_json(source_path(evidence, "edge_report"))
    scan_path = evidence["artifacts"]["scan_summary"]["path"]
    edge_path = evidence["artifacts"]["edge_summary"]["path"]
    parameter_path = evidence["artifacts"]["parameter_audit"]["path"]
    scan_rows = []
    for row in read_csv(REPO / scan_path):
        if float(row["refdens_cm3"]) <= 0.0:
            continue
        scan_rows.append({
            **row,
            "refdens_log10_cm3": math.log10(float(row["refdens_cm3"])),
            "condition": f"{row['device']}, Vd={float(row['drain_voltage_V']):g} V",
        })
    edge_rows = []
    for row in read_csv(REPO / edge_path):
        edge_rows.append({
            **row,
            "condition": f"{row['device']}, Vd={float(row['drain_voltage_V']):g} V",
        })
    parameter_rows = parameter["parameter_rows"]
    minimum_electron_limiter = min(float(row[
        "minimum_electron_mobility_limiter"]) for row in edge_rows)
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")
    findings = evidence["findings"]
    title = "SimpleMOS M9 HighFieldSaturation 参考密度与边级诊断报告"
    headline = [{
        "sentaurus_curves": evidence["execution"]["sentaurus_curves"],
        "direct_bias_points": evidence["execution"]["total_direct_bias_points"],
        "material_response_refdens_cm3": findings[
            "material_response_from_refdens_cm3"],
        "maximum_scan_response_dex": findings["maximum_scan_response_dex"],
        "minimum_electron_limiter": minimum_electron_limiter,
        "baseline_cases": evidence["baseline_guard"]["case_count"],
    }]
    summary_source = {"id": "refdens_scan", "label": "M9 冻结参考密度扫描",
                      "path": scan_path}
    edge_source = {"id": "edge_diagnostics", "label": "M9 冻结边级 HFS 摘要",
                   "path": edge_path}
    parameter_source = {"id": "parameter_audit", "label": "M9 HFS 参数审计",
                        "path": parameter_path}
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report", "title": title,
            "description": "SDevice-only HFS density interpolation and Vela per-edge limiter diagnosis.",
            "generatedAt": generated,
            "cards": [
                {"id": "curve_count", "description": "六个配置、四个器件/漏压条件。",
                 "dataset": "headline", "sourceId": "refdens_scan",
                 "metrics": [{"label": "Sentaurus 曲线", "field": "sentaurus_curves", "format": "number"}]},
                {"id": "point_count", "description": "全部为直接求解，禁止插值。",
                 "dataset": "headline", "sourceId": "refdens_scan",
                 "metrics": [{"label": "直接偏置点", "field": "direct_bias_points", "format": "number"}]},
                {"id": "material_density", "description": "四工况均达到 1e-4 dex 响应的最低密度。",
                 "dataset": "headline", "sourceId": "refdens_scan",
                 "metrics": [{"label": "物质性起点", "field": "material_response_refdens_cm3", "format": "number", "unit": "cm^-3"}]},
                {"id": "max_response", "description": "扫描中最大的 Sentaurus 电流响应。",
                 "dataset": "headline", "sourceId": "refdens_scan",
                 "metrics": [{"label": "最大响应", "field": "maximum_scan_response_dex", "format": "number", "unit": "dex"}]},
                {"id": "edge_limiter", "description": "16 个状态中最强的电子 HFS 限制。",
                 "dataset": "headline", "sourceId": "edge_diagnostics",
                 "metrics": [{"label": "最小 μ/μ0", "field": "minimum_electron_limiter", "format": "number"}]},
                {"id": "baseline", "description": "原始 M8 曲线及哈希保持不变。",
                 "dataset": "headline", "sourceId": "refdens_scan",
                 "metrics": [{"label": "冻结基线曲线", "field": "baseline_cases", "format": "number"}]},
            ],
            "charts": [
                {"id": "scan_chart", "title": "参考密度对最大比较误差的影响",
                 "subtitle": "负值为改善；横轴为 log10(RefDens/cm^-3)",
                 "type": "line", "dataset": "scan_rows", "sourceId": "refdens_scan",
                 "valueFormat": "number", "encodings": {
                     "x": {"field": "refdens_log10_cm3", "type": "quantitative", "label": "log10 参考密度"},
                     "y": {"field": "maximum_error_change_dex", "type": "quantitative", "label": "最大误差变化, dex"},
                     "color": {"field": "condition", "type": "nominal", "label": "器件与漏压"},
                     "tooltip": [
                         {"field": "maximum_sentaurus_response_dex", "type": "quantitative", "label": "最大 Sentaurus 响应, dex"},
                         {"field": "refdens_cm3", "type": "quantitative", "label": "参考密度, cm^-3"},
                     ]}},
                {"id": "limiter_chart", "title": "Vela 边级电子 HFS 等效限制因子",
                 "subtitle": "每个状态取传输边最小 μ/μ0；越低表示速度饱和越强",
                 "type": "line", "dataset": "edge_rows", "sourceId": "edge_diagnostics",
                 "valueFormat": "number", "encodings": {
                     "x": {"field": "gate_voltage_V", "type": "quantitative", "label": "栅压, V"},
                     "y": {"field": "minimum_electron_mobility_limiter", "type": "quantitative", "label": "最小 μ/μ0"},
                     "color": {"field": "condition", "type": "nominal", "label": "器件与漏压"},
                     "tooltip": [
                         {"field": "maximum_electron_drive_V_per_cm", "type": "quantitative", "label": "最大电子驱动力, V/cm"},
                         {"field": "maximum_electron_high_field_ratio", "type": "quantitative", "label": "最大 μ0F/vsat"},
                     ]}},
            ],
            "tables": [
                {"id": "parameter_table", "title": "300 K HFS 参数逐项审计",
                 "subtitle": "Sentaurus sdevice -P 导出与 Vela 默认值",
                 "dataset": "parameter_rows", "sourceId": "parameter_audit",
                 "defaultSort": {"field": "carrier", "direction": "asc"},
                 "columns": [
                     {"field": "carrier", "label": "载流子", "type": "text"},
                     {"field": "parameter", "label": "参数", "type": "text"},
                     {"field": "sentaurus", "label": "Sentaurus", "format": "number"},
                     {"field": "vela", "label": "Vela", "format": "number"},
                     {"field": "status", "label": "结论", "type": "text"},
                     {"field": "relevance_at_300K", "label": "300 K 相关性", "type": "text"},
                 ]},
                {"id": "scan_table", "title": "20 组正参考密度逐工况结果",
                 "subtitle": "固定 Vela；只改变 Sentaurus HFS 参考密度",
                 "dataset": "scan_rows", "sourceId": "refdens_scan",
                 "defaultSort": {"field": "maximum_error_change_dex", "direction": "desc"},
                 "columns": [
                     {"field": "device", "label": "器件", "type": "text"},
                     {"field": "drain_voltage_V", "label": "Vd, V", "format": "number"},
                     {"field": "refdens_cm3", "label": "RefDens, cm^-3", "format": "number"},
                     {"field": "maximum_sentaurus_response_dex", "label": "最大响应, dex", "format": "number"},
                     {"field": "maximum_error_change_dex", "label": "最大误差变化, dex", "format": "number"},
                 ]},
            ],
            "sources": [summary_source, edge_source, parameter_source],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {"id": "summary", "type": "markdown", "body":
                 "## 技术摘要\n\n公开的 300 K Caughey–Thomas 核心公式及 `vsat0/beta0/alpha` 与 Vela 一致；剩余差异不来自这些核心常数。参考密度扫描显示，`1e2 cm⁻³` 的平均改善仅为 -2.25e-6 dex，低于物质性阈值；从 `1e6 cm⁻³` 起四个工况全部恶化，`1e10 cm⁻³` 的平均最大误差增加 0.488 dex。因此参考密度插值是可辨识的影响源，但不应作为 Vela 默认参数拟合。原 16 条/816 点基线哈希保持不变。"},
                {"id": "metrics", "type": "metric-strip", "cardIds": ["curve_count", "point_count", "material_density", "max_response", "edge_limiter", "baseline"]},
                {"id": "finding_scan", "type": "markdown", "sourceId": "refdens_scan", "body":
                 "## 参考密度只在接近零时近似空对照，增大后系统性恶化\n\n扫描覆盖 `0、1e2、1e4、1e6、1e8、1e10 cm⁻³`。显式 GradQuasiFermi 与默认曲线完全相同；`1e2` 的微小改善没有达到 1e-4 dex 物质性门槛。`1e6、1e8、1e10` 在四个工况均提高最大误差，且响应随密度单调放大。"},
                {"id": "scan_chart_block", "type": "chart", "chartId": "scan_chart"},
                {"id": "parameters", "type": "markdown", "sourceId": "parameter_audit", "body":
                 "## 核心饱和速度与指数不是当前残差来源\n\nSentaurus 导出电子/空穴 `vsat0=1.07e7/8.37e6 cm/s`、`beta0=1.109/1.213`、`alpha=0`；换算后与 Vela 默认值及公式逐点完全一致。`betaexp` 与 `vsatexp` 在固定 300 K 时为中性，`K_dT` 属于载流子温度分支。`E0_TrEf/Ksmooth_TrEf` 是否触发额外私有分支仍无法由导出文件单独证明。"},
                {"id": "parameter_table_block", "type": "table", "tableId": "parameter_table"},
                {"id": "edge", "type": "markdown", "sourceId": "edge_diagnostics", "body":
                 "## 边级诊断把驱动力、饱和参数和限制过程连成可审计链\n\n16 个状态覆盖 n17/n21、两档漏压及 Vg=0/0.05/0.8/2.5 V。诊断现在输出逐边载流子密度、GradQF 驱动力、低场迁移率、饱和速度、beta、`μ0F/vsat`、等效限制因子和重构误差。生产等效限制重构误差最大 2.22e-16；相邻单元非线性限制后再平均造成的最大均值近似残差为 0.619%，已单独标识。"},
                {"id": "limiter_chart_block", "type": "chart", "chartId": "limiter_chart"},
                {"id": "scope", "type": "markdown", "body":
                 "## 范围、数据与方法\n\n范围仅限 SDevice。Sentaurus 使用 T-2022.03-SP2、两份固定 TDR、Vd=0.05/1 V 与 Vg=0:0.05:2.5 V 的 51 点精确栅压网格；禁止插值。比较固定 M8-A Vela 曲线，只改变 Sentaurus 的 HFS 参考密度。误差为电流绝对值对数比，电流门槛 1e-18 A/μm。Vela 边级输出为只读诊断，不参与方程或默认模型。"},
                {"id": "detail_heading", "type": "markdown", "body":
                 "## 逐工况扫描结果\n\n最大响应衡量 Sentaurus 选项本身的可辨识度；最大误差变化衡量它对固定 Vela 比较的方向。"},
                {"id": "scan_table_block", "type": "table", "tableId": "scan_table"},
                {"id": "limitations", "type": "markdown", "body":
                 "## 局限性、稳健性与未解析项\n\n本轮证明参考密度插值能够稳定改变结果，却不能反演 Sentaurus 私有的 GradQuasiFermi 构造、平滑或极限公式。高密度工况的强响应主要集中在深关断/弱反型，亚飞安电流更敏感。边级最小限制因子是全器件极值，不等同于沟道加权的端电流贡献。参数审计仅针对 300 K 等温漂移扩散。"},
                {"id": "recommendations", "type": "markdown", "body":
                 "## 建议的下一步\n\n1. 保持 Vela 默认 HFS 参数不变，并把 16 条原始曲线继续作为回归基线。\n2. 使用新增逐边输出定位对端电流贡献最大的沟道边，再比较驱动力离散方式，而不是调大 RefDens。\n3. 若继续审查 Sentaurus 私有行为，应优先设计温度扫描或空间探针来辨识 `betaexp/vsatexp` 与 `E0_TrEf/Ksmooth_TrEf` 的激活条件。\n4. 将本轮 `1e6–1e10 cm⁻³` 系统性恶化作为防止误用参考密度拟合的负向回归证据。"},
                {"id": "questions", "type": "markdown", "body":
                 "## 后续问题\n\n剩余基线误差中，GradQF 的空间离散、接触/边界处极限处理和私有低密度平滑各占多少；以及沟道加权边诊断能否把约 0.1 dex 峰值误差缩小到单一空间区域，仍需进一步验证。"},
            ],
        },
        "snapshot": {"version": 1, "generatedAt": generated, "status": "ready",
                     "datasets": {"headline": headline, "scan_rows": scan_rows,
                                  "edge_rows": edge_rows, "parameter_rows": parameter_rows},
                     "accessIssues": []},
        "sources": [
            {"id": "refdens_scan", "query": {"language": "sql", "engine": "duckdb",
                "sql": f"SELECT * FROM read_csv_auto('{scan_path}', header=true)",
                "description": "Twenty positive-reference-density comparisons plus frozen null control in the comparison report.",
                "executed_at": generated, "tables_used": [scan_path],
                "metric_definitions": ["maximum_error_change_dex is controlled maximum absolute log-current error minus the explicit-GradQF baseline", "maximum_sentaurus_response_dex is max abs(log10(abs(Id_RefDens)/abs(Id_default))) over 51 exact gate points"],
                "filters": ["n17 and n21", "Vd=0.05 and 1 V", "Vg=0:0.05:2.5 V", "no interpolation", "current floor 1e-18 A/um"]}},
            {"id": "edge_diagnostics", "query": {"language": "sql", "engine": "duckdb",
                "sql": f"SELECT * FROM read_csv_auto('{edge_path}', header=true)",
                "description": "Sixteen Vela per-edge HFS diagnostic state summaries.",
                "executed_at": generated, "tables_used": [edge_path],
                "metric_definitions": ["minimum_electron_mobility_limiter is min production edge mobility divided by zero-drive mobility", "reconstruction error compares the exact effective edge denominator against production mobility"],
                "filters": ["n17 and n21", "Vd=0.05 and 1 V", "Vg in 0, 0.05, 0.8, 2.5 V", "transport-capable edges only"]}},
            {"id": "parameter_audit", "query": {"language": "sql", "engine": "duckdb",
                "sql": f"SELECT * FROM read_csv_auto('{str(Path(parameter_path).with_name('hfs_parameter_comparison.csv')).replace(chr(92), '/')}', header=true)",
                "description": "Sentaurus sdevice -P HighFieldDependence parameters compared with Vela defaults.",
                "executed_at": generated, "tables_used": [parameter_path],
                "metric_definitions": ["vsat0 is converted from cm/s to m/s before comparison", "core formula equivalence is evaluated at 300 K over a broad dimensionless-field grid"],
                "filters": ["Sentaurus T-2022.03-SP2", "SimpleMOS n17 explicit GradQuasiFermi", "300 K isothermal DD"]}},
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "built", "output": str(output),
                      "datasets": {key: len(value) for key, value in artifact[
                          "snapshot"]["datasets"].items()}}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
