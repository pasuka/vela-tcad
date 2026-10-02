#!/usr/bin/env python3
"""Build the canonical portable-report input for SimpleMOS M10."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "simplemos_m10_fixed_state_replay_evidence.json")
DEFAULT_OUTPUT = REPO / "docs/validation/reports/simplemos_m10/artifact.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for raw in csv.DictReader(stream):
            row: dict[str, Any] = dict(raw)
            for key, value in raw.items():
                if value not in (None, ""):
                    try:
                        row[key] = float(value)
                    except ValueError:
                        pass
            rows.append(row)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    evidence = read_json(args.evidence.resolve())
    aggregate = evidence["findings"]
    rows = read_csv(REPO / evidence["artifacts"]["summary"]["path"])
    for row in rows:
        row["condition"] = (f"{row['device']}, Vd={row['drain_voltage_V']:g} V, "
                            f"Vg={row['gate_voltage_V']:g} V")
    mobility_rows = []
    for row in rows:
        mobility_rows.extend([
            {"state": row["state"], "condition": row["condition"],
             "stage": "Vela drive", "p95_error_dex": row[
                 "mobility_vela_drive_active_p95_error_dex"]},
            {"state": row["state"], "condition": row["condition"],
             "stage": "Sentaurus drive", "p95_error_dex": row[
                 "mobility_sentaurus_drive_active_p95_error_dex"]},
        ])
    terminal_rows = []
    for row in rows:
        terminal_rows.extend([
            {"state": row["state"], "condition": row["condition"],
             "stage": "Vela drive", "error_dex": row[
                 "terminal_vela_drive_error_dex"]},
            {"state": row["state"], "condition": row["condition"],
             "stage": "Sentaurus drive", "error_dex": row[
                 "terminal_sentaurus_drive_error_dex"]},
            {"state": row["state"], "condition": row["condition"],
             "stage": "Sentaurus final mobility", "error_dex": row[
                 "terminal_sentaurus_mobility_error_dex"]},
        ])
    finding = aggregate["mobility_improvement"]
    strong = aggregate["terminal_current_strong_inversion"]["vela_drive_vela_hfs"]
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")
    title = "SimpleMOS M10 冻结状态中间物理量交叉回放报告"
    headline = [{
        "states": evidence["execution"]["state_count"],
        "median_reduction_percent": 100.0 * finding["median_reduction_fraction"],
        "p95_reduction_percent": 100.0 * finding["p95_reduction_fraction"],
        "mapping_sensitivity_dex": aggregate[
            "mapping_sensitivity_max_p95_difference_dex"],
        "strong_inversion_max_terminal_error_dex": strong["maximum_dex"],
    }]
    summary_path = evidence["artifacts"]["summary"]["path"]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report", "title": title,
            "description": "SDevice-only fixed-state replay of Sentaurus GradQuasiFermi and mobility in Vela.",
            "generatedAt": generated,
            "cards": [
                {"id": "states", "description": "2 个网格、2 档漏压、4 个栅压。",
                 "dataset": "headline", "sourceId": "m10_summary",
                 "metrics": [{"label": "精确冻结状态", "field": "states", "format": "number"}]},
                {"id": "median", "description": "活跃边迁移率误差的平均中位数降幅。",
                 "dataset": "headline", "sourceId": "m10_summary",
                 "metrics": [{"label": "中位误差降幅", "field": "median_reduction_percent", "format": "number", "unit": "%"}]},
                {"id": "p95", "description": "活跃边迁移率误差的平均 P95 降幅。",
                 "dataset": "headline", "sourceId": "m10_summary",
                 "metrics": [{"label": "P95 误差降幅", "field": "p95_reduction_percent", "format": "number", "unit": "%"}]},
                {"id": "mapping", "description": "两种节点到边映射的最大 P95 差。",
                 "dataset": "headline", "sourceId": "m10_summary",
                 "metrics": [{"label": "映射敏感度", "field": "mapping_sensitivity_dex", "format": "number", "unit": "dex"}]},
                {"id": "terminal", "description": "Vg≥0.8 V 的冻结状态 SG 端电流。",
                 "dataset": "headline", "sourceId": "m10_summary",
                 "metrics": [{"label": "强反型最大端电流误差", "field": "strong_inversion_max_terminal_error_dex", "format": "number", "unit": "dex"}]},
            ],
            "charts": [
                {"id": "mobility", "title": "逐状态活跃边电子迁移率 P95 误差",
                 "subtitle": "相同 Sentaurus 状态，只替换 HFS 驱动力",
                 "type": "bar", "dataset": "mobility_rows", "sourceId": "m10_summary",
                 "valueFormat": "number", "encodings": {
                     "x": {"field": "condition", "type": "nominal", "label": "冻结状态"},
                     "y": {"field": "p95_error_dex", "type": "quantitative", "label": "P95 误差, dex"},
                     "color": {"field": "stage", "type": "nominal", "label": "回放阶段"}}},
                {"id": "terminal", "title": "端电流回放误差",
                 "subtitle": "深关断点反映极小电流条件数，不用于 HFS 定量归因",
                 "type": "bar", "dataset": "terminal_rows", "sourceId": "m10_summary",
                 "valueFormat": "number", "encodings": {
                     "x": {"field": "condition", "type": "nominal", "label": "冻结状态"},
                     "y": {"field": "error_dex", "type": "quantitative", "label": "端电流误差, dex"},
                     "color": {"field": "stage", "type": "nominal", "label": "回放阶段"}}},
            ],
            "tables": [{
                "id": "states_table", "title": "16 个冻结状态的核心指标",
                "subtitle": "所有栅压均来自直接求解快照，未插值",
                "dataset": "state_rows", "sourceId": "m10_summary",
                "columns": [
                    {"field": "device", "label": "网格", "type": "text"},
                    {"field": "drain_voltage_V", "label": "Vd, V", "format": "number"},
                    {"field": "gate_voltage_V", "label": "Vg, V", "format": "number"},
                    {"field": "mobility_vela_drive_active_p95_error_dex", "label": "Vela 驱动 P95, dex", "format": "number"},
                    {"field": "mobility_sentaurus_drive_active_p95_error_dex", "label": "Sentaurus 驱动 P95, dex", "format": "number"},
                    {"field": "terminal_vela_drive_error_dex", "label": "端电流误差, dex", "format": "number"},
                ]}],
            "sources": [{"id": "m10_summary", "label": "M10 冻结状态汇总", "path": summary_path}],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {"id": "summary", "type": "markdown", "body":
                 "## 技术摘要\n\n16 个冻结状态全部完成直接快照、字段导出和交叉回放。用 Sentaurus 的 GradQuasiFermi 驱动力替换 Vela 自身驱动力后，活跃边电子迁移率的平均中位误差从 0.0230 dex 降至 0.00838 dex，平均 P95 从 0.2365 dex 降至 0.1641 dex，且 16/16 状态均改善。这证明驱动力构造/离散是剩余误差的稳定贡献源，但不能解释全部残差。"},
                {"id": "metrics", "type": "metric-strip", "cardIds": ["states", "median", "p95", "mapping", "terminal"]},
                {"id": "finding", "type": "markdown", "sourceId": "m10_summary", "body":
                 "## 替换驱动力产生稳定、可重复的迁移率改善\n\n回放保留 Sentaurus 电势、准费米势和载流子浓度状态，并调用 Vela 生产 HFS 公式。两种节点到边映射的最大 P95 差低于 3.3e-4 dex，说明改善不是某一种映射选择的偶然结果。"},
                {"id": "mobility_chart", "type": "chart", "chartId": "mobility"},
                {"id": "terminal_note", "type": "markdown", "body":
                 "## 端电流提供一致性检查，而非全偏置 HFS 归因\n\nVg≥0.8 V 的八个强反型状态最大端电流误差低于 4.3e-4 dex。深关断状态的 Sentaurus 电流处于亚飞安量级，净电流由大分量抵消得到，节点到边映射和离散差异会被放大至 0.2–2.8 dex；这些点不适合衡量 HFS 参数优劣。"},
                {"id": "terminal_chart", "type": "chart", "chartId": "terminal"},
                {"id": "scope", "type": "markdown", "body":
                 "## 范围、数据与方法\n\n范围仅限 SDevice。使用 n17/n21 两个 SimpleMOS 网格、Vd=0.05/1 V、Vg=0/0.05/0.8/2.5 V。Sentaurus T-2022.03-SP2 依次求解并保存 16 个精确状态，导出电势、准费米势、载流子、GradQuasiFermi、迁移率、速度与电流密度；Vela 的外部驱动力接口仅用于只读诊断，不改变默认模型或方程。"},
                {"id": "table", "type": "table", "tableId": "states_table"},
                {"id": "limitations", "type": "markdown", "body":
                 "## 局限性与剩余误差\n\nSentaurus 输出的是节点后处理量，Vela 使用有限体积边量，因此必须做节点到边映射。Sentaurus 未单独导出 HFS 前的低场迁移率，当前无法把剩余迁移率误差严格拆成 PhuMob/Enormal 与 HFS 私有平滑两部分。投影的 Sentaurus 节点电流仅用于空间定位，不是严格的离散恒等式。"},
                {"id": "next", "type": "markdown", "body":
                 "## 结论与下一步\n\nGradQuasiFermi 驱动力构造/离散已被确认是稳定贡献源；替换最终 Sentaurus 迁移率并未改善深关断端电流，说明这些点的主导问题不在 HFS。若继续缩小迁移率残差，应增加 Sentaurus 低场迁移率对照输出或关闭 HFS 的同状态配对快照，再将 PhuMob、Enormal 和 HFS 限制逐层拆分。"},
            ],
        },
        "snapshot": {"version": 1, "generatedAt": generated, "status": "ready",
                     "datasets": {"headline": headline, "state_rows": rows,
                                  "mobility_rows": mobility_rows,
                                  "terminal_rows": terminal_rows},
                     "accessIssues": []},
        "sources": [{"id": "m10_summary", "query": {
            "language": "sql", "engine": "duckdb",
            "sql": f"SELECT * FROM read_csv_auto('{summary_path}', header=true)",
            "description": "Sixteen exact SimpleMOS fixed-state replay summaries.",
            "executed_at": generated, "tables_used": [summary_path],
            "metric_definitions": [
                "mobility errors are absolute log10 ratios on active transport edges",
                "active edges exceed 0.1 percent of the state maximum projected Sentaurus current",
                "terminal-current error is absolute log10 ratio in A/um"],
            "filters": ["n17 and n21", "Vd=0.05 and 1 V",
                        "Vg=0,0.05,0.8,2.5 V", "no interpolation"]}}],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8", newline="\n")
    chart_map = output.with_name("CHART_MAP.md")
    chart_map.write_text(
        "# SimpleMOS M10 chart map\n\n"
        "| Report chart | Dataset | Metric | Interpretation |\n"
        "|---|---|---|---|\n"
        "| `mobility` | `mobility_rows` | active-edge electron mobility P95 absolute log10 error | isolates response to GradQuasiFermi drive substitution |\n"
        "| `terminal` | `terminal_rows` | drain-cut total-current absolute log10 error | strong-inversion consistency check; deep-off conditioning diagnostic |\n",
        encoding="utf-8", newline="\n")
    print(json.dumps({"status": "built", "output": str(output),
                      "states": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
