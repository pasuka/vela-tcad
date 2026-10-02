#!/usr/bin/env python3
"""Build the canonical portable-report input for SimpleMOS M11."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "simplemos_m11_mobility_factorial_evidence.json")
DEFAULT_OUTPUT = REPO / "docs/validation/reports/simplemos_m11/artifact.json"


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
    findings = evidence["findings"]
    self_effects = findings["self_consistent_aggregate_effects_dex"]
    frozen_effects = findings["frozen_vela_aggregate_effects"]
    variant_rows = findings["self_consistent_variant_aggregate_dex"]

    labels = {
        "phumob": "PhuMob", "enormal": "Enormal", "hfs": "HFS",
        "phumob:enormal": "PhuMob×Enormal",
        "phumob:hfs": "PhuMob×HFS", "enormal:hfs": "Enormal×HFS",
        "phumob:enormal:hfs": "三因素交互",
    }
    aggregate_rows = [{
        "term": labels[term], "effect_dex": value,
        "direction": "误差增加" if value > 0 else "误差降低",
    } for term, value in self_effects[
        "maximum_absolute_log10_ratio_above_floor"].items()]
    condition_rows = []
    for condition in findings["self_consistent_condition_effects_dex"]:
        for term in ("phumob", "enormal", "hfs"):
            condition_rows.append({
                "condition": (f"{condition['device']}, "
                              f"Vd={condition['drain_voltage_V']:g} V"),
                "term": labels[term],
                "effect_dex": condition["effects"]
                    ["maximum_absolute_log10_ratio_above_floor"][term],
            })
    frozen_rows = [{
        "term": labels[term], "effect_dex": value,
        "direction": "误差增加" if value > 0 else "误差降低",
    } for term, value in frozen_effects[
        "mobility_active_p95_error_dex"].items()]
    for row in variant_rows:
        row["models"] = row["variant"].replace("p", "P=").replace(
            "_e", ", E=").replace("_h", ", H=")

    main_effect = self_effects["maximum_absolute_log10_ratio_above_floor"]
    best = min(variant_rows,
               key=lambda row: row["maximum_absolute_log10_ratio_above_floor"])
    full = next(row for row in variant_rows if row["variant"] == "p1_e1_h1")
    headline = [{
        "curves": evidence["data_quality"]["curve_count"],
        "points": evidence["data_quality"]["direct_bias_points_per_solver"],
        "hfs_effect_dex": main_effect["hfs"],
        "enormal_effect_dex": main_effect["enormal"],
        "best_max_error_dex": best["maximum_absolute_log10_ratio_above_floor"],
        "best_vs_full_dex": (full["maximum_absolute_log10_ratio_above_floor"]
                             - best["maximum_absolute_log10_ratio_above_floor"]),
    }]
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")
    title = "SimpleMOS M11 PhuMob × Enormal × HFS 三因素矩阵报告"
    source_path = evidence["artifacts"]["self_consistent_effects"]["path"]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report", "title": title,
            "description": "SDevice-only full factorial mobility-model attribution.",
            "generatedAt": generated,
            "cards": [
                {"id": "curves", "description": "2 个网格、2 档漏压、8 种开关组合。", "dataset": "headline", "sourceId": "m11", "metrics": [{"label": "自洽曲线", "field": "curves", "format": "number"}]},
                {"id": "points", "description": "每个求解器均使用直接栅压点。", "dataset": "headline", "sourceId": "m11", "metrics": [{"label": "每求解器偏置点", "field": "points", "format": "number"}]},
                {"id": "hfs", "description": "正值表示开启 HFS 后最大 Id 误差增加。", "dataset": "headline", "sourceId": "m11", "metrics": [{"label": "HFS 主效应", "field": "hfs_effect_dex", "format": "number", "unit": "dex"}]},
                {"id": "enormal", "description": "负值表示 Enormal 改善最大 Id 一致性。", "dataset": "headline", "sourceId": "m11", "metrics": [{"label": "Enormal 主效应", "field": "enormal_effect_dex", "format": "number", "unit": "dex"}]},
                {"id": "best", "description": "八种组合中跨四工况的最低平均最大误差。", "dataset": "headline", "sourceId": "m11", "metrics": [{"label": "最优组合最大误差", "field": "best_max_error_dex", "format": "number", "unit": "dex"}]},
            ],
            "charts": [
                {"id": "aggregate", "title": "自洽 Id 最大误差的主效应与交互效应", "subtitle": "负值改善一致性；正值增加误差", "type": "bar", "dataset": "aggregate_rows", "sourceId": "m11", "valueFormat": "number", "encodings": {"x": {"field": "term", "type": "nominal", "label": "因子或交互项"}, "y": {"field": "effect_dex", "type": "quantitative", "label": "效应, dex"}, "color": {"field": "direction", "type": "nominal", "label": "方向"}}},
                {"id": "condition", "title": "四个工况的三项主效应", "subtitle": "HFS 在四个工况均增加最大 Id 误差，Enormal 均降低", "type": "bar", "dataset": "condition_rows", "sourceId": "m11", "valueFormat": "number", "encodings": {"x": {"field": "condition", "type": "nominal", "label": "工况"}, "y": {"field": "effect_dex", "type": "quantitative", "label": "最大误差效应, dex"}, "color": {"field": "term", "type": "nominal", "label": "模型"}}},
                {"id": "frozen", "title": "冻结全物理状态的活跃边迁移率 P95 误差效应", "subtitle": "仅在同一 Sentaurus 状态上替换 Vela 迁移率公式", "type": "bar", "dataset": "frozen_rows", "sourceId": "m11", "valueFormat": "number", "encodings": {"x": {"field": "term", "type": "nominal", "label": "因子或交互项"}, "y": {"field": "effect_dex", "type": "quantitative", "label": "P95 误差效应, dex"}, "color": {"field": "direction", "type": "nominal", "label": "方向"}}},
            ],
            "tables": [{"id": "variants", "title": "八种自洽组合的跨工况误差", "subtitle": "最大值和中位值均为四个器件/漏压工况的平均", "dataset": "variant_rows", "sourceId": "m11", "columns": [
                {"field": "variant", "label": "组合", "type": "text"},
                {"field": "maximum_absolute_log10_ratio_above_floor", "label": "平均最大误差, dex", "format": "number"},
                {"field": "median_absolute_log10_ratio_above_floor", "label": "平均中位误差, dex", "format": "number"},
            ]}],
            "sources": [{"id": "m11", "label": "M11 冻结证据", "path": source_path}],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {"id": "summary", "type": "markdown", "body": "## 技术摘要\n\n完整 2³ 矩阵的 32 条 Sentaurus/Vela 自洽曲线全部收敛并逐点对齐。以最大 Id 误差为响应，HFS 是唯一在四个器件/漏压工况中均增加误差的主因子（聚合 +0.00679 dex）；Enormal 在四个工况均改善一致性（聚合 -0.00576 dex）；PhuMob 的聚合效应为 -0.00663 dex，主要来自 n21、Vd=0.05 V。"},
                {"id": "metrics", "type": "metric-strip", "cardIds": ["curves", "points", "hfs", "enormal", "best"]},
                {"id": "self", "type": "markdown", "body": "## 重新自洽求解：HFS 是稳定的误差放大项\n\n自洽矩阵同时包含模型的直接迁移率响应以及电势、载流子和电流的非线性反馈。PhuMob+Enormal、关闭 HFS 的 p1_e1_h0 具有最低的跨工况平均最大误差 0.04693 dex；全模型 p1_e1_h1 为 0.05145 dex。该差异是归因证据，不足以直接改变默认模型。"},
                {"id": "aggregate_chart", "type": "chart", "chartId": "aggregate"},
                {"id": "condition_chart", "type": "chart", "chartId": "condition"},
                {"id": "variant_table", "type": "table", "tableId": "variants"},
                {"id": "frozen_note", "type": "markdown", "body": "## 冻结状态回放：三项模型都能改善全物理迁移率复现\n\n在相同的全模型 Sentaurus 状态上只更换 Vela 迁移率组合时，PhuMob、Enormal 和 HFS 对活跃边迁移率 P95 误差的主效应均为负。PhuMob 主导低场基线，HFS 对高场区继续改善。冻结结果不能解释为 Sentaurus 在同一状态下关闭模型后的响应。"},
                {"id": "frozen_chart", "type": "chart", "chartId": "frozen"},
                {"id": "interpretation", "type": "markdown", "body": "## 综合判断\n\nHFS 在冻结迁移率回放中有益、在自洽跨求解器 Id 对比中却稳定放大误差。这种符号反转说明剩余差异更可能位于 HFS 驱动力的空间构造、边界/参考密度限制、离散映射及其状态反馈，而不是简单的“缺少 HFS”。该结论与 M10 中替换 Sentaurus GradQuasiFermi 后迁移率残差显著下降相互印证。"},
                {"id": "scope", "type": "markdown", "body": "## 范围、方法与限制\n\n范围仅限 SDevice。自洽矩阵使用 n17/n21、Vd=0.05/1 V、Vg=0:0.05:2.5 V，共 1,632 个直接点/求解器，不插值。冻结矩阵使用 M10 的 16 个全物理状态并执行 128 个 Vela 公式评估。Sentaurus 未提供同一冻结状态下的八模型重算，也未公开所有内部 HFS 平滑和极限处理，因此不能据此锁定单一私有公式。"},
                {"id": "next", "type": "markdown", "body": "## 结论与后续\n\n三因素独立消融已完成：Enormal 不是主要误差源；PhuMob 总体改善一致性；HFS 是自洽 Id 误差的稳定放大项。后续应保持当前 16 条全模型曲线为回归基线，针对 GradQuasiFermi 边驱动力、参考密度和平滑限制增加逐边诊断，而不直接修改默认模型。"},
            ],
        },
        "snapshot": {"version": 1, "generatedAt": generated, "status": "ready",
                     "datasets": {"headline": headline,
                                  "aggregate_rows": aggregate_rows,
                                  "condition_rows": condition_rows,
                                  "frozen_rows": frozen_rows,
                                  "variant_rows": variant_rows},
                     "accessIssues": []},
        "sources": [{"id": "m11", "query": {
            "language": "sql", "engine": "duckdb",
            "sql": f"SELECT * FROM read_csv_auto('{source_path}', header=true)",
            "description": "M11 orthogonal factorial effects and frozen evidence.",
            "executed_at": generated, "tables_used": [source_path],
            "metric_definitions": [
                "factorial effect is mean response at coded product +1 minus mean response at -1",
                "negative error effect improves Vela-to-Sentaurus agreement",
                "self-consistent response includes nonlinear state feedback"],
            "filters": ["SDevice only", "n17 and n21", "Vd=0.05 and 1 V",
                        "Vg=0:0.05:2.5 V", "no interpolation"]}}],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8", newline="\n")
    output.with_name("CHART_MAP.md").write_text(
        "# SimpleMOS M11 chart map\n\n"
        "| Report chart | Dataset | Metric | Interpretation |\n"
        "|---|---|---|---|\n"
        "| `aggregate` | `aggregate_rows` | self-consistent maximum Id error effect | main and interaction attribution |\n"
        "| `condition` | `condition_rows` | condition-specific maximum Id error main effect | checks cross-condition sign stability |\n"
        "| `frozen` | `frozen_rows` | frozen-state active-edge mobility P95 error effect | isolates Vela formula response on one state |\n",
        encoding="utf-8", newline="\n")
    print(json.dumps({"status": "built", "output": str(output),
                      "variants": len(variant_rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
