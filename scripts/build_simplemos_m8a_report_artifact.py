#!/usr/bin/env python3
"""Build the canonical portable-report input for SimpleMOS M8-A."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8a_model_ablation_evidence.json"
)
DEFAULT_OUTPUT = (
    REPO / "docs/validation/reports/simplemos_m8a/artifact.json"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = []
        for row in csv.DictReader(handle):
            parsed = dict(row)
            for key, value in row.items():
                if key in {"variant", "device", "status"}:
                    continue
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    pass
            rows.append(parsed)
        return rows


def source_path(evidence: dict[str, Any], key: str) -> Path:
    return REPO / evidence["artifacts"][key]["path"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    evidence_path = args.evidence.resolve()
    evidence = read_json(evidence_path)
    first = read_csv(source_path(evidence, "first_round_summary"))
    confirmation = read_csv(source_path(evidence, "confirmation_summary"))
    first_by_variant = {row["variant"]: row for row in first}
    no_hfs_controls = [row for row in confirmation if row["variant"] == "no_hfs"]
    phumob_controls = [row for row in confirmation if row["variant"] == "phumob_only"]
    consistent_count = sum(row["maximum_change_vs_full_dex"] < 0.0
                           for row in no_hfs_controls)
    full = first_by_variant["full"]
    no_hfs = first_by_variant["no_hfs"]
    no_bgn = first_by_variant["no_bgn"]
    headline = [{
        "full_max_dex": full["maximum_absolute_log10_ratio"],
        "no_hfs_max_dex": no_hfs["maximum_absolute_log10_ratio"],
        "no_hfs_improvement_dex": -no_hfs["maximum_change_vs_full_dex"],
        "no_bgn_max_dex": no_bgn["maximum_absolute_log10_ratio"],
        "consistent_conditions": consistent_count,
        "condition_count": len(no_hfs_controls),
    }]
    label_map = {
        "no_srh": "No SRH", "plain_srh": "Plain SRH",
        "no_bgn": "No OldSlotboom", "no_enormal": "No Enormal",
        "no_hfs": "No HFS", "phumob_only": "PhuMob only",
        "constant_mu": "Constant mobility",
    }
    first_zoom = [{
        **row, "variant_label": label_map[row["variant"]],
    } for row in first if row["variant"] not in {"full", "no_bgn"}]
    confirmation_chart = []
    for row in confirmation:
        if row["variant"] not in {"no_hfs", "phumob_only"}:
            continue
        confirmation_chart.append({
            **row,
            "condition": f"{row['device']}, Vd={row['drain_voltage_V']:g} V",
            "variant_label": "No HFS" if row["variant"] == "no_hfs"
            else "PhuMob only",
        })
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")
    first_source_path = evidence["artifacts"]["first_round_summary"]["path"]
    confirmation_source_path = evidence["artifacts"]["confirmation_summary"]["path"]
    title = "SimpleMOS M8-A 器件模型误差归因报告"
    hfs_all = consistent_count == len(no_hfs_controls)
    conclusion = (
        "HFS 关闭在四个控制工况中均降低最大误差，支持其为剩余误差的稳定贡献源。"
        if hfs_all else
        f"HFS 关闭仅在 {consistent_count}/{len(no_hfs_controls)} 个控制工况降低最大误差，"
        "因此它是 n23 残差的贡献因素，但尚不能视为跨工况唯一根因。"
    )
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report", "title": title,
            "description": "SimpleMOS SDevice-only paired model-ablation diagnosis.",
            "generatedAt": generated,
            "cards": [
                {"id": "full_error", "description": "n23、Vd=0.05 V 全物理基线的最大逐点误差。",
                 "dataset": "headline", "sourceId": "first_round",
                 "metrics": [{"label": "Full 最大误差", "field": "full_max_dex",
                              "format": "number", "unit": "dex"}]},
                {"id": "no_hfs_error", "description": "成对关闭 HFS 后的 n23 最大误差。",
                 "dataset": "headline", "sourceId": "first_round",
                 "metrics": [{"label": "No HFS 最大误差", "field": "no_hfs_max_dex",
                              "format": "number", "unit": "dex"}]},
                {"id": "hfs_gain", "description": "相对 full 的最大误差下降量。",
                 "dataset": "headline", "sourceId": "first_round",
                 "metrics": [{"label": "HFS 关闭改善", "field": "no_hfs_improvement_dex",
                              "format": "number", "unit": "dex"}]},
                {"id": "no_bgn_error", "description": "关闭 OldSlotboom 后的最大跨求解器误差。",
                 "dataset": "headline", "sourceId": "first_round",
                 "metrics": [{"label": "No OldSlotboom", "field": "no_bgn_max_dex",
                              "format": "number", "unit": "dex"}]},
            ],
            "charts": [
                {"id": "first_round_chart", "title": "第一轮模型消融的误差变化",
                 "subtitle": "n23、Vd=0.05 V；为显示小差异，图中不含 1.600 dex 的 No OldSlotboom 离群项",
                 "type": "bar", "dataset": "first_round_zoom", "sourceId": "first_round",
                 "valueFormat": "number",
                 "encodings": {
                     "x": {"field": "variant_label", "type": "nominal", "label": "模型变体"},
                     "y": {"field": "maximum_change_vs_full_dex", "type": "quantitative",
                           "label": "最大绝对对数误差相对 full 的变化, dex"},
                     "tooltip": [
                         {"field": "maximum_absolute_log10_ratio", "type": "quantitative",
                          "label": "最大误差, dex"},
                         {"field": "maximum_relative_error", "type": "quantitative",
                          "label": "最大相对误差"},
                     ],
                 }},
                {"id": "confirmation_chart", "title": "HFS 候选源的跨工况确认",
                 "subtitle": "n17/n21 与 Vd=0.05/1 V；负值表示相对 full 改善",
                 "type": "bar", "dataset": "confirmation_chart",
                 "sourceId": "confirmation", "valueFormat": "number",
                 "encodings": {
                     "x": {"field": "condition", "type": "nominal", "label": "器件与漏压"},
                     "y": {"field": "maximum_change_vs_full_dex", "type": "quantitative",
                           "label": "最大绝对对数误差变化, dex"},
                     "color": {"field": "variant_label", "type": "nominal", "label": "变体"},
                     "tooltip": [
                         {"field": "maximum_absolute_log10_ratio", "type": "quantitative",
                          "label": "最大误差, dex"},
                         {"field": "maximum_relative_error", "type": "quantitative",
                          "label": "最大相对误差"},
                     ],
                 }},
            ],
            "tables": [
                {"id": "confirmation_table", "title": "第二轮逐工况结果",
                 "subtitle": "12 条曲线，均为 51 个直接栅压点且禁止插值",
                 "dataset": "confirmation", "sourceId": "confirmation",
                 "defaultSort": {"field": "maximum_absolute_log10_ratio", "direction": "desc"},
                 "columns": [
                     {"field": "device", "label": "器件", "type": "text"},
                     {"field": "drain_voltage_V", "label": "Vd, V", "format": "number"},
                     {"field": "variant", "label": "变体", "type": "text"},
                     {"field": "maximum_absolute_log10_ratio", "label": "最大误差, dex",
                      "format": "number"},
                     {"field": "maximum_relative_error", "label": "最大相对误差",
                      "format": "number"},
                     {"field": "maximum_change_vs_full_dex", "label": "相对 full 变化, dex",
                      "format": "number"},
                 ]},
            ],
            "sources": [
                {"id": "first_round", "label": "M8-A 第一轮冻结汇总",
                 "path": first_source_path},
                {"id": "confirmation", "label": "M8-A 第二轮冻结汇总",
                 "path": confirmation_source_path},
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {"id": "technical_summary", "type": "markdown", "body":
                 "## 技术摘要\n\n本实验在不改变网格、器件结构与偏置格点的前提下，对 Sentaurus Device 与 Vela 同时关闭或简化同一组物理模型。"
                 f"第一轮识别出 HFS 为唯一降低 n23 最大误差的单项关闭项；{conclusion}"},
                {"id": "metrics", "type": "metric-strip",
                 "cardIds": ["full_error", "no_hfs_error", "hfs_gain", "no_bgn_error"]},
                {"id": "first_heading", "type": "markdown", "body":
                 "## 第一轮：模型开关定位\n\n八种成对配置覆盖 SRH、OldSlotboom、Enormal、HighFieldSaturation、PhuMob 与常数迁移率。负值代表误差下降。"},
                {"id": "first_chart", "type": "chart", "chartId": "first_round_chart"},
                {"id": "first_interpretation", "type": "markdown", "body":
                 "No HFS 的改善最大，PhuMob only 次之；关闭 SRH、Enormal 或采用常数迁移率均使误差增大。关闭 OldSlotboom 后误差剧增，说明该模型是当前一致性的必要条件，而不是可通过关闭消除的残差来源。"},
                {"id": "confirmation_heading", "type": "markdown", "body":
                 "## 第二轮：跨器件与漏压确认\n\n控制矩阵固定 n17/n21 两种 NWell 掺杂，并在 Vd=0.05 V 与 1 V 下比较 full、No HFS 和 PhuMob only。"},
                {"id": "confirmation_chart_block", "type": "chart",
                 "chartId": "confirmation_chart"},
                {"id": "confirmation_interpretation", "type": "markdown", "body":
                 f"{conclusion} 这一区分了“n23 单点改善”与“跨工况稳定贡献”，但仍不能把公开模型名称直接等同于 Sentaurus 内部某一条私有公式。"},
                {"id": "scope_method", "type": "markdown", "body":
                 "## 范围、数据与方法\n\n范围仅限 SDevice 器件仿真，不包含 SProcess。全部比较使用相同上游 TDR、相同接触、相同 Vg=0:0.05:2.5 V 的 51 点格点、绝对电流幅值和 1e-18 A/um 电流门槛；禁止插值。单项消融用于诊断，不单独构成因果证明。"},
                {"id": "confirmation_table_block", "type": "table",
                 "tableId": "confirmation_table"},
                {"id": "limitations", "type": "markdown", "body":
                 "## 局限性与不确定性\n\nSentaurus 的内部平滑、默认参数组合、驱动力计算与极限处理没有完全公开。成对关闭 HFS 只能说明 HFS 模型族相关差异对误差有贡献，无法进一步区分是驱动力离散、参数、平滑函数还是极限处理。No OldSlotboom 变体超出原数值一致性阈值，因此只用于敏感性诊断。"},
                {"id": "recommendations", "type": "markdown", "body":
                 "## 建议的下一步\n\n1. 若 No HFS 在四个控制工况均改善，优先细分 HFS 驱动力与速度饱和参数实验。\n2. 保持 OldSlotboom 与已修正 SRH 默认寿命，不将其作为下一轮关闭项。\n3. 对 HFS 进一步使用固定迁移率场或参数扫描时，继续采用双端成对配置与精确偏置格点。"},
                {"id": "questions", "type": "markdown", "body":
                 "## 后续问题\n\n需要进一步回答：误差改善主要来自 Vela 的 quasi-Fermi-gradient 驱动力离散，还是来自 Sentaurus 未公开的 HFS 参数与平滑；以及该贡献在更高漏压或不同沟道长度下是否保持。"},
            ],
        },
        "snapshot": {
            "version": 1, "generatedAt": generated, "status": "ready",
            "datasets": {"headline": headline, "first_round": first,
                         "first_round_zoom": first_zoom,
                         "confirmation": confirmation,
                         "confirmation_chart": confirmation_chart},
            "accessIssues": [],
        },
        "sources": [
            {"id": "first_round", "query": {
                "language": "sql", "engine": "duckdb",
                "sql": ("SELECT * FROM read_csv_auto('" + first_source_path
                        + "', header=true)"),
                "description": "Paired 51-point comparison of eight n23 model variants.",
                "executed_at": generated,
                "tables_used": [first_source_path],
                "metric_definitions": [
                    "maximum_absolute_log10_ratio is max abs(log10(abs(Id,Vela)/abs(Id,Sentaurus))) above 1e-18 A/um",
                    "maximum_change_vs_full_dex subtracts the full-physics maximum from the variant maximum",
                ],
                "filters": ["device n23", "Vd=0.05 V", "Vg=0:0.05:2.5 V",
                            "no interpolation"]}},
            {"id": "confirmation", "query": {
                "language": "sql", "engine": "duckdb",
                "sql": ("SELECT * FROM read_csv_auto('"
                        + confirmation_source_path + "', header=true)"),
                "description": "Paired confirmation for n17/n21, two drain biases, and three mobility variants.",
                "executed_at": generated,
                "tables_used": [confirmation_source_path],
                "metric_definitions": [
                    "maximum_absolute_log10_ratio is max abs(log10(abs(Id,Vela)/abs(Id,Sentaurus))) above 1e-18 A/um",
                    "negative maximum_change_vs_full_dex indicates improved paired agreement",
                ],
                "filters": ["devices n17 and n21", "Vd in 0.05 and 1 V",
                            "Vg=0:0.05:2.5 V", "no interpolation"]}},
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "built", "output": str(output),
                      "hfs_consistent_conditions": consistent_count}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
