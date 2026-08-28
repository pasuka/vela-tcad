#!/usr/bin/env python3
"""Build the canonical portable-report input for SimpleMOS M8-B."""

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
    / "simplemos_m8b_hfs_controls_evidence.json"
)
DEFAULT_OUTPUT = REPO / "docs/validation/reports/simplemos_m8b/artifact.json"
LABELS = {
    "explicit_gradqf": "Explicit GradQF (null)",
    "eparallel": "Eparallel",
    "qf_at_contacts": "QF at contacts",
    "no_parallel_boundary": "No parallel boundary",
    "refdens_efield_1e8": "RefDens E-field 1e8",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = []
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = dict(row)
            for key, value in row.items():
                if key in {"device", "variant", "dimension"}:
                    continue
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    pass
            rows.append(parsed)
        return rows


def source_path(evidence: dict[str, Any], key: str) -> Path:
    return REPO / evidence["artifacts"][key]["path"]


def average(values: list[float]) -> float:
    return sum(values) / len(values)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    evidence = read_json(args.evidence.resolve())
    report = read_json(source_path(evidence, "comparison_report"))
    convergence = read_json(source_path(evidence, "convergence_report"))
    summary_path = evidence["artifacts"]["summary"]["path"]
    summary = read_csv(REPO / summary_path)
    summaries = {item["variant"]: item
                 for item in report["variant_summaries"]}
    non_null = [item for item in report["variant_summaries"]
                if item["variant"] != "explicit_gradqf"]
    best = evidence["findings"]["best_mean_error_change_control"]
    best_summary = summaries[best]
    mean_best_change = average(best_summary["maximum_error_changes_dex"])
    material = evidence["findings"]["material_in_all_conditions"]
    improving = evidence["findings"]["improves_all_conditions"]
    null_response = summaries["explicit_gradqf"]["maximum_response_dex"]
    selected_damping = float(convergence["selected_probe"][
        "line_search_damping"])
    damping_robustness = max(
        (float(item["maximum_absolute_log10_current_difference_dex"])
         for item in convergence["robustness_pairs"]), default=0.0)

    effect_rows = []
    for row in summary:
        effect_rows.append({
            **row,
            "variant_label": LABELS[row["variant"]],
            "condition": f"{row['device']}, Vd={row['drain_voltage_V']:g} V",
        })
    best_response_rows = []
    for case in report["cases"]:
        if case["variant"] != best:
            continue
        for row in read_csv(REPO / case["response_csv"]):
            best_response_rows.append({
                **row,
                "device": case["device"],
                "drain_voltage_V": case["drain_voltage_V"],
                "variant": best,
                "condition": (f"{case['device']}, "
                              f"Vd={case['drain_voltage_V']:g} V"),
            })
    headline = [{
        "new_curves": evidence["execution"]["curves"],
        "direct_bias_points": evidence["execution"]["total_direct_bias_points"],
        "null_response_dex": null_response,
        "best_mean_change_dex": mean_best_change,
        "material_control_count": len(material),
        "stable_improving_control_count": len(improving),
        "selected_line_search_damping": selected_damping,
        "damping_robustness_dex": damping_robustness,
    }]
    generated = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z")
    title = "SimpleMOS M8-B HighFieldSaturation 控制项归因报告"
    best_label = LABELS[best]
    if improving:
        conclusion = (
            "跨四个工况稳定改善的控制项为 "
            + "、".join(LABELS[item] for item in improving) + "。")
    else:
        conclusion = "没有控制项在四个工况中都降低最大误差。"
    material_text = ("、".join(LABELS[item] for item in material)
                     if material else "无")
    summary_source = {
        "id": "control_summary", "label": "M8-B 冻结控制矩阵汇总",
        "path": summary_path,
    }
    response_source_path = evidence["artifacts"]["comparison_report"]["path"]
    response_glob = str(
        Path(response_source_path).parent
        / f"*_{best}_*_response.csv"
    ).replace("\\", "/")
    response_source = {
        "id": "control_response", "label": "M8-B 冻结逐点响应",
        "path": response_source_path,
    }
    convergence_source_path = evidence["artifacts"]["convergence_report"]["path"]
    convergence_source = {
        "id": "convergence_probe", "label": "M8-B Eparallel 收敛探针",
        "path": convergence_source_path,
    }
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report", "title": title,
            "description": "Sentaurus-only HFS control diagnosis against fixed Vela curves.",
            "generatedAt": generated,
            "cards": [
                {"id": "curve_count", "description": "新增 Sentaurus 单变量扫描曲线数。",
                 "dataset": "headline", "sourceId": "control_summary",
                 "metrics": [{"label": "新增曲线", "field": "new_curves",
                              "format": "number"}]},
                {"id": "point_count", "description": "全部为直接求解的精确栅压点。",
                 "dataset": "headline", "sourceId": "control_summary",
                 "metrics": [{"label": "直接偏置点", "field": "direct_bias_points",
                              "format": "number"}]},
                {"id": "null_response", "description": "显式 GradQF 相对默认设置的最大响应。",
                 "dataset": "headline", "sourceId": "control_summary",
                 "metrics": [{"label": "空对照响应", "field": "null_response_dex",
                              "format": "number", "unit": "dex"}]},
                {"id": "best_change", "description": f"{best_label} 在四个工况的最大误差变化均值。",
                 "dataset": "headline", "sourceId": "control_summary",
                 "metrics": [{"label": "最佳均值变化", "field": "best_mean_change_dex",
                              "format": "number", "unit": "dex", "signed": True}]},
                {"id": "selected_damping", "description": "使 n17/n21 高漏压 Eparallel 均收敛的统一线搜索阻尼。",
                 "dataset": "headline", "sourceId": "convergence_probe",
                 "metrics": [{"label": "统一阻尼", "field": "selected_line_search_damping",
                              "format": "number"}]},
                {"id": "damping_robustness", "description": "n21 两个成功阻尼解之间的最大逐点电流差。",
                 "dataset": "headline", "sourceId": "convergence_probe",
                 "metrics": [{"label": "阻尼敏感度", "field": "damping_robustness_dex",
                              "format": "number", "unit": "dex"}]},
            ],
            "charts": [
                {"id": "effect_chart", "title": "四个工况的最大误差变化",
                 "subtitle": "负值表示相对 Sentaurus 默认 HFS 配置改善",
                 "type": "bar", "dataset": "effect_rows",
                 "sourceId": "control_summary", "valueFormat": "number",
                 "encodings": {
                     "x": {"field": "condition", "type": "nominal",
                           "label": "器件与漏压"},
                     "y": {"field": "maximum_error_change_dex",
                           "type": "quantitative",
                           "label": "最大绝对对数误差变化, dex"},
                     "color": {"field": "variant_label", "type": "nominal",
                               "label": "控制项"},
                     "tooltip": [
                         {"field": "maximum_sentaurus_response_dex",
                          "type": "quantitative", "label": "Sentaurus 最大响应, dex"},
                         {"field": "maximum_absolute_log10_ratio_above_floor",
                          "type": "quantitative", "label": "控制后最大误差, dex"},
                     ],
                 }},
                {"id": "response_chart", "title": f"{best_label} 的 Sentaurus 响应",
                 "subtitle": "相对默认 HFS 的逐栅压电流对数变化",
                 "type": "line", "dataset": "best_response",
                 "sourceId": "control_response", "valueFormat": "number",
                 "encodings": {
                     "x": {"field": "gate_voltage_V", "type": "quantitative",
                           "label": "栅压, V"},
                     "y": {"field": "sentaurus_log10_response",
                           "type": "quantitative",
                           "label": "log10(|Id,control|/|Id,default|), dex"},
                     "color": {"field": "condition", "type": "nominal",
                               "label": "器件与漏压"},
                     "tooltip": [
                         {"field": "controlled_absolute_log10_error",
                          "type": "quantitative", "label": "控制后绝对误差, dex"},
                         {"field": "controlled_sentaurus_current_A_per_um",
                          "type": "quantitative", "label": "Sentaurus Id, A/um"},
                     ],
                 }},
            ],
            "tables": [
                {"id": "case_table", "title": "20 组逐工况诊断指标",
                 "subtitle": "固定 Vela；只改变一个 Sentaurus HFS 控制项",
                 "dataset": "effect_rows", "sourceId": "control_summary",
                 "defaultSort": {"field": "maximum_sentaurus_response_dex",
                                 "direction": "desc"},
                 "columns": [
                     {"field": "device", "label": "器件", "type": "text"},
                     {"field": "drain_voltage_V", "label": "Vd, V", "format": "number"},
                     {"field": "variant_label", "label": "控制项", "type": "text"},
                     {"field": "maximum_sentaurus_response_dex", "label": "最大响应, dex",
                      "format": "number"},
                     {"field": "maximum_error_change_dex", "label": "最大误差变化, dex",
                      "format": "number"},
                     {"field": "maximum_absolute_log10_ratio_above_floor",
                      "label": "控制后最大误差, dex", "format": "number"},
                 ]},
            ],
            "sources": [summary_source, response_source, convergence_source],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {"id": "summary", "type": "markdown",
                 "body": "## 技术摘要\n\n"
                 f"空对照最大响应为 {null_response:.3e} dex。{conclusion} "
                 f"按四工况最大误差变化均值，变化最有利的是 {best_label}（{mean_best_change:+.6f} dex），但其方向并不一致，不能作为统一修正。"
                 f"在全部工况达到 1e-4 dex 响应门槛的控制项：{material_text}。"
                 f"两个高漏压 Eparallel 工况使用统一 LineSearchDamping={selected_damping:g} 恢复收敛，"
                 f"成功阻尼解间最大差为 {damping_robustness:.3e} dex。"},
                {"id": "metrics", "type": "metric-strip",
                 "cardIds": ["curve_count", "point_count", "null_response",
                             "best_change", "selected_damping", "damping_robustness"]},
                {"id": "scope", "type": "markdown", "body":
                 "## 范围、数据与指标定义\n\n范围仅包含 SDevice 器件仿真，不包含 SProcess。每个控制在 n17/n21、Vd=0.05/1 V 下与复用的 M8-A 默认 Sentaurus 曲线比较，同时固定 Vela M8-A full 曲线。Sentaurus 响应定义为 log10(|Id,control|/|Id,default|)；最大误差变化定义为控制后与默认配置的最大 |log10(|Id,Vela|/|Id,Sentaurus|)| 之差。电流门槛为 1e-18 A/um，禁止偏置插值。"},
                {"id": "design", "type": "markdown", "body":
                 "## 实验设计与数值稳健性\n\n矩阵覆盖五类 Sentaurus 配置：显式 GradQuasiFermi 空对照、Eparallel 驱动力、接触处准费米梯度替换、关闭边界层界面平行处理、以及 1e8 cm^-3 的 HFS 电场参考密度。每条曲线固定 Vg=0:0.05:2.5 V 共 51 点；上游 TDR、接触和物理模型不变。两个 Vd=1 V 的 Eparallel 默认连续扫描触及 MinStep；六组阻尼探针表明 0.5 均失败、0.1 仅 n21 成功、0.01 在两器件均成功。正式重试只提高迭代上限、收紧 MinStep 并对栅压 Coupled 步使用统一 0.01 阻尼。"},
                {"id": "robustness", "type": "markdown",
                 "sourceId": "convergence_probe", "body":
                 f"## 阻尼改变收敛路径，但没有实质改变解\n\nn21 在 LineSearchDamping=0.1 与 0.01 下均输出 51 个精确栅压点，两条电流曲线最大差仅 {damping_robustness:.6f} dex。该差异远小于约 0.1 dex 的跨求解器残差，因此阻尼用于恢复数值收敛，不作为物理模型控制项参与归因。"},
                {"id": "effect_heading", "type": "markdown", "body":
                 "## 哪个控制项稳定改变剩余误差\n\n图中同时报告方向与跨工况稳定性；单一工况改善不能升级为稳定归因。"},
                {"id": "effect", "type": "chart", "chartId": "effect_chart"},
                {"id": "response_heading", "type": "markdown", "body":
                 f"## 平均误差变化最有利控制的偏置位置\n\n{best_label} 只是平均指标最有利、并非跨工况稳定改善；逐点响应用于判断变化集中在深关断、弱反型还是强反型，避免只用一个最大值掩盖偏置依赖。"},
                {"id": "response", "type": "chart", "chartId": "response_chart"},
                {"id": "table_heading", "type": "markdown", "body":
                 "## 逐工况结果\n\n最大响应衡量 Sentaurus 选项本身的可辨识度；最大误差变化衡量它对固定 Vela 比较的方向。"},
                {"id": "table", "type": "table", "tableId": "case_table"},
                {"id": "limitations", "type": "markdown", "body":
                 "## 局限性与不确定性\n\n这些控制能把残差关联到驱动力、接触替换、边界层处理或参考密度插值，但不能反演 Sentaurus 私有内部公式。参考密度只测试 1e8 cm^-3 一个值；控制之间可能存在交互，本轮采用的是单变量而非全因子设计。亚飞安点仍受电流门槛与非线性收敛容差影响。"},
                {"id": "recommendations", "type": "markdown", "body":
                 "## 建议的下一步\n\n1. 不把 Eparallel 或 RefDens_ElectricField=1e8 作为 Vela 的误差修正：两者在四个工况中都使最大误差增大。\n2. 下一轮优先排查默认 GradQuasiFermi HFS 的参数、平滑与极限处理，并可对参考密度做更细的对数扫描，而不是只测试 1e8 cm^-3。\n3. 接触处准费米梯度替换方向混合，关闭边界层界面平行处理低于响应门槛，可降低排查优先级。\n4. 保留显式 GradQF 空对照和 LineSearchDamping=0.01 的高漏压恢复策略，分别作为模型漂移守卫和数值收敛策略。"},
                {"id": "questions", "type": "markdown", "body":
                 "## 后续问题\n\n控制效应是否随沟道长度、漏压继续升高或网格细化而保持；以及多个控制同时开启时是否呈可加性，仍需后续实验回答。"},
            ],
        },
        "snapshot": {
            "version": 1, "generatedAt": generated, "status": "ready",
            "datasets": {"headline": headline, "effect_rows": effect_rows,
                         "best_response": best_response_rows},
            "accessIssues": [],
        },
        "sources": [
            {"id": "control_summary", "query": {
                "language": "sql", "engine": "duckdb",
                "sql": "SELECT * FROM read_csv_auto('" + summary_path
                       + "', header=true)",
                "description": "Twenty fixed-Vela HFS-control comparisons.",
                "executed_at": generated, "tables_used": [summary_path],
                "metric_definitions": [
                    "maximum_error_change_dex subtracts the reused M8-A default-HFS maximum absolute log10 error from the controlled case",
                    "maximum_sentaurus_response_dex is max abs(log10(abs(Id_control)/abs(Id_default))) over 51 exact gate points",
                ],
                "filters": ["devices n17 and n21", "Vd in 0.05 and 1 V",
                            "Vg=0:0.05:2.5 V", "current floor 1e-18 A/um",
                            "no interpolation"]}},
            {"id": "control_response", "query": {
                "language": "sql", "engine": "duckdb",
                "sql": "SELECT * FROM read_csv_auto('" + response_glob
                       + "', header=true, filename=true)",
                "description": "Exact-bias Sentaurus response and fixed-Vela residual rows.",
                "executed_at": generated, "tables_used": [response_source_path],
                "metric_definitions": [
                    "sentaurus_log10_response = log10(abs(Id_control)/abs(Id_default))",
                    "controlled_absolute_log10_error = abs(log10(abs(Id_Vela)/abs(Id_control)))",
                ],
                "filters": [f"variant {best}", "all four device/drain conditions",
                            "51 direct gate points per condition"]}},
            {"id": "convergence_probe", "query": {
                "language": "sql", "engine": "duckdb",
                "sql": "SELECT * FROM read_json_auto('" + convergence_source_path.replace("\\", "/") + "')",
                "description": "Six Eparallel damping probes for the two Vd=1 V failures.",
                "executed_at": generated,
                "tables_used": [convergence_source_path],
                "metric_definitions": [
                    "success requires exit code zero, no failure marker, and exactly 51 direct gate points",
                    "damping robustness is the maximum pointwise absolute log10 current difference between two successful n21 damping probes",
                ],
                "filters": ["devices n17 and n21", "Vd=1 V",
                            "LineSearchDamping in 0.5, 0.1, and 0.01"]}},
        ],
    }
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "built", "output": str(output),
                      "best_control": best,
                      "stable_improving_controls": improving}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
