#!/usr/bin/env python3
"""Build the bounded report artifact for the Genius BJT recovery audit."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--result", type=Path,
        default=FIXTURE / "reports" / "hole_current_recovery_topology.json",
    )
    parser.add_argument(
        "--output", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "hole_current_recovery_topology" / "report_artifact.json",
    )
    args = parser.parse_args()
    result = json.loads(args.result.read_text(encoding="utf-8"))
    groups = {item["name"]: item for item in result["groups"]}

    comparisons = []
    for key, label in (
        ("dominant_current", "主电流区"),
        ("weak_main_current", "弱电流区"),
        ("original_tail", "原始超限尾部"),
    ):
        item = groups[key]
        for method, metric_key, failure_key in (
            ("节点直接拟合", "sdevice_vs_direct_node_fit", "direct_node_fit"),
            ("先单元后节点", "sdevice_vs_cell_first_area_projection", "cell_first_area_projection"),
        ):
            metric = item[metric_key]
            comparisons.append(
                {
                    "population": label,
                    "population_key": key,
                    "recovery_method": method,
                    "p95_error_decade": metric["magnitude_abs_dex_p95"],
                    "median_error_decade": metric["magnitude_abs_dex_median"],
                    "normalized_vector_rmse": metric["normalized_vector_rmse"],
                    "cosine_similarity": metric["cosine_similarity"],
                    "node_count": item["node_count"],
                    "nodes_above_0p5_decade": item["nodes_above_0p5_decade"][failure_key],
                }
            )

    topology = []
    for key, label in (
        ("dominant_current", "主电流区"),
        ("weak_main_current", "弱电流区"),
        ("original_tail", "原始超限尾部"),
    ):
        item = groups[key]
        topology.append(
            {
                "population": label,
                "node_count": item["node_count"],
                "fit_residual_median": item["direct_fit_normalized_residual"]["median"],
                "fit_residual_p95": item["direct_fit_normalized_residual"]["p95"],
                "cell_dispersion_median": item["incident_cell_vector_dispersion"]["median"],
                "cell_cancellation_median": item["incident_cell_cancellation_factor"]["median"],
                "condition_number_median": item["nodal_fit_condition_number"]["median"],
            }
        )

    title = "Genius NPN BJT 弱空穴电流节点恢复定位"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    comparison_csv = args.output.parent / "report_recovery_comparison.csv"
    topology_csv = args.output.parent / "report_topology_summary.csv"
    for output_path, records in (
        (comparison_csv, comparisons),
        (topology_csv, topology),
    ):
        with output_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    result_source = {
        "id": "topology_audit",
        "label": "弱空穴电流恢复拓扑审计",
        "path": "reference_tcad/genius_bjt_sentaurus2022/reports/hole_current_recovery_topology.json",
    }
    comparison_source = {
        "id": "recovery_comparison_source",
        "label": "节点恢复误差分组结果",
        "path": comparison_csv.relative_to(REPO).as_posix(),
        "query": {
            "engine": "duckdb",
            "language": "sql",
            "description": "读取相同 SG 边通量的直接节点与单元优先恢复误差。",
            "sql": "SELECT * FROM read_csv_auto('report_recovery_comparison.csv', header=true) ORDER BY population, recovery_method",
            "tables_used": ["report_recovery_comparison.csv"],
            "filters": ["VBE=0.70 V", "VCE=3.00 V", "coarse common mesh"],
            "metric_definitions": ["p95_error_decade is the 95th percentile of absolute log10 current-magnitude error within each population"],
        },
    }
    topology_source = {
        "id": "topology_summary_source",
        "label": "节点邻域拓扑汇总",
        "path": topology_csv.relative_to(REPO).as_posix(),
        "query": {
            "engine": "duckdb",
            "language": "sql",
            "description": "读取节点群体的恢复拓扑统计。",
            "sql": "SELECT * FROM read_csv_auto('report_topology_summary.csv', header=true) ORDER BY cell_dispersion_median DESC",
            "tables_used": ["report_topology_summary.csv"],
            "filters": ["VBE=0.70 V", "VCE=3.00 V", "coarse common mesh"],
            "metric_definitions": ["cell_dispersion_median is the median area-weighted RMS spread of incident cell vectors normalized by the area-projected node-vector magnitude"],
        },
    }
    sources = [result_source, comparison_source, topology_source]
    tail = groups["original_tail"]
    main = groups["main_current"]
    generated = datetime.now(timezone.utc).isoformat()
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "generatedAt": generated,
            "cards": [],
            "charts": [
                {
                    "id": "recovery_error_chart",
                    "title": "不同节点恢复顺序的空穴电流 P95 误差",
                    "description": "VBE=0.70 V、VCE=3.00 V 粗网格；误差单位为 decade，0.5 decade 为当前诊断门槛。",
                    "type": "bar",
                    "dataset": "recovery_comparison",
                    "sourceId": "recovery_comparison_source",
                    "encodings": {
                        "x": {"field": "population", "type": "nominal", "title": "节点群体"},
                        "y": {"field": "p95_error_decade", "type": "quantitative", "title": "P95 绝对对数误差（decade）"},
                        "color": {"field": "recovery_method", "type": "nominal"},
                    },
                }
            ],
            "tables": [
                {
                    "id": "topology_table",
                    "title": "节点邻域拓扑诊断",
                    "description": "同一 SG 边通量在不同节点群体中的拟合残差、单元方向离散度、抵消因子和条件数。",
                    "dataset": "topology_summary",
                    "sourceId": "topology_summary_source",
                    "defaultSort": {"field": "cell_dispersion_median", "direction": "desc"},
                    "columns": [
                        {"field": "population", "label": "节点群体", "format": "text"},
                        {"field": "node_count", "label": "节点数", "format": "number"},
                        {"field": "fit_residual_median", "label": "拟合残差中位数", "format": "number"},
                        {"field": "cell_dispersion_median", "label": "单元矢量离散度", "format": "number"},
                        {"field": "cell_cancellation_median", "label": "幅值抵消因子", "format": "number"},
                        {"field": "condition_number_median", "label": "条件数", "format": "number"},
                    ],
                }
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {
                    "id": "summary", "type": "markdown", "sourceId": "topology_audit",
                    "body": (
                        "## 技术摘要\n\n"
                        "粗网格空穴节点电流的 0.828-decade P95 失败已经定位到边通量向节点矢量的恢复支撑与运算顺序，而不是器件状态、迁移率或 SG 通量本身。"
                        f"用同一边通量改为先单元后节点后，完整掩膜 P95 降至 `{main['sdevice_vs_cell_first_area_projection']['magnitude_abs_dex_p95']:.4f} decade`，"
                        "180 个超限节点全部消失。该候选目前仅用于诊断，尚未替换生产默认值。"
                    ),
                },
                {
                    "id": "comparison_text", "type": "markdown", "sourceId": "topology_audit",
                    "body": (
                        "## 恢复顺序解释了弱电流尾部\n\n"
                        f"生产公式由保存的 SG 边数据重放到 `{result['production_replay']['maximum_absolute_vector_difference_A_per_cm2']:.2e} A/cm²`。"
                        f"原始尾部 180 个节点的 P95 从 `{tail['sdevice_vs_direct_node_fit']['magnitude_abs_dex_p95']:.4f}` 降至 "
                        f"`{tail['sdevice_vs_cell_first_area_projection']['magnitude_abs_dex_p95']:.4f} decade`；"
                        "主电流区两种恢复都已经低于 0.05 decade，因此收益集中在弱电流基区。"
                    ),
                },
                {"id": "comparison_chart", "type": "chart", "chartId": "recovery_error_chart"},
                {
                    "id": "definitions", "type": "markdown",
                    "body": (
                        "## 范围与指标定义\n\n"
                        "对象为 VBE=0.70 V、VCE=3.00 V、5611 节点/10940 三角形的共同粗网格。"
                        "主掩膜保留 SDevice 电流幅值不低于峰值 1e-6 的节点；主电流区阈值为峰值 1e-2；"
                        "原始尾部指直接节点拟合误差超过 0.5 decade 的 180 个节点。P95 为群体绝对 log10 幅值误差的第 95 百分位。"
                    ),
                },
                {
                    "id": "topology_text", "type": "markdown", "sourceId": "topology_audit",
                    "body": (
                        "## 拓扑证据排除了病态矩阵和简单抵消\n\n"
                        "尾部拟合条件数中位数为 2.0，幅值抵消因子中位数为 1.00009，均不足以解释约 22 倍的幅值差。"
                        "单元电流矢量离散度中位数则由主电流区的 0.108 增至尾部的 1.847，表明弱电流顶点邻域明显偏离单一平滑矢量假设。"
                    ),
                },
                {"id": "topology", "type": "table", "tableId": "topology_table"},
                {
                    "id": "method", "type": "markdown",
                    "body": (
                        "## 方法与限制\n\n"
                        "直接节点路径仅拟合目标节点的入射边；单元优先路径先用每个三角形的三条边重构单元矢量，再按面积投影到节点，因此包含目标节点对边的信息。"
                        "同一状态、迁移率和 SG 边通量使该 A/B 只改变表示支撑。SDevice 未公开精确的单元到顶点权重，因此当前结论能锁定差异类别，但不能宣称复刻了其专有公式。"
                    ),
                },
                {
                    "id": "next", "type": "markdown",
                    "body": (
                        "## 建议的下一步\n\n"
                        "将单元优先恢复实现为默认关闭的诊断候选，在 0、1、2、3 V 及局部加密网格上回归。"
                        "若端口守恒保持不变且节点场门槛稳定通过，再讨论是否将其设为绘图默认值；生产有限体积通量与求解残差无需改变。"
                    ),
                },
                {
                    "id": "questions", "type": "markdown",
                    "body": (
                        "## 尚待确认\n\n"
                        "需要进一步确认该恢复方式在其他器件、非结构化网格和材料界面附近是否仍保持优势，以及是否需要区域内投影以避免跨材料平滑。"
                    ),
                },
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated,
            "status": "ready",
            "datasets": {
                "recovery_comparison": comparisons,
                "topology_summary": topology,
            },
        },
        "sources": sources,
    }
    args.output.write_text(
        json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
