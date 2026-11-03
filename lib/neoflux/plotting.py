"""
SITC-quality figure generation for NeoFlux.

Produces publication-ready figures suitable for SITC 2026 oral
presentation and poster. Uses matplotlib with a custom style for
consistent, professional appearance.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Non-interactive backend for pipeline use
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger("neoflux")

# NeoFlux color palette -- designed for accessibility and print
NEOFLUX_COLORS = {
    "SNV": "#4C72B0",
    "insertion": "#55A868",
    "deletion": "#C44E52",
    "frameshift": "#8172B3",
    "fusion": "#CCB974",
    "splice": "#64B5CD",
    "non_canonical_orf": "#E5AE38",
    "primary": "#2C3E50",
    "secondary": "#E74C3C",
    "accent": "#3498DB",
    "positive": "#27AE60",
    "negative": "#95A5A6",
}

# Publication figure parameters
FIGURE_DPI = 300
FONT_SIZE_TITLE = 14
FONT_SIZE_LABEL = 12
FONT_SIZE_TICK = 10


def set_neoflux_style():
    """Apply NeoFlux publication style to matplotlib."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": FONT_SIZE_TICK,
        "axes.titlesize": FONT_SIZE_TITLE,
        "axes.labelsize": FONT_SIZE_LABEL,
        "xtick.labelsize": FONT_SIZE_TICK,
        "ytick.labelsize": FONT_SIZE_TICK,
        "legend.fontsize": FONT_SIZE_TICK,
        "figure.dpi": FIGURE_DPI,
        "savefig.dpi": FIGURE_DPI,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    sns.set_palette(list(NEOFLUX_COLORS.values())[:7])


def plot_neoantigen_burden(
    burden_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Neoantigen Burden by Variant Type",
) -> Path:
    """
    Stacked bar chart of neoantigen burden per patient by variant type.

    Args:
        burden_df: DataFrame from CohortAnalyzer.compute_neoantigen_burden().
        output_path: Path to save the figure.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    set_neoflux_style()

    variant_cols = [
        "SNV", "insertion", "deletion", "frameshift",
        "fusion", "splice", "non_canonical_orf",
    ]
    # Only plot columns that exist and have non-zero values
    plot_cols = [c for c in variant_cols if c in burden_df.columns]
    colors = [NEOFLUX_COLORS.get(c, "#999999") for c in plot_cols]

    fig, ax = plt.subplots(figsize=(max(8, len(burden_df) * 0.5), 5))

    burden_df_sorted = burden_df.sort_values("total", ascending=False)
    x = np.arange(len(burden_df_sorted))
    bottom = np.zeros(len(burden_df_sorted))

    for col, color in zip(plot_cols, colors):
        values = burden_df_sorted[col].values.astype(float)
        ax.bar(x, values, bottom=bottom, color=color, label=col, width=0.7)
        bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(
        burden_df_sorted["sample_id"], rotation=45, ha="right"
    )
    ax.set_xlabel("Patient")
    ax.set_ylabel("Number of Neoantigen Candidates")
    ax.set_title(title)
    ax.legend(title="Variant Type", bbox_to_anchor=(1.02, 1), loc="upper left")

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved neoantigen burden plot: %s", output_path)
    return output_path


def plot_score_distribution(
    candidates_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "NeoFlux Score Distribution",
) -> Path:
    """
    Violin plot of NeoFlux scores stratified by MS detection status.

    For SITC: demonstrates that MS-detected neoantigens receive higher
    ML scores, validating the model.

    Args:
        candidates_df: Cohort-wide candidate DataFrame.
        output_path: Path to save the figure.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    set_neoflux_style()

    fig, ax = plt.subplots(figsize=(6, 5))

    if "ms_detected" in candidates_df.columns and "neoflux_score" in candidates_df.columns:
        plot_df = candidates_df.dropna(subset=["neoflux_score", "ms_detected"])
        plot_df["MS Detection"] = plot_df["ms_detected"].map(
            {True: "Detected", False: "Not Detected"}
        )

        if len(plot_df) > 0:
            sns.violinplot(
                data=plot_df,
                x="MS Detection",
                y="neoflux_score",
                palette=[NEOFLUX_COLORS["positive"], NEOFLUX_COLORS["negative"]],
                ax=ax,
                inner="box",
                cut=0,
            )
    else:
        # If no MS data, just plot score distribution
        scores = candidates_df.get("neoflux_score", pd.Series(dtype=float)).dropna()
        if len(scores) > 0:
            sns.histplot(scores, ax=ax, bins=50, color=NEOFLUX_COLORS["primary"])

    ax.set_ylabel("NeoFlux Prioritization Score")
    ax.set_title(title)

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved score distribution plot: %s", output_path)
    return output_path


def plot_score_vs_ms_correlation(
    correlation_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "NeoFlux Score vs MS Detection Rate",
) -> Path:
    """
    Bar chart showing MS detection rate by NeoFlux score decile.

    For SITC: shows positive correlation between model confidence and
    empirical evidence of presentation.

    Args:
        correlation_df: DataFrame from CohortAnalyzer.compute_score_vs_ms_correlation().
        output_path: Path to save the figure.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    set_neoflux_style()

    fig, ax1 = plt.subplots(figsize=(8, 5))

    x = np.arange(len(correlation_df))

    # Bar chart: detection rate
    ax1.bar(
        x,
        correlation_df["detection_rate"],
        color=NEOFLUX_COLORS["accent"],
        alpha=0.8,
        width=0.6,
        label="MS Detection Rate",
    )
    ax1.set_ylabel("MS Detection Rate", color=NEOFLUX_COLORS["accent"])
    ax1.set_ylim(0, 1)
    ax1.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))

    # Secondary axis: number of candidates
    ax2 = ax1.twinx()
    ax2.plot(
        x,
        correlation_df["n_candidates"],
        color=NEOFLUX_COLORS["primary"],
        marker="o",
        linewidth=2,
        label="N Candidates",
    )
    ax2.set_ylabel("Number of Candidates", color=NEOFLUX_COLORS["primary"])

    ax1.set_xticks(x)
    ax1.set_xticklabels(
        correlation_df["score_bin"], rotation=45, ha="right"
    )
    ax1.set_xlabel("NeoFlux Score Bin")
    ax1.set_title(title)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved score vs MS correlation plot: %s", output_path)
    return output_path


def plot_hla_loh_heatmap(
    loh_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "HLA LOH Impact on Neoantigen Presentation",
) -> Path:
    """
    Heatmap showing fraction of neoantigen candidates lost to HLA LOH.

    Args:
        loh_df: DataFrame from CohortAnalyzer.compute_hla_loh_impact().
        output_path: Path to save the figure.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    set_neoflux_style()

    fig, ax = plt.subplots(figsize=(8, 4))

    patients = loh_df["sample_id"].values
    fractions = loh_df["fraction_lost"].values

    # Horizontal bar chart
    colors = [
        NEOFLUX_COLORS["secondary"] if f > 0.1 else NEOFLUX_COLORS["negative"]
        for f in fractions
    ]

    y_pos = np.arange(len(patients))
    ax.barh(y_pos, fractions, color=colors, height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(patients)
    ax.set_xlabel("Fraction of Candidates Affected by HLA LOH")
    ax.set_title(title)
    ax.set_xlim(0, max(fractions.max() * 1.2, 0.3) if len(fractions) > 0 else 0.3)
    ax.xaxis.set_major_formatter(ticker.PercentFormatter(1.0))
    ax.invert_yaxis()

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved HLA LOH heatmap: %s", output_path)
    return output_path


def plot_feature_importance(
    importance_df: pd.DataFrame,
    output_path: str | Path,
    top_n: int = 15,
    title: str = "NeoFlux Model Feature Importance",
) -> Path:
    """
    Horizontal bar chart of top ML model feature importances.

    Args:
        importance_df: DataFrame with columns [feature, importance].
        output_path: Path to save the figure.
        top_n: Number of top features to display.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    set_neoflux_style()

    top_df = importance_df.head(top_n).sort_values("importance")

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.barh(
        top_df["feature"],
        top_df["importance"],
        color=NEOFLUX_COLORS["accent"],
        height=0.6,
    )
    ax.set_xlabel("Feature Importance (Gain)")
    ax.set_title(title)

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved feature importance plot: %s", output_path)
    return output_path


def plot_roc_pr_curves(
    y_true: np.ndarray,
    y_scores: dict[str, np.ndarray],
    output_path: str | Path,
    title: str = "NeoFlux Model Performance",
) -> Path:
    """
    Side-by-side ROC and Precision-Recall curves for model benchmarking.

    Args:
        y_true: True binary labels.
        y_scores: Dict mapping model name to predicted scores.
        output_path: Path to save the figure.
        title: Figure title.

    Returns:
        Path to the saved figure.
    """
    from sklearn.metrics import (
        average_precision_score,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )

    set_neoflux_style()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    colors = list(NEOFLUX_COLORS.values())

    for i, (name, scores) in enumerate(y_scores.items()):
        color = colors[i % len(colors)]

        # ROC curve
        fpr, tpr, _ = roc_curve(y_true, scores)
        auroc = roc_auc_score(y_true, scores)
        ax1.plot(fpr, tpr, color=color, linewidth=2,
                 label=f"{name} (AUROC={auroc:.3f})")

        # PR curve
        precision, recall, _ = precision_recall_curve(y_true, scores)
        auprc = average_precision_score(y_true, scores)
        ax2.plot(recall, precision, color=color, linewidth=2,
                 label=f"{name} (AUPRC={auprc:.3f})")

    # ROC plot formatting
    ax1.plot([0, 1], [0, 1], "k--", alpha=0.3, linewidth=1)
    ax1.set_xlabel("False Positive Rate")
    ax1.set_ylabel("True Positive Rate")
    ax1.set_title("ROC Curve")
    ax1.legend(loc="lower right")

    # PR plot formatting
    baseline = y_true.sum() / len(y_true)
    ax2.axhline(y=baseline, color="k", linestyle="--", alpha=0.3, linewidth=1)
    ax2.set_xlabel("Recall")
    ax2.set_ylabel("Precision")
    ax2.set_title("Precision-Recall Curve")
    ax2.legend(loc="upper right")

    fig.suptitle(title, fontsize=FONT_SIZE_TITLE + 2, y=1.02)

    output_path = Path(output_path)
    fig.savefig(output_path)
    plt.close(fig)
    logger.info("Saved ROC/PR curves: %s", output_path)
    return output_path
