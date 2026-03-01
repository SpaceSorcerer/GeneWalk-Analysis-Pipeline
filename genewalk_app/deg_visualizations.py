"""DESeq2 / Differential Expression visualizations.

Standard plots for exploring differential expression results before
and alongside GeneWalk/GSEA analyses.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# MA plot
# ---------------------------------------------------------------------------

def ma_plot(
    deg: pd.DataFrame,
    basemean_col: str = "baseMean",
    log2fc_col: str = "log2fc",
    padj_col: str = "padj",
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> go.Figure:
    """MA plot: log2 fold-change vs mean expression (baseMean).

    The classic DESeq2 diagnostic.  X-axis is log10(baseMean), Y-axis is
    log2FC.  Points are colored by significance.
    """
    if basemean_col not in deg.columns or log2fc_col not in deg.columns:
        return go.Figure().update_layout(
            title="MA Plot requires baseMean and log2FC columns"
        )

    df = deg.dropna(subset=[basemean_col, log2fc_col]).copy()
    df["log10_basemean"] = np.log10(df[basemean_col].clip(lower=0.1))
    df["category"] = "Not significant"
    if padj_col in df.columns:
        df.loc[
            (df[padj_col] <= padj_threshold) & (df[log2fc_col] >= fc_threshold),
            "category",
        ] = "Up-regulated"
        df.loc[
            (df[padj_col] <= padj_threshold) & (df[log2fc_col] <= -fc_threshold),
            "category",
        ] = "Down-regulated"

    fig = px.scatter(
        df,
        x="log10_basemean",
        y=log2fc_col,
        color="category",
        color_discrete_map={
            "Up-regulated": "#D9534F",
            "Down-regulated": "#4A90D9",
            "Not significant": "#B0BEC5",
        },
        hover_data=["gene"] if "gene" in df.columns else None,
        labels={
            "log10_basemean": "log10(Mean Expression)",
            log2fc_col: "log2 Fold Change",
        },
        opacity=0.6,
    )
    fig.add_hline(y=0, line_color="#333", line_width=1)
    fig.add_hline(y=fc_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.add_hline(y=-fc_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.update_layout(
        title="MA Plot (Expression vs Fold Change)",
        template="plotly_white",
        height=500,
    )
    return fig


# ---------------------------------------------------------------------------
# P-value distribution
# ---------------------------------------------------------------------------

def pvalue_histogram(
    deg: pd.DataFrame,
    pval_col: str = "pvalue",
) -> go.Figure:
    """Histogram of raw p-values.

    A well-calibrated DE analysis produces a uniform distribution with
    a spike near zero (true positives). Anti-conservative (all near zero)
    or bimodal shapes suggest problems.
    """
    if pval_col not in deg.columns:
        return go.Figure().update_layout(title=f"Column '{pval_col}' not found")

    df = deg.dropna(subset=[pval_col])

    fig = px.histogram(
        df,
        x=pval_col,
        nbins=50,
        color_discrete_sequence=["#4A90D9"],
        labels={pval_col: "Raw P-value"},
    )
    fig.update_layout(
        title="Raw P-value Distribution",
        template="plotly_white",
        height=400,
        xaxis_title="P-value",
        yaxis_title="Count",
    )
    return fig


def padj_histogram(
    deg: pd.DataFrame,
    padj_col: str = "padj",
) -> go.Figure:
    """Histogram of adjusted p-values."""
    if padj_col not in deg.columns:
        return go.Figure().update_layout(title=f"Column '{padj_col}' not found")

    df = deg.dropna(subset=[padj_col])

    fig = px.histogram(
        df,
        x=padj_col,
        nbins=50,
        color_discrete_sequence=["#5CB85C"],
        labels={padj_col: "Adjusted P-value (FDR)"},
    )
    fig.update_layout(
        title="Adjusted P-value Distribution",
        template="plotly_white",
        height=400,
        xaxis_title="Adjusted P-value",
        yaxis_title="Count",
    )
    return fig


# ---------------------------------------------------------------------------
# Top genes ranked bar chart
# ---------------------------------------------------------------------------

def top_genes_bar(
    deg: pd.DataFrame,
    log2fc_col: str = "log2fc",
    padj_col: str = "padj",
    padj_threshold: float = 0.05,
    top_n: int = 30,
) -> go.Figure:
    """Horizontal bar chart of the top N most significant genes by |log2FC|."""
    if log2fc_col not in deg.columns:
        return go.Figure().update_layout(title="Missing log2FC column")

    df = deg.copy()
    if padj_col in df.columns:
        df = df[df[padj_col] <= padj_threshold]
    if df.empty:
        return go.Figure().update_layout(title="No significant genes at this threshold")

    df["abs_log2fc"] = df[log2fc_col].abs()
    df = df.sort_values("abs_log2fc", ascending=False).head(top_n)
    df = df.sort_values(log2fc_col, ascending=True)
    df["direction"] = df[log2fc_col].apply(
        lambda x: "Up-regulated" if x > 0 else "Down-regulated"
    )

    fig = px.bar(
        df,
        x=log2fc_col,
        y="gene" if "gene" in df.columns else df.index,
        color="direction",
        color_discrete_map={
            "Up-regulated": "#D9534F",
            "Down-regulated": "#4A90D9",
        },
        orientation="h",
        hover_data=[c for c in [padj_col, "baseMean"] if c in df.columns],
        labels={log2fc_col: "log2 Fold Change", "gene": ""},
    )
    fig.update_layout(
        title=f"Top {min(top_n, len(df))} Differentially Expressed Genes",
        template="plotly_white",
        height=max(450, len(df) * 22),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


# ---------------------------------------------------------------------------
# log2FC vs Wald statistic scatter
# ---------------------------------------------------------------------------

def lfc_vs_stat_scatter(
    deg: pd.DataFrame,
    log2fc_col: str = "log2fc",
    stat_col: str = "stat",
    padj_col: str = "padj",
    padj_threshold: float = 0.05,
) -> go.Figure:
    """Scatter of log2FC vs Wald test statistic (DESeq2 'stat' column)."""
    if stat_col not in deg.columns or log2fc_col not in deg.columns:
        return go.Figure().update_layout(title="Missing log2FC or stat column")

    df = deg.dropna(subset=[log2fc_col, stat_col]).copy()
    df["significant"] = "Not significant"
    if padj_col in df.columns:
        df.loc[df[padj_col] <= padj_threshold, "significant"] = "Significant"

    fig = px.scatter(
        df,
        x=log2fc_col,
        y=stat_col,
        color="significant",
        color_discrete_map={
            "Significant": "#D9534F",
            "Not significant": "#B0BEC5",
        },
        hover_data=["gene"] if "gene" in df.columns else None,
        labels={
            log2fc_col: "log2 Fold Change",
            stat_col: "Wald Statistic",
        },
        opacity=0.6,
    )
    fig.update_layout(
        title="log2FC vs Wald Statistic",
        template="plotly_white",
        height=500,
    )
    return fig


# ---------------------------------------------------------------------------
# LFC standard error plot
# ---------------------------------------------------------------------------

def lfc_se_scatter(
    deg: pd.DataFrame,
    basemean_col: str = "baseMean",
    lfcse_col: str = "lfcSE",
) -> go.Figure:
    """Scatter of log10(baseMean) vs lfcSE — shows shrinkage effects."""
    if basemean_col not in deg.columns or lfcse_col not in deg.columns:
        return go.Figure().update_layout(
            title="Missing baseMean or lfcSE columns"
        )

    df = deg.dropna(subset=[basemean_col, lfcse_col]).copy()
    df["log10_basemean"] = np.log10(df[basemean_col].clip(lower=0.1))

    fig = px.scatter(
        df,
        x="log10_basemean",
        y=lfcse_col,
        hover_data=["gene"] if "gene" in df.columns else None,
        labels={
            "log10_basemean": "log10(Mean Expression)",
            lfcse_col: "log2FC Standard Error",
        },
        opacity=0.5,
        color_discrete_sequence=["#4A90D9"],
    )
    fig.update_layout(
        title="Expression Level vs Fold Change Uncertainty",
        template="plotly_white",
        height=500,
    )
    return fig


# ---------------------------------------------------------------------------
# Summary counts
# ---------------------------------------------------------------------------

def deg_summary_counts(
    deg: pd.DataFrame,
    log2fc_col: str = "log2fc",
    padj_col: str = "padj",
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> dict:
    """Compute summary statistics for the DEG table.

    "Significant" means passing BOTH padj AND |log2FC| thresholds,
    which matches the actual genes used for GeneWalk and ORA.
    """
    total = len(deg)
    has_padj = padj_col in deg.columns
    has_lfc = log2fc_col in deg.columns

    if has_padj and has_lfc:
        sig_padj = deg[deg[padj_col] <= padj_threshold]
        up = sig_padj[sig_padj[log2fc_col] >= fc_threshold]
        down = sig_padj[sig_padj[log2fc_col] <= -fc_threshold]
        deg_count = len(up) + len(down)
    elif has_padj:
        up = pd.DataFrame()
        down = pd.DataFrame()
        deg_count = len(deg[deg[padj_col] <= padj_threshold])
    else:
        up = pd.DataFrame()
        down = pd.DataFrame()
        deg_count = 0

    return {
        "total_genes": total,
        "significant": deg_count,
        "up_regulated": len(up),
        "down_regulated": len(down),
        "not_significant": total - deg_count,
        "pct_significant": round(100 * deg_count / total, 1) if total > 0 else 0,
    }


def deg_category_pie(
    deg: pd.DataFrame,
    log2fc_col: str = "log2fc",
    padj_col: str = "padj",
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> go.Figure:
    """Pie chart showing proportions of up/down/not-significant genes."""
    counts = deg_summary_counts(deg, log2fc_col, padj_col, fc_threshold, padj_threshold)

    labels = ["Up-regulated", "Down-regulated", "Not significant"]
    values = [counts["up_regulated"], counts["down_regulated"], counts["not_significant"]]
    colors = ["#D9534F", "#4A90D9", "#B0BEC5"]

    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        marker_colors=colors,
        textinfo="label+value+percent",
        hole=0.35,
    )])
    fig.update_layout(
        title="Gene Classification Summary",
        template="plotly_white",
        height=400,
    )
    return fig


def basemean_distribution(
    deg: pd.DataFrame,
    basemean_col: str = "baseMean",
) -> go.Figure:
    """Histogram of log10(baseMean) to show expression level distribution."""
    if basemean_col not in deg.columns:
        return go.Figure().update_layout(title=f"Column '{basemean_col}' not found")

    df = deg.dropna(subset=[basemean_col]).copy()
    df["log10_basemean"] = np.log10(df[basemean_col].clip(lower=0.1))

    fig = px.histogram(
        df,
        x="log10_basemean",
        nbins=50,
        color_discrete_sequence=["#F0AD4E"],
        labels={"log10_basemean": "log10(Mean Expression)"},
    )
    fig.update_layout(
        title="Expression Level Distribution",
        template="plotly_white",
        height=400,
        yaxis_title="Count",
    )
    return fig
