"""Visualization functions for GSEA results and cross-method comparisons."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from genewalk_app.visualizations import GO_DOMAIN_COLORS


# ---------------------------------------------------------------------------
# GSEA-specific plots
# ---------------------------------------------------------------------------

def nes_bar_chart(
    gsea_df: pd.DataFrame,
    fdr_threshold: float = 0.25,
    top_n: int = 30,
) -> go.Figure:
    """Bar chart of Normalized Enrichment Scores, colored by direction."""
    if gsea_df.empty or "nes" not in gsea_df.columns or "term" not in gsea_df.columns:
        return go.Figure().update_layout(title="No GSEA results to display")

    df = gsea_df.copy()
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold]
    if df.empty:
        return go.Figure().update_layout(
            title=f"No significant pathways (FDR <= {fdr_threshold})"
        )

    df = df.sort_values("nes", key=abs, ascending=False).head(top_n)
    df = df.sort_values("nes", ascending=True)
    df["direction"] = df["nes"].apply(lambda x: "Up-regulated" if x > 0 else "Down-regulated")

    fig = px.bar(
        df,
        x="nes",
        y="term",
        color="direction",
        color_discrete_map={
            "Up-regulated": "#D9534F",
            "Down-regulated": "#4A90D9",
        },
        orientation="h",
        hover_data=[c for c in ["fdr", "pval", "gene_set_library"] if c in df.columns],
        labels={"nes": "Normalized Enrichment Score", "term": ""},
    )
    fig.update_layout(
        title="GSEA: Top Enriched Pathways by NES",
        template="plotly_white",
        height=max(450, len(df) * 25),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


def gsea_dot_plot(
    gsea_df: pd.DataFrame,
    fdr_threshold: float = 0.25,
    top_n: int = 40,
) -> go.Figure:
    """Dot plot: NES vs term, sized by gene count, colored by FDR."""
    if gsea_df.empty or "nes" not in gsea_df.columns:
        return go.Figure().update_layout(title="No GSEA results to display")

    df = gsea_df.copy()
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold]
    if df.empty:
        return go.Figure().update_layout(
            title=f"No significant pathways (FDR <= {fdr_threshold})"
        )

    df = df.sort_values("nes", key=abs, ascending=False).head(top_n)

    # Estimate gene set size from tag_pct if available
    if "tag_pct" in df.columns:
        df["tag_pct_clean"] = pd.to_numeric(
            df["tag_pct"].astype(str).str.rstrip("%"), errors="coerce"
        ).fillna(10)
    else:
        df["tag_pct_clean"] = 10

    if "fdr" in df.columns:
        df["neg_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))
        color_col = "neg_log10_fdr"
        color_label = "-log10(FDR)"
    else:
        color_col = "nes"
        color_label = "NES"

    fig = px.scatter(
        df,
        x="nes",
        y="term",
        size="tag_pct_clean",
        color=color_col,
        color_continuous_scale="RdBu_r",
        hover_data=[c for c in ["fdr", "pval", "gene_set_library"] if c in df.columns],
        labels={
            "nes": "Normalized Enrichment Score",
            "term": "",
            color_col: color_label,
            "tag_pct_clean": "Gene %",
        },
    )
    fig.update_layout(
        title="GSEA Dot Plot",
        template="plotly_white",
        height=max(450, len(df) * 22),
    )
    return fig


# ---------------------------------------------------------------------------
# GSEA leading edge
# ---------------------------------------------------------------------------

def gsea_leading_edge_table(
    gsea_df: pd.DataFrame,
    fdr_threshold: float = 0.25,
    top_n: int = 20,
) -> pd.DataFrame:
    """Extract a summary table of leading edge genes from GSEA results.

    Returns a DataFrame with columns: term, nes, fdr, lead_genes,
    n_lead_genes, gene_set_library.  This is the most biologically
    actionable output from GSEA — the core genes driving each enrichment.
    """
    if gsea_df.empty or "nes" not in gsea_df.columns:
        return pd.DataFrame()

    df = gsea_df.copy()
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold]
    if df.empty:
        return pd.DataFrame()

    df = df.sort_values("nes", key=abs, ascending=False).head(top_n)

    # Identify the lead genes column (varies by GSEApy version)
    lead_col = None
    for candidate in ["lead_genes", "Lead_genes", "genes", "ledge_genes"]:
        if candidate in df.columns:
            lead_col = candidate
            break

    cols = ["term", "nes"]
    if "fdr" in df.columns:
        cols.append("fdr")
    if lead_col:
        df["lead_genes"] = df[lead_col].astype(str)
        df["n_lead_genes"] = df["lead_genes"].apply(
            lambda x: len(x.split(";")) if x and x != "nan" else 0
        )
        cols.extend(["lead_genes", "n_lead_genes"])
    if "gene_set_library" in df.columns:
        cols.append("gene_set_library")

    return df[[c for c in cols if c in df.columns]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# ORA plots
# ---------------------------------------------------------------------------

def ora_dot_plot(
    ora_df: pd.DataFrame,
    label: str = "",
    fdr_threshold: float = 0.05,
    top_n: int = 25,
) -> go.Figure:
    """Dot plot for ORA results: sized by overlap ratio, colored by FDR.

    Uses the parsed overlap_ratio and overlap_count columns from
    _clean_ora_results for proper numeric sizing.
    """
    if ora_df.empty or "term" not in ora_df.columns:
        return go.Figure().update_layout(
            title=f"No ORA results{' (' + label + ')' if label else ''}"
        )

    df = ora_df.copy()
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold]
    if df.empty:
        return go.Figure().update_layout(
            title=f"No significant terms (FDR <= {fdr_threshold})"
        )

    if "combined_score" in df.columns:
        df = df.sort_values("combined_score", ascending=False).head(top_n)
    elif "fdr" in df.columns:
        df = df.sort_values("fdr", ascending=True).head(top_n)

    # Size by overlap count if available
    size_col = "overlap_count" if "overlap_count" in df.columns else None
    if size_col is None:
        df["_dot_size"] = 10
        size_col = "_dot_size"

    if "fdr" in df.columns:
        df["neg_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))
        color_col = "neg_log10_fdr"
        color_label = "-log10(FDR)"
    elif "combined_score" in df.columns:
        color_col = "combined_score"
        color_label = "Combined Score"
    else:
        df["_color"] = 1
        color_col = "_color"
        color_label = ""

    hover_cols = [c for c in ["fdr", "overlap", "overlap_ratio", "genes"] if c in df.columns]

    fig = px.scatter(
        df,
        x="combined_score" if "combined_score" in df.columns else "neg_log10_fdr",
        y="term",
        size=size_col,
        color=color_col,
        color_continuous_scale="Reds",
        hover_data=hover_cols,
        labels={
            "combined_score": "Combined Score",
            "neg_log10_fdr": "-log10(FDR)",
            "term": "",
            color_col: color_label,
            "overlap_count": "Genes in Overlap",
        },
    )
    title = f"ORA Dot Plot{' (' + label + ')' if label else ''}"
    fig.update_layout(
        title=title,
        template="plotly_white",
        height=max(450, len(df) * 22),
    )
    return fig


def ora_bar_chart(
    ora_df: pd.DataFrame,
    label: str = "",
    fdr_threshold: float = 0.05,
    top_n: int = 20,
) -> go.Figure:
    """Horizontal bar chart of ORA results."""
    if ora_df.empty or "term" not in ora_df.columns:
        return go.Figure().update_layout(title=f"No ORA results{' (' + label + ')' if label else ''}")

    df = ora_df.copy()
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold]
    if df.empty:
        return go.Figure().update_layout(title=f"No significant terms (FDR <= {fdr_threshold})")

    if "combined_score" in df.columns:
        df = df.sort_values("combined_score", ascending=False).head(top_n)
        x_col = "combined_score"
        x_label = "Combined Score"
    elif "fdr" in df.columns:
        df = df.sort_values("fdr", ascending=True).head(top_n)
        df["neg_log10_fdr"] = -np.log10(df["fdr"].clip(lower=1e-300))
        x_col = "neg_log10_fdr"
        x_label = "-log10(FDR)"
    else:
        return go.Figure().update_layout(title="Missing score columns")

    fig = px.bar(
        df,
        x=x_col,
        y="term",
        orientation="h",
        color="gene_set_library" if "gene_set_library" in df.columns else None,
        hover_data=[c for c in ["fdr", "pval", "overlap", "genes"] if c in df.columns],
        labels={x_col: x_label, "term": ""},
    )
    title = f"ORA: Top Enriched Terms{' (' + label + ')' if label else ''}"
    fig.update_layout(
        title=title,
        template="plotly_white",
        height=max(400, len(df) * 25),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


# ---------------------------------------------------------------------------
# GeneWalk comparison plots
# ---------------------------------------------------------------------------

def direction_volcano(
    gw_up: pd.DataFrame,
    gw_down: pd.DataFrame,
    padj_col: str = "global_padj",
    padj_threshold: float = 0.05,
) -> go.Figure:
    """Combined volcano plot with gene-GO pairs from both directions."""
    required = {"hgnc_symbol", "go_name", "sim"}
    for label, gw in [("up", gw_up), ("down", gw_down)]:
        missing = required - set(gw.columns)
        if missing:
            return go.Figure().update_layout(
                title=f"Missing columns in {label}-regulated GeneWalk results: {', '.join(missing)}"
            )
        if padj_col not in gw.columns:
            return go.Figure().update_layout(title=f"Missing {padj_col} column in {label}-regulated results")

    keep_cols = ["hgnc_symbol", "go_name", "sim", padj_col]
    up = gw_up[keep_cols].copy()
    up["direction"] = "Up-regulated"
    down = gw_down[keep_cols].copy()
    down["direction"] = "Down-regulated"

    combined = pd.concat([up, down], ignore_index=True).dropna(subset=[padj_col, "sim"])
    combined["neg_log10_padj"] = -np.log10(combined[padj_col].clip(lower=1e-300))

    fig = px.scatter(
        combined,
        x="sim",
        y="neg_log10_padj",
        color="direction",
        color_discrete_map={
            "Up-regulated": "#D9534F",
            "Down-regulated": "#4A90D9",
        },
        hover_data=["hgnc_symbol", "go_name"],
        labels={
            "sim": "Similarity Score",
            "neg_log10_padj": f"-log10({padj_col})",
        },
        opacity=0.6,
    )
    fig.add_hline(
        y=-np.log10(padj_threshold), line_dash="dash", line_color="#888",
        annotation_text=f"padj = {padj_threshold}",
    )
    fig.update_layout(
        title="GeneWalk: Up vs Down Regulated Gene-GO Associations",
        template="plotly_white",
        height=550,
    )
    return fig


def shared_terms_bar(
    shared_df: pd.DataFrame,
    top_n: int = 25,
) -> go.Figure:
    """Bar chart of GO terms found significant in both up and down analyses."""
    if shared_df.empty:
        return go.Figure().update_layout(title="No shared significant GO terms found")

    padj_cols = [c for c in shared_df.columns if c.startswith("best_padj")]
    if len(padj_cols) < 2:
        return go.Figure().update_layout(title="Insufficient data for shared terms plot")

    df = shared_df.head(top_n).copy()
    up_col = padj_cols[0]
    down_col = padj_cols[1]
    df["neg_log10_up"] = -np.log10(df[up_col].clip(lower=1e-300))
    df["neg_log10_down"] = -np.log10(df[down_col].clip(lower=1e-300))

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=df["go_name"], x=df["neg_log10_up"],
        name="Up-regulated", orientation="h",
        marker_color="#D9534F", opacity=0.8,
    ))
    fig.add_trace(go.Bar(
        y=df["go_name"], x=df["neg_log10_down"],
        name="Down-regulated", orientation="h",
        marker_color="#4A90D9", opacity=0.8,
    ))
    fig.update_layout(
        title="GO Terms Significant in Both Directions",
        barmode="group",
        template="plotly_white",
        height=max(400, len(df) * 30),
        xaxis_title="-log10(adjusted p-value)",
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


def concordance_bar(
    concordance_df: pd.DataFrame,
    top_n: int = 30,
) -> go.Figure:
    """Bar chart showing terms found by multiple analysis methods."""
    if concordance_df.empty:
        return go.Figure().update_layout(title="No concordant terms across methods")

    df = concordance_df.head(top_n).copy()

    fig = px.bar(
        df,
        x="n_methods",
        y="term",
        color="found_in",
        orientation="h",
        hover_data=[c for c in ["gsea_fdr", "ora_fdr"] if c in df.columns],
        labels={"n_methods": "Number of Methods", "term": "", "found_in": "Found in"},
    )
    fig.update_layout(
        title="Terms Found by Multiple Analysis Methods",
        template="plotly_white",
        height=max(400, len(df) * 25),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


def deg_overview_volcano(
    deg: pd.DataFrame,
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> go.Figure:
    """Standard DEG volcano plot: log2FC vs -log10(padj)."""
    if "log2fc" not in deg.columns or "padj" not in deg.columns:
        return go.Figure().update_layout(title="Missing log2fc/padj columns")

    df = deg.dropna(subset=["log2fc", "padj"]).copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["category"] = "Not significant"
    df.loc[
        (df["padj"] <= padj_threshold) & (df["log2fc"] >= fc_threshold), "category"
    ] = "Up-regulated"
    df.loc[
        (df["padj"] <= padj_threshold) & (df["log2fc"] <= -fc_threshold), "category"
    ] = "Down-regulated"

    fig = px.scatter(
        df,
        x="log2fc",
        y="neg_log10_padj",
        color="category",
        color_discrete_map={
            "Up-regulated": "#D9534F",
            "Down-regulated": "#4A90D9",
            "Not significant": "#B0BEC5",
        },
        hover_data=["gene"] if "gene" in df.columns else None,
        labels={
            "log2fc": "log2 Fold Change",
            "neg_log10_padj": "-log10(adjusted p-value)",
        },
        opacity=0.6,
    )
    fig.add_vline(x=fc_threshold, line_dash="dash", line_color="#888")
    fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="#888")
    fig.add_hline(
        y=-np.log10(padj_threshold), line_dash="dash", line_color="#888",
        annotation_text=f"padj = {padj_threshold}",
    )
    fig.update_layout(
        title="Differential Expression Volcano Plot",
        template="plotly_white",
        height=550,
    )
    return fig


def gene_set_summary_metrics(
    up_genes: list[str],
    down_genes: list[str],
    gw_up: pd.DataFrame | None,
    gw_down: pd.DataFrame | None,
    gsea_df: pd.DataFrame | None,
    ora_up: pd.DataFrame | None,
    ora_down: pd.DataFrame | None,
    padj_col: str = "global_padj",
    padj_threshold: float = 0.05,
    gsea_fdr: float = 0.25,
    ora_fdr: float = 0.05,
) -> dict:
    """Compute summary metrics for the overview tab."""
    metrics: dict = {
        "n_up": len(up_genes),
        "n_down": len(down_genes),
    }

    for label, gw in [("up", gw_up), ("down", gw_down)]:
        if gw is not None and not gw.empty and padj_col in gw.columns:
            sig = gw[gw[padj_col] <= padj_threshold]
            metrics[f"gw_{label}_sig_pairs"] = len(sig)
            metrics[f"gw_{label}_sig_genes"] = (
                sig["hgnc_symbol"].nunique() if "hgnc_symbol" in sig.columns else 0
            )
            metrics[f"gw_{label}_sig_terms"] = (
                sig["go_name"].nunique() if "go_name" in sig.columns else 0
            )
        else:
            metrics[f"gw_{label}_sig_pairs"] = 0
            metrics[f"gw_{label}_sig_genes"] = 0
            metrics[f"gw_{label}_sig_terms"] = 0

    if gsea_df is not None and not gsea_df.empty and "fdr" in gsea_df.columns:
        sig_gsea = gsea_df[gsea_df["fdr"] <= gsea_fdr]
        if "nes" in sig_gsea.columns:
            metrics["gsea_sig_up"] = int((sig_gsea["nes"] > 0).sum())
            metrics["gsea_sig_down"] = int((sig_gsea["nes"] < 0).sum())
        else:
            metrics["gsea_sig_up"] = len(sig_gsea)
            metrics["gsea_sig_down"] = 0
    else:
        metrics["gsea_sig_up"] = 0
        metrics["gsea_sig_down"] = 0

    for label, ora in [("up", ora_up), ("down", ora_down)]:
        if ora is not None and not ora.empty and "fdr" in ora.columns:
            metrics[f"ora_{label}_sig"] = len(ora[ora["fdr"] <= ora_fdr])
        else:
            metrics[f"ora_{label}_sig"] = 0

    return metrics
