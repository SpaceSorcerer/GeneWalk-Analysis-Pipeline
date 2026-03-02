"""Splicing analysis visualizations.

All functions return Plotly figures that can be rendered with
``st.plotly_chart(fig, width="stretch")``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# ───────────────────────────────────────────────────────────────────────
# Delta-PSI volcano
# ───────────────────────────────────────────────────────────────────────

def dpsi_volcano(
    df: pd.DataFrame,
    dpsi_threshold: float = 0.1,
    fdr_threshold: float = 0.05,
) -> go.Figure:
    """Volcano plot: delta-PSI vs -log10(p-value/FDR).

    Points colored by significance and direction of splicing change.
    """
    if df.empty or "dpsi" not in df.columns:
        return go.Figure().update_layout(title="No splicing data to display")

    plot = df.dropna(subset=["dpsi"]).copy()

    stat_col = "fdr" if "fdr" in plot.columns and plot["fdr"].notna().any() else "pvalue"
    if stat_col not in plot.columns or plot[stat_col].isna().all():
        return go.Figure().update_layout(
            title="No p-value/FDR available for volcano plot"
        )

    plot["neg_log10_p"] = -np.log10(plot[stat_col].clip(lower=1e-300))

    plot["category"] = "Not significant"
    sig_mask = plot[stat_col] <= fdr_threshold
    plot.loc[sig_mask & (plot["dpsi"] >= dpsi_threshold), "category"] = "Inclusion increased"
    plot.loc[sig_mask & (plot["dpsi"] <= -dpsi_threshold), "category"] = "Inclusion decreased"

    hover_cols = [c for c in ["gene", "event_id", "event_type", "dpsi", stat_col]
                  if c in plot.columns]

    fig = px.scatter(
        plot,
        x="dpsi",
        y="neg_log10_p",
        color="category",
        color_discrete_map={
            "Inclusion increased": "#D9534F",
            "Inclusion decreased": "#4A90D9",
            "Not significant": "#B0BEC5",
        },
        hover_data=hover_cols,
        labels={
            "dpsi": "\u0394PSI (Condition 2 \u2212 Condition 1)",
            "neg_log10_p": f"-log10({stat_col.upper()})",
        },
        opacity=0.6,
    )
    fig.add_vline(x=dpsi_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.add_vline(x=-dpsi_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.add_hline(
        y=-np.log10(fdr_threshold), line_dash="dash", line_color="#888", line_width=0.8,
    )
    fig.update_layout(
        title="Splicing Volcano Plot (\u0394PSI vs Significance)",
        template="plotly_white",
        height=550,
    )
    return fig


# ───────────────────────────────────────────────────────────────────────
# Event type distribution
# ───────────────────────────────────────────────────────────────────────

def event_type_summary(
    df: pd.DataFrame,
    show_all: bool = True,
) -> go.Figure:
    """Bar chart of splicing event counts by type (SE, RI, A3SS, etc.)."""
    if df.empty or "event_type" not in df.columns:
        return go.Figure().update_layout(title="No event type data")

    counts = df["event_type"].value_counts().reset_index()
    counts.columns = ["Event Type", "Count"]

    type_colors = {
        "SE": "#4A90D9",
        "RI": "#D9534F",
        "A3SS": "#5CB85C",
        "A5SS": "#F0AD4E",
        "MXE": "#9B59B6",
        "MIC": "#E67E22",
        "other": "#95A5A6",
    }

    fig = px.bar(
        counts,
        x="Event Type",
        y="Count",
        color="Event Type",
        color_discrete_map=type_colors,
        text="Count",
    )
    fig.update_layout(
        title="Splicing Events by Type",
        template="plotly_white",
        height=400,
        showlegend=False,
    )
    fig.update_traces(textposition="outside")
    return fig


def event_type_pie(df: pd.DataFrame) -> go.Figure:
    """Pie chart of splicing event type distribution."""
    if df.empty or "event_type" not in df.columns:
        return go.Figure().update_layout(title="No event type data")

    counts = df["event_type"].value_counts()

    type_colors = {
        "SE": "#4A90D9",
        "RI": "#D9534F",
        "A3SS": "#5CB85C",
        "A5SS": "#F0AD4E",
        "MXE": "#9B59B6",
        "MIC": "#E67E22",
        "other": "#95A5A6",
    }
    colors = [type_colors.get(et, "#95A5A6") for et in counts.index]

    fig = go.Figure(data=[go.Pie(
        labels=counts.index.tolist(),
        values=counts.values.tolist(),
        marker_colors=colors,
        textinfo="label+value+percent",
        hole=0.35,
    )])
    fig.update_layout(
        title="Splicing Event Types",
        template="plotly_white",
        height=400,
    )
    return fig


# ───────────────────────────────────────────────────────────────────────
# PSI distribution
# ───────────────────────────────────────────────────────────────────────

def psi_distribution(
    df: pd.DataFrame,
) -> go.Figure:
    """Overlaid histograms of PSI values for condition 1 vs condition 2."""
    cols = [c for c in ("psi_1", "psi_2") if c in df.columns and df[c].notna().any()]
    if not cols:
        return go.Figure().update_layout(title="No PSI data available")

    fig = go.Figure()
    names = {"psi_1": "Condition 1", "psi_2": "Condition 2"}
    colors = {"psi_1": "#4A90D9", "psi_2": "#D9534F"}

    for col in cols:
        vals = df[col].dropna()
        fig.add_trace(go.Histogram(
            x=vals,
            name=names.get(col, col),
            marker_color=colors.get(col, "#888"),
            opacity=0.6,
            nbinsx=50,
        ))

    fig.update_layout(
        barmode="overlay",
        title="PSI Distribution by Condition",
        xaxis_title="PSI (Percent Spliced In)",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
    )
    return fig


def dpsi_distribution(df: pd.DataFrame) -> go.Figure:
    """Histogram of delta-PSI values."""
    if df.empty or "dpsi" not in df.columns:
        return go.Figure().update_layout(title="No \u0394PSI data")

    vals = df["dpsi"].dropna()

    fig = px.histogram(
        vals,
        nbins=60,
        color_discrete_sequence=["#4A90D9"],
        labels={"value": "\u0394PSI"},
    )
    fig.add_vline(x=0, line_color="#333", line_width=1)
    fig.update_layout(
        title="\u0394PSI Distribution",
        template="plotly_white",
        height=400,
        xaxis_title="\u0394PSI (Condition 2 \u2212 Condition 1)",
        yaxis_title="Count",
    )
    return fig


# ───────────────────────────────────────────────────────────────────────
# Top events bar chart
# ───────────────────────────────────────────────────────────────────────

def top_splicing_events_bar(
    df: pd.DataFrame,
    top_n: int = 30,
    sort_by: str = "dpsi",
) -> go.Figure:
    """Horizontal bar chart of top splicing events by |delta-PSI| or FDR."""
    if df.empty or "dpsi" not in df.columns:
        return go.Figure().update_layout(title="No splicing data")

    plot = df.dropna(subset=["dpsi"]).copy()

    if sort_by == "fdr" and "fdr" in plot.columns:
        plot = plot.sort_values("fdr", ascending=True).head(top_n)
    else:
        plot["abs_dpsi"] = plot["dpsi"].abs()
        plot = plot.sort_values("abs_dpsi", ascending=False).head(top_n)

    # Label: gene (event_type)
    if "event_type" in plot.columns:
        plot["label"] = plot["gene"] + " (" + plot["event_type"] + ")"
    else:
        plot["label"] = plot["gene"]

    plot = plot.sort_values("dpsi", ascending=True)

    plot["direction"] = plot["dpsi"].apply(
        lambda x: "Inclusion increased" if x > 0 else "Inclusion decreased"
    )

    fig = px.bar(
        plot,
        x="dpsi",
        y="label",
        color="direction",
        color_discrete_map={
            "Inclusion increased": "#D9534F",
            "Inclusion decreased": "#4A90D9",
        },
        orientation="h",
        hover_data=[c for c in ["fdr", "psi_1", "psi_2", "event_id"] if c in plot.columns],
        labels={"dpsi": "\u0394PSI", "label": ""},
    )
    fig.update_layout(
        title=f"Top {min(top_n, len(plot))} Differential Splicing Events",
        template="plotly_white",
        height=max(450, len(plot) * 22),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


# ───────────────────────────────────────────────────────────────────────
# Per-gene splicing detail
# ───────────────────────────────────────────────────────────────────────

def gene_splicing_detail(
    df: pd.DataFrame,
    gene: str,
) -> go.Figure:
    """Dot plot of all splicing events for a single gene.

    Shows PSI in condition 1 and condition 2 as paired dots connected
    by lines, grouped by event type.
    """
    gdf = df[df["gene"].str.upper() == gene.upper()].copy()
    if gdf.empty:
        return go.Figure().update_layout(title=f"No splicing events for {gene}")

    # Create event labels
    if "event_type" in gdf.columns:
        gdf["label"] = gdf["event_type"] + ": " + gdf["event_id"]
    else:
        gdf["label"] = gdf["event_id"]

    fig = go.Figure()

    has_psi1 = "psi_1" in gdf.columns and gdf["psi_1"].notna().any()
    has_psi2 = "psi_2" in gdf.columns and gdf["psi_2"].notna().any()

    if has_psi1 and has_psi2:
        # Paired dot plot
        for _, row in gdf.iterrows():
            fig.add_trace(go.Scatter(
                x=[row.get("psi_1"), row.get("psi_2")],
                y=[row["label"], row["label"]],
                mode="lines+markers",
                marker=dict(size=10),
                line=dict(color="#888", width=1),
                showlegend=False,
                hoverinfo="text",
                text=[
                    f"Cond 1: PSI={row.get('psi_1', float('nan')):.2f}"
                    if pd.notna(row.get("psi_1")) else "Cond 1: PSI=N/A",
                    f"Cond 2: PSI={row.get('psi_2', float('nan')):.2f}"
                    if pd.notna(row.get("psi_2")) else "Cond 2: PSI=N/A",
                ],
            ))

        # Add condition markers on top
        fig.add_trace(go.Scatter(
            x=gdf["psi_1"],
            y=gdf["label"],
            mode="markers",
            name="Condition 1",
            marker=dict(size=10, color="#4A90D9"),
        ))
        fig.add_trace(go.Scatter(
            x=gdf["psi_2"],
            y=gdf["label"],
            mode="markers",
            name="Condition 2",
            marker=dict(size=10, color="#D9534F"),
        ))

        fig.update_layout(
            title=f"Splicing Events for {gene}",
            xaxis_title="PSI (Percent Spliced In)",
            template="plotly_white",
            height=max(350, len(gdf) * 40),
        )
    elif "dpsi" in gdf.columns:
        # Fallback: just show delta-PSI as bars
        fig = px.bar(
            gdf,
            x="dpsi",
            y="label",
            orientation="h",
            color="dpsi",
            color_continuous_scale="RdBu_r",
            color_continuous_midpoint=0,
            labels={"dpsi": "\u0394PSI", "label": ""},
        )
        fig.update_layout(
            title=f"Splicing Events for {gene}",
            template="plotly_white",
            height=max(350, len(gdf) * 40),
        )

    return fig


# ───────────────────────────────────────────────────────────────────────
# Genes with most splicing events
# ───────────────────────────────────────────────────────────────────────

def genes_by_event_count(
    df: pd.DataFrame,
    top_n: int = 30,
) -> go.Figure:
    """Bar chart: genes ranked by number of differential splicing events."""
    if df.empty or "gene" not in df.columns:
        return go.Figure().update_layout(title="No splicing data")

    counts = df["gene"].value_counts().head(top_n).reset_index()
    counts.columns = ["gene", "n_events"]
    counts = counts.sort_values("n_events", ascending=True)

    fig = px.bar(
        counts,
        x="n_events",
        y="gene",
        orientation="h",
        color_discrete_sequence=["#4A90D9"],
        labels={"n_events": "Differential Splicing Events", "gene": ""},
    )
    fig.update_layout(
        title=f"Top {min(top_n, len(counts))} Genes by Splicing Event Count",
        template="plotly_white",
        height=max(450, len(counts) * 22),
    )
    return fig
