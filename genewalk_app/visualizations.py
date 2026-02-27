"""Interactive Plotly visualizations for GeneWalk results."""

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import networkx as nx

# Unified GO domain color palette used across all charts.
GO_DOMAIN_COLORS = {
    "biological_process": "#4A90D9",
    "molecular_function": "#F0AD4E",
    "cellular_component": "#5CB85C",
}


def volcano_plot(df: pd.DataFrame, padj_col: str = "gene_padj", padj_threshold: float = 0.05) -> go.Figure:
    """Volcano-style plot: similarity score vs -log10(adjusted p-value)."""
    if padj_col not in df.columns or "sim" not in df.columns:
        return go.Figure().update_layout(title="Missing required columns for volcano plot")

    plot_df = df.dropna(subset=[padj_col, "sim"]).copy()
    plot_df["neg_log10_padj"] = -np.log10(plot_df[padj_col].clip(lower=1e-300))
    plot_df["significant"] = plot_df[padj_col].le(padj_threshold).map(
        {True: "Significant", False: "Not significant"}
    )

    fig = px.scatter(
        plot_df,
        x="sim",
        y="neg_log10_padj",
        color="significant",
        color_discrete_map={"Significant": "#D9534F", "Not significant": "#B0BEC5"},
        hover_data=["hgnc_symbol", "go_name", padj_col, "sim"],
        labels={
            "sim": "Similarity Score",
            "neg_log10_padj": "-log10(Adjusted P-value)",
            "significant": "",
        },
    )

    fig.add_hline(
        y=-np.log10(padj_threshold),
        line_dash="dash",
        line_color="#888",
        annotation_text=f"padj = {padj_threshold}",
    )

    fig.update_layout(
        title="Volcano Plot: Gene-GO Term Associations",
        template="plotly_white",
        height=550,
    )
    return fig


def gene_bar_chart(
    df: pd.DataFrame,
    gene: str,
    padj_col: str = "gene_padj",
    top_n: int = 20,
) -> go.Figure:
    """Horizontal bar chart of top GO terms for a single gene."""
    gene_df = df[df["hgnc_symbol"] == gene].copy()
    if gene_df.empty:
        return go.Figure().update_layout(title=f"No results for {gene}")

    if padj_col in gene_df.columns:
        gene_df = gene_df.sort_values(padj_col, ascending=True).head(top_n)
        gene_df["neg_log10_padj"] = -np.log10(gene_df[padj_col].clip(lower=1e-300))
    else:
        gene_df = gene_df.head(top_n)

    # Color by GO domain if available
    color_col = "go_domain" if "go_domain" in gene_df.columns else None

    fig = px.bar(
        gene_df,
        x="neg_log10_padj" if padj_col in df.columns else "sim",
        y="go_name",
        color=color_col,
        orientation="h",
        hover_data=[c for c in ["go_id", padj_col, "sim"] if c in gene_df.columns],
        color_discrete_map=GO_DOMAIN_COLORS,
        labels={
            "neg_log10_padj": "-log10(Adjusted P-value)",
            "go_name": "GO Term",
            "sim": "Similarity",
        },
    )

    fig.update_layout(
        title=f"Top GO Terms for {gene}",
        template="plotly_white",
        height=max(400, top_n * 25),
        yaxis={"categoryorder": "total ascending"},
    )
    return fig


def go_domain_pie(df: pd.DataFrame, padj_col: str = "gene_padj", padj_threshold: float = 0.05) -> go.Figure:
    """Pie chart of GO domain distribution among significant results."""
    if "go_domain" not in df.columns:
        return go.Figure().update_layout(title="GO domain column not found")

    sig = df[df[padj_col] <= padj_threshold] if padj_col in df.columns else df
    counts = sig["go_domain"].value_counts().reset_index()
    counts.columns = ["go_domain", "count"]

    fig = px.pie(
        counts,
        values="count",
        names="go_domain",
        color="go_domain",
        color_discrete_map=GO_DOMAIN_COLORS,
    )
    fig.update_layout(
        title="GO Domain Distribution (Significant Terms)",
        template="plotly_white",
    )
    return fig


def gene_similarity_heatmap(
    df: pd.DataFrame,
    genes: list[str] | None = None,
    padj_col: str = "gene_padj",
    padj_threshold: float = 0.05,
    max_terms: int = 30,
) -> go.Figure:
    """Heatmap of similarity scores: genes vs GO terms."""
    sig = df[df[padj_col] <= padj_threshold] if padj_col in df.columns else df

    if genes:
        sig = sig[sig["hgnc_symbol"].isin(genes)]

    if sig.empty:
        return go.Figure().update_layout(title="No significant results to display")

    # Pick top GO terms by frequency
    top_terms = sig["go_name"].value_counts().head(max_terms).index.tolist()
    sig = sig[sig["go_name"].isin(top_terms)]

    pivot = sig.pivot_table(
        index="hgnc_symbol",
        columns="go_name",
        values="sim",
        aggfunc="mean",
    ).fillna(0)

    fig = px.imshow(
        pivot,
        color_continuous_scale="Blues",
        labels={"color": "Similarity"},
        aspect="auto",
    )
    fig.update_layout(
        title="Gene-GO Term Similarity Heatmap",
        template="plotly_white",
        height=max(400, len(pivot) * 28),
    )
    return fig


def summary_bar(summary_df: pd.DataFrame) -> go.Figure:
    """Bar chart of significant GO term counts per gene."""
    if summary_df.empty:
        return go.Figure().update_layout(title="No summary data")

    fig = px.bar(
        summary_df.head(30),
        x="hgnc_symbol",
        y="significant_go_terms",
        color="mean_similarity",
        color_continuous_scale="Blues",
        hover_data=["top_go_term", "best_padj"],
        labels={
            "hgnc_symbol": "Gene",
            "significant_go_terms": "Significant GO Terms",
            "mean_similarity": "Mean Similarity",
        },
    )
    fig.update_layout(
        title="Genes Ranked by Number of Significant GO Associations",
        template="plotly_white",
        height=450,
        xaxis_tickangle=-45,
    )
    return fig


def pvalue_distribution(df: pd.DataFrame, padj_col: str = "gene_padj") -> go.Figure:
    """Histogram of adjusted p-value distribution."""
    if padj_col not in df.columns:
        return go.Figure().update_layout(title="P-value column not found")

    fig = px.histogram(
        df,
        x=padj_col,
        nbins=50,
        labels={padj_col: "Adjusted P-value"},
        color_discrete_sequence=["#4A90D9"],
    )
    fig.update_layout(
        title="Distribution of Adjusted P-values",
        template="plotly_white",
        height=400,
    )
    return fig


def gene_go_network(
    df: pd.DataFrame,
    genes: list[str] | None = None,
    padj_col: str = "gene_padj",
    padj_threshold: float = 0.05,
    max_edges: int = 200,
) -> go.Figure:
    """Interactive network graph: genes connected to their significant GO terms."""
    sig = df[df[padj_col] <= padj_threshold] if padj_col in df.columns else df

    if genes:
        sig = sig[sig["hgnc_symbol"].isin(genes)]

    if sig.empty:
        return go.Figure().update_layout(title="No significant results to display")

    # Limit edges for readability
    if len(sig) > max_edges:
        sig = sig.sort_values(padj_col, ascending=True).head(max_edges)

    G = nx.Graph()

    for _, row in sig.iterrows():
        gene = row.get("hgnc_symbol", "")
        go_term = row.get("go_name", "")
        sim = row.get("sim", 0)
        padj = row.get(padj_col, 1)
        domain = row.get("go_domain", "unknown")
        if not gene or not go_term:
            continue
        G.add_node(gene, node_type="gene")
        G.add_node(go_term, node_type="go", domain=domain)
        G.add_edge(gene, go_term, sim=sim, padj=padj)

    if len(G.nodes) == 0:
        return go.Figure().update_layout(title="No nodes to display")

    pos = nx.spring_layout(G, k=1.8 / (len(G.nodes) ** 0.3), seed=42, iterations=60)

    # Edges
    edge_x, edge_y = [], []
    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        mode="lines",
        line={"width": 0.5, "color": "#C0C0C0"},
        hoverinfo="none",
    )

    # Gene nodes
    gene_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "gene"]
    gene_x = [pos[n][0] for n in gene_nodes]
    gene_y = [pos[n][1] for n in gene_nodes]
    gene_sizes = [8 + 2 * G.degree(n) for n in gene_nodes]
    gene_hover = [f"<b>{n}</b><br>Connections: {G.degree(n)}" for n in gene_nodes]

    gene_trace = go.Scatter(
        x=gene_x, y=gene_y,
        mode="markers+text",
        marker={"size": gene_sizes, "color": "#D9534F", "line": {"width": 1, "color": "#fff"}},
        text=gene_nodes,
        textposition="top center",
        textfont={"size": 9, "color": "#333"},
        hovertext=gene_hover,
        hoverinfo="text",
        name="Genes",
    )

    # GO term nodes -- color by domain
    domain_color = {**GO_DOMAIN_COLORS, "unknown": "#999"}
    go_nodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "go"]
    go_x = [pos[n][0] for n in go_nodes]
    go_y = [pos[n][1] for n in go_nodes]
    go_colors = [domain_color.get(G.nodes[n].get("domain", "unknown"), "#999") for n in go_nodes]
    go_sizes = [6 + 1.5 * G.degree(n) for n in go_nodes]
    go_hover = [
        f"<b>{n}</b><br>Domain: {G.nodes[n].get('domain', '?')}<br>Connections: {G.degree(n)}"
        for n in go_nodes
    ]

    go_trace = go.Scatter(
        x=go_x, y=go_y,
        mode="markers",
        marker={"size": go_sizes, "color": go_colors, "line": {"width": 0.5, "color": "#fff"}},
        hovertext=go_hover,
        hoverinfo="text",
        name="GO Terms",
    )

    fig = go.Figure(data=[edge_trace, go_trace, gene_trace])
    fig.update_layout(
        title="Gene-GO Term Network",
        template="plotly_white",
        showlegend=True,
        height=650,
        xaxis={"visible": False},
        yaxis={"visible": False},
        hovermode="closest",
    )
    return fig
