from __future__ import annotations

"""Shared dashboard rendering for the GeneWalk Analysis Pipeline.

Both the web visualizer (app.py) and the desktop app (desktop.py) use this
module to render the interactive results dashboard, keeping the visualization
code in one place.
"""

import io

import pandas as pd
import streamlit as st

from genewalk_app.runner import filter_results, get_gene_summary
from genewalk_app.visualizations import (
    gene_bar_chart,
    gene_go_network,
    gene_similarity_heatmap,
    go_domain_pie,
    pvalue_distribution,
    summary_bar,
    volcano_plot,
)


def _sorted_genes(
    df: pd.DataFrame, method: str, padj_col: str, padj_threshold: float,
) -> list[str]:
    """Return gene list ordered by the chosen method."""
    if "hgnc_symbol" not in df.columns:
        return []

    if method == "Alphabetical":
        return sorted(df["hgnc_symbol"].dropna().unique().tolist())

    if method == "Input order":
        # Preserve the row order from the original CSV / GeneWalk output.
        seen: set[str] = set()
        ordered: list[str] = []
        for g in df["hgnc_symbol"].dropna():
            if g not in seen:
                seen.add(g)
                ordered.append(g)
        return ordered

    if method == "Best p-value" and padj_col in df.columns:
        best = (
            df.dropna(subset=[padj_col])
            .groupby("hgnc_symbol")[padj_col]
            .min()
            .sort_values()
        )
        return best.index.tolist()

    if method == "Most significant GO terms":
        sig = df[df[padj_col] <= padj_threshold] if padj_col in df.columns else df
        counts = sig.groupby("hgnc_symbol").size().sort_values(ascending=False)
        return counts.index.tolist()

    if method == "Mean similarity" and "sim" in df.columns:
        mean_sim = (
            df.dropna(subset=["sim"])
            .groupby("hgnc_symbol")["sim"]
            .mean()
            .sort_values(ascending=False)
        )
        return mean_sim.index.tolist()

    return sorted(df["hgnc_symbol"].dropna().unique().tolist())


def render_dashboard(
    df: pd.DataFrame,
    run_log: str | None = None,
    key_prefix: str = "",
):
    """Render the full interactive results dashboard.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results dataframe.
    run_log : str, optional
        Markdown-formatted run log to display in an expander.
    key_prefix : str
        Prefix for all Streamlit widget keys, allowing multiple dashboard
        instances on the same page (e.g. up-regulated vs down-regulated).
    """
    if run_log:
        with st.expander("Run log", expanded=False):
            st.markdown(run_log)

    # --- Filters ---
    st.markdown('<p class="section-header">Filters</p>', unsafe_allow_html=True)

    with st.container():
        filter_cols = st.columns([1, 1, 1])

        padj_col = filter_cols[0].selectbox(
            "P-value column",
            [c for c in ["gene_padj", "global_padj"] if c in df.columns] or ["(none)"],
            help="Which adjusted p-value to use for filtering and visualization.",
            key=f"{key_prefix}padj_col",
        )

        padj_threshold = filter_cols[1].slider(
            "Significance threshold",
            0.001, 1.0, 0.05, 0.001,
            format="%.3f",
            help="Gene-GO pairs with adjusted p-value below this threshold are "
                 "considered significant.",
            key=f"{key_prefix}padj_threshold",
        )

        domain_options = (
            ["All"] + sorted(df["go_domain"].dropna().unique().tolist())
            if "go_domain" in df.columns
            else ["All"]
        )
        domain_filter = filter_cols[2].selectbox(
            "GO domain",
            domain_options,
            help="Filter to a specific Gene Ontology domain.",
            key=f"{key_prefix}domain_filter",
        )

    effective_padj_col = padj_col if padj_col != "(none)" else "gene_padj"

    filtered = filter_results(
        df,
        padj_col=effective_padj_col,
        padj_threshold=padj_threshold,
        go_domain=domain_filter if domain_filter != "All" else None,
    )

    # --- Summary metrics ---
    st.markdown("")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Total Pairs", f"{len(df):,}")
    m2.metric("Significant", f"{len(filtered):,}")
    m3.metric(
        "Unique Genes",
        f"{filtered['hgnc_symbol'].nunique() if 'hgnc_symbol' in filtered.columns else 'N/A'}",
    )
    m4.metric(
        "GO Terms",
        f"{filtered['go_name'].nunique() if 'go_name' in filtered.columns else 'N/A'}",
    )

    # --- Interpretive guide ---
    with st.expander("How to interpret these results", expanded=False):
        st.markdown("""
<div class="guide-section">

<h4>What GeneWalk measures</h4>
<p>GeneWalk does <strong>not</strong> measure how much a gene changes
(fold-change, expression level, etc.). Instead, it asks: <em>"Given the
genes you provided, which GO functions are most relevant to this specific
biological context?"</em></p>

<p>It works by building a network of your genes and all their known GO
annotations, then using DeepWalk (a graph-embedding algorithm) to learn
which gene&ndash;function associations are <strong>specific to your gene
set</strong> versus generic background knowledge.</p>

<h4>Key columns explained</h4>
<ul>
<li><strong>Similarity score (sim)</strong> &mdash; How closely a gene
and a GO term are embedded in the context-specific network. Higher =
stronger association <em>in this context</em>. A gene may have a known
GO annotation that scores low here because it&rsquo;s not relevant to
your particular set of genes.</li>
<li><strong>gene_padj</strong> &mdash; Per-gene adjusted p-value. Tests
whether this gene&ndash;GO similarity is higher than expected by chance,
correcting for multiple GO terms tested within each gene. Good for
exploring what&rsquo;s relevant <em>for a specific gene</em>.</li>
<li><strong>global_padj</strong> &mdash; Global adjusted p-value. Same
test but corrected across <em>all genes and all GO terms</em>
simultaneously. Much more stringent &mdash; use this for high-confidence
hits you would report in a paper.</li>
</ul>

<h4>What does "most significant" mean here?</h4>
<p>A gene with many significant GO terms has many functions that are
specifically relevant to your biological context &mdash; it&rsquo;s a
<strong>functional hub</strong> in your gene set. This does NOT mean
it&rsquo;s the most differentially expressed or the most important
biologically. It means GeneWalk found strong, context-specific evidence
connecting it to multiple functions.</p>

<h4>Important caveats</h4>
<ul>
<li><strong>No directionality:</strong> GeneWalk does not know whether
your genes are up- or down-regulated. It treats all input genes equally.
If directionality matters, consider running separate analyses for up-
and down-regulated gene lists.</li>
<li><strong>Annotation bias:</strong> Well-studied genes (e.g. TP53,
EGFR) have more GO annotations and may appear more "significant" simply
because more is known about them. Sparse annotations for less-studied
genes can lead to fewer significant results, not necessarily less
biological importance.</li>
<li><strong>Context is your gene list:</strong> Results change depending
on which genes you include. The same gene may have different significant
GO terms when analyzed with different gene sets.</li>
</ul>
</div>
        """, unsafe_allow_html=True)

    # --- Gene sort order (shared across tabs) ---
    sort_methods = [
        "Alphabetical",
        "Input order",
        "Best p-value",
        "Most significant GO terms",
        "Mean similarity",
    ]
    gene_sort = st.selectbox(
        "Sort genes by",
        sort_methods,
        index=0,
        help="Controls gene ordering in the Per-Gene Explorer, Network, and "
             "Heatmap selectors. 'Input order' preserves the order from your "
             "original gene list / GeneWalk output.",
        key=f"{key_prefix}gene_sort",
    )
    available_genes = _sorted_genes(
        df, gene_sort, effective_padj_col, padj_threshold,
    )

    # --- Visualization tabs ---
    st.markdown("")
    tab_overview, tab_gene, tab_network, tab_heatmap, tab_table = st.tabs(
        ["Overview", "Per-Gene Explorer", "Network", "Heatmap", "Data Table"]
    )

    # ---- Tab: Overview ----------------------------------------------------
    with tab_overview:
        col_left, col_right = st.columns(2)
        with col_left:
            st.plotly_chart(
                volcano_plot(df, padj_col=effective_padj_col,
                             padj_threshold=padj_threshold),
                width="stretch",
                key=f"{key_prefix}volcano",
            )
        with col_right:
            st.plotly_chart(
                go_domain_pie(df, padj_col=effective_padj_col,
                              padj_threshold=padj_threshold),
                width="stretch",
                key=f"{key_prefix}domain_pie",
            )

        if padj_col != "(none)":
            gene_summary = get_gene_summary(
                df, padj_col=effective_padj_col,
                padj_threshold=padj_threshold,
            )
            if not gene_summary.empty:
                st.plotly_chart(summary_bar(gene_summary),
                               width="stretch",
                               key=f"{key_prefix}summary_bar")

        st.plotly_chart(
            pvalue_distribution(df, padj_col=effective_padj_col),
            width="stretch",
            key=f"{key_prefix}pval_dist",
        )

    # ---- Tab: Per-Gene Explorer -------------------------------------------
    with tab_gene:
        if available_genes:
            gene_col1, gene_col2 = st.columns([2, 1])
            selected_gene = gene_col1.selectbox(
                "Select a gene", available_genes,
                key=f"{key_prefix}selected_gene",
            )
            top_n = gene_col2.slider(
                "GO terms to show", 5, 50, 20,
                key=f"{key_prefix}top_n",
            )

            st.plotly_chart(
                gene_bar_chart(df, selected_gene,
                               padj_col=effective_padj_col, top_n=top_n),
                width="stretch",
                key=f"{key_prefix}gene_bar",
            )

            with st.expander(f"Raw data for {selected_gene}", expanded=False):
                gene_data = df[df["hgnc_symbol"] == selected_gene]
                if effective_padj_col in gene_data.columns:
                    gene_data = gene_data.sort_values(effective_padj_col)
                st.dataframe(gene_data, width="stretch", height=300)
        else:
            st.info("No gene symbols found in results.")

    # ---- Tab: Network -----------------------------------------------------
    with tab_network:
        if available_genes:
            st.markdown(
                '<div class="info-tip">Select genes and adjust the edge limit '
                "to explore how your genes connect to GO terms. Larger networks "
                "take longer to render.</div>",
                unsafe_allow_html=True,
            )
            net_col1, net_col2 = st.columns([3, 1])
            net_genes = net_col1.multiselect(
                "Select genes for network",
                available_genes,
                default=(available_genes[:6]
                         if len(available_genes) >= 6 else available_genes),
                key=f"{key_prefix}net_genes",
            )
            max_edges = net_col2.slider("Max edges", 50, 500, 200,
                                        key=f"{key_prefix}net_edges")

            st.plotly_chart(
                gene_go_network(
                    df,
                    genes=net_genes or None,
                    padj_col=effective_padj_col,
                    padj_threshold=padj_threshold,
                    max_edges=max_edges,
                ),
                width="stretch",
                key=f"{key_prefix}network",
            )
            st.caption(
                "Red = genes | Blue = biological process | "
                "Orange = molecular function | Green = cellular component. "
                "Node size scales with connections."
            )
        else:
            st.info("No gene symbols found in results.")

    # ---- Tab: Heatmap -----------------------------------------------------
    with tab_heatmap:
        if available_genes:
            hm_col1, hm_col2 = st.columns([3, 1])
            heatmap_genes = hm_col1.multiselect(
                "Select genes for heatmap (leave empty for all significant)",
                available_genes,
                default=(available_genes[:10]
                         if len(available_genes) >= 10 else available_genes),
                key=f"{key_prefix}heatmap_genes",
            )
            max_terms = hm_col2.slider("Max GO terms", 10, 100, 30,
                                       key=f"{key_prefix}max_terms")

            st.plotly_chart(
                gene_similarity_heatmap(
                    df,
                    genes=heatmap_genes or None,
                    padj_col=effective_padj_col,
                    padj_threshold=padj_threshold,
                    max_terms=max_terms,
                ),
                width="stretch",
                key=f"{key_prefix}heatmap",
            )
        else:
            st.info("No gene symbols found in results.")

    # ---- Tab: Data Table --------------------------------------------------
    with tab_table:
        st.markdown(f"Showing **{len(filtered):,}** filtered rows")
        st.dataframe(filtered, width="stretch", height=500)

        csv_buf = io.StringIO()
        filtered.to_csv(csv_buf, index=False)
        st.download_button(
            "Download filtered results as CSV",
            csv_buf.getvalue(),
            file_name="genewalk_filtered_results.csv",
            mime="text/csv",
            key=f"{key_prefix}download_csv",
        )

    # Footer
    st.markdown(
        '<div class="app-footer">Built with '
        '<a href="https://streamlit.io">Streamlit</a> &middot; '
        'Analysis by <a href="https://github.com/churchmanlab/genewalk">'
        "GeneWalk</a> &middot; "
        "Ietswaart et al., <em>Genome Biology</em> 2021</div>",
        unsafe_allow_html=True,
    )
