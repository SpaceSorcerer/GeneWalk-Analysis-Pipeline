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


def render_dashboard(df: pd.DataFrame, run_log: str | None = None):
    """Render the full interactive results dashboard.

    Parameters
    ----------
    df : pd.DataFrame
        GeneWalk results dataframe.
    run_log : str, optional
        Markdown-formatted run log to display in an expander.
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
        )

        padj_threshold = filter_cols[1].slider(
            "Significance threshold",
            0.001, 1.0, 0.05, 0.001,
            format="%.3f",
            help="Gene-GO pairs with adjusted p-value below this threshold are "
                 "considered significant.",
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

    # --- Available genes (shared across tabs) ---
    available_genes = (
        sorted(df["hgnc_symbol"].dropna().unique().tolist())
        if "hgnc_symbol" in df.columns
        else []
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
                use_container_width=True, theme=None,
            )
        with col_right:
            st.plotly_chart(
                go_domain_pie(df, padj_col=effective_padj_col,
                              padj_threshold=padj_threshold),
                use_container_width=True, theme=None,
            )

        if padj_col != "(none)":
            gene_summary = get_gene_summary(
                df, padj_col=effective_padj_col,
                padj_threshold=padj_threshold,
            )
            if not gene_summary.empty:
                st.plotly_chart(summary_bar(gene_summary),
                               use_container_width=True, theme=None)

        st.plotly_chart(
            pvalue_distribution(df, padj_col=effective_padj_col),
            use_container_width=True, theme=None,
        )

    # ---- Tab: Per-Gene Explorer -------------------------------------------
    with tab_gene:
        if available_genes:
            gene_col1, gene_col2 = st.columns([2, 1])
            selected_gene = gene_col1.selectbox("Select a gene", available_genes)
            top_n = gene_col2.slider("GO terms to show", 5, 50, 20)

            st.plotly_chart(
                gene_bar_chart(df, selected_gene,
                               padj_col=effective_padj_col, top_n=top_n),
                use_container_width=True, theme=None,
            )

            with st.expander(f"Raw data for {selected_gene}", expanded=False):
                gene_data = df[df["hgnc_symbol"] == selected_gene]
                if effective_padj_col in gene_data.columns:
                    gene_data = gene_data.sort_values(effective_padj_col)
                st.dataframe(gene_data, use_container_width=True, height=300)
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
                key="net_genes",
            )
            max_edges = net_col2.slider("Max edges", 50, 500, 200,
                                        key="net_edges")

            st.plotly_chart(
                gene_go_network(
                    df,
                    genes=net_genes or None,
                    padj_col=effective_padj_col,
                    padj_threshold=padj_threshold,
                    max_edges=max_edges,
                ),
                use_container_width=True, theme=None,
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
            )
            max_terms = hm_col2.slider("Max GO terms", 10, 100, 30)

            st.plotly_chart(
                gene_similarity_heatmap(
                    df,
                    genes=heatmap_genes or None,
                    padj_col=effective_padj_col,
                    padj_threshold=padj_threshold,
                    max_terms=max_terms,
                ),
                use_container_width=True, theme=None,
            )
        else:
            st.info("No gene symbols found in results.")

    # ---- Tab: Data Table --------------------------------------------------
    with tab_table:
        st.markdown(f"Showing **{len(filtered):,}** filtered rows")
        st.dataframe(filtered, use_container_width=True, height=500)

        csv_buf = io.StringIO()
        filtered.to_csv(csv_buf, index=False)
        st.download_button(
            "Download filtered results as CSV",
            csv_buf.getvalue(),
            file_name="genewalk_filtered_results.csv",
            mime="text/csv",
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
