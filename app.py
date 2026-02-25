"""GeneWalk Analysis Pipeline -- Streamlit App.

Launch with:  streamlit run app.py
"""

import io
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from genewalk_app.runner import (
    filter_results,
    find_results_csv,
    get_gene_summary,
    load_results,
    run_genewalk,
    save_gene_list,
)
from genewalk_app.styles import get_custom_css
from genewalk_app.visualizations import (
    gene_bar_chart,
    gene_go_network,
    gene_similarity_heatmap,
    go_domain_pie,
    pvalue_distribution,
    summary_bar,
    volcano_plot,
)

SAMPLE_GENES_PATH = Path(__file__).parent / "sample_data" / "sample_genes.txt"
SAMPLE_RESULTS_PATH = Path(__file__).parent / "sample_data" / "sample_genewalk_results.csv"

# ---------------------------------------------------------------------------
# Page config & custom styles
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GeneWalk Analysis Pipeline",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(get_custom_css(), unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Session state (must be initialized before any code that reads it)
# ---------------------------------------------------------------------------
if "results_df" not in st.session_state:
    st.session_state.results_df = None
if "run_log" not in st.session_state:
    st.session_state.run_log = None

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## GeneWalk Analysis")
    st.caption("Upload genes, configure, and explore results.")

    st.markdown('<p class="section-header">Gene List Input</p>', unsafe_allow_html=True)

    input_method = st.radio(
        "Choose input method",
        ["Upload file", "Paste genes", "Use sample data (with demo results)"],
        label_visibility="collapsed",
    )

    genes: list[str] = []

    if input_method == "Upload file":
        uploaded = st.file_uploader(
            "Upload gene list (.txt, one gene per line)",
            type=["txt", "csv", "tsv"],
            help="Plain text file with one gene identifier per line.",
        )
        if uploaded:
            content = uploaded.getvalue().decode("utf-8")
            genes = [g.strip() for g in content.splitlines() if g.strip()]
    elif input_method == "Paste genes":
        pasted = st.text_area(
            "Paste genes (one per line)",
            height=180,
            placeholder="BRAF\nMAP2K1\nMAPK1\n...",
        )
        genes = [g.strip() for g in pasted.splitlines() if g.strip()]
    else:
        if SAMPLE_GENES_PATH.exists():
            genes = [g.strip() for g in SAMPLE_GENES_PATH.read_text().splitlines() if g.strip()]
        if SAMPLE_RESULTS_PATH.exists() and st.session_state.results_df is None:
            st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
            st.session_state.run_log = "Demo results loaded from sample data."

    if genes:
        st.success(f"{len(genes)} genes loaded")
        with st.expander("Preview gene list", expanded=False):
            st.code("\n".join(genes[:50]) + ("\n..." if len(genes) > 50 else ""), language=None)

    st.markdown("---")
    st.markdown('<p class="section-header">Analysis Parameters</p>', unsafe_allow_html=True)

    id_type = st.selectbox(
        "Gene ID type",
        ["hgnc_symbol", "hgnc_id", "ensembl_id", "mgi_id", "rgd_id", "entrez"],
        index=0,
        help="Identifier type used in your gene list.",
    )

    col1, col2 = st.columns(2)
    nproc = col1.number_input("CPU cores", min_value=1, max_value=32, value=4,
                              help="Parallel workers for GeneWalk.")
    nreps_graph = col2.number_input("Graph reps", min_value=1, max_value=20, value=3,
                                    help="DeepWalk repetitions on the condition-specific network.")
    nreps_null = st.number_input("Null reps", min_value=1, max_value=20, value=3,
                                 help="DeepWalk repetitions on randomized networks.")
    alpha_fdr = st.slider("FDR threshold", 0.01, 1.0, 1.0, 0.01,
                          help="FDR threshold for GeneWalk. Use 1.0 to keep all results.")
    project_name = st.text_input("Project name", value="genewalk_analysis")

    st.markdown("---")
    run_clicked = st.button("Run GeneWalk", type="primary", disabled=len(genes) == 0,
                            use_container_width=True)

    st.markdown("---")
    st.markdown('<p class="section-header">Or Upload Existing Results</p>', unsafe_allow_html=True)
    uploaded_results = st.file_uploader(
        "Upload genewalk_results.csv",
        type=["csv"],
        key="results_upload",
        help="Upload a CSV from a previous GeneWalk run to skip the analysis step.",
    )

# ---------------------------------------------------------------------------
# Handle uploaded / run results
# ---------------------------------------------------------------------------
if uploaded_results is not None:
    st.session_state.results_df = pd.read_csv(uploaded_results)
    st.session_state.run_log = "Results loaded from uploaded CSV."

if run_clicked and genes:
    with st.status("Running GeneWalk analysis...", expanded=True) as status:
        st.write("Saving gene list...")
        tmp = Path(tempfile.mkdtemp())
        gene_file = save_gene_list(genes, tmp / "genes.txt")

        st.write(f"Starting GeneWalk with {nproc} cores, {nreps_graph} graph reps, "
                 f"{nreps_null} null reps...")
        result = run_genewalk(
            gene_file=gene_file,
            project=project_name,
            id_type=id_type,
            nproc=nproc,
            nreps_graph=nreps_graph,
            nreps_null=nreps_null,
            alpha_fdr=alpha_fdr,
        )

        st.session_state.run_log = (
            f"**Return code:** {result['return_code']}\n\n"
            f"**Output dir:** `{result['output_dir']}`\n\n"
            f"```\n{result['stderr'][-3000:] if result['stderr'] else '(no stderr)'}\n```"
        )

        if result["return_code"] == 0:
            csv_path = find_results_csv(result["output_dir"])
            if csv_path:
                st.session_state.results_df = load_results(csv_path)
                status.update(label="Analysis complete!", state="complete")
            else:
                status.update(label="Results CSV not found", state="error")
        else:
            status.update(label="GeneWalk exited with an error", state="error")

# ---------------------------------------------------------------------------
# Main content area
# ---------------------------------------------------------------------------

# Hero banner
st.markdown("""
<div class="hero-banner">
    <div class="hero-badge">Genomics Tool</div>
    <h1>GeneWalk Analysis Pipeline</h1>
    <p>Identify relevant gene functions for your biological context using network
    representation learning.
    Powered by <a href="https://github.com/churchmanlab/genewalk"
    style="color: #93c5fd; text-decoration: none; font-weight: 500;">GeneWalk</a>.</p>
</div>
""", unsafe_allow_html=True)

if st.session_state.results_df is None:
    # ----- Landing page -----
    col_a, col_b, col_c = st.columns(3, gap="medium")
    with col_a:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">1</div>
            <h3>Load Genes</h3>
            <p>Upload a text file with one gene per line, paste them directly,
            or try the built-in sample data (MAPK/ERK pathway).</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Try Sample Data", key="load_sample_btn", use_container_width=True,
                      type="primary"):
            if SAMPLE_RESULTS_PATH.exists():
                st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
                st.session_state.run_log = "Demo results loaded from sample data."
                st.rerun()
            else:
                st.error("Sample data files not found.")
    with col_b:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">2</div>
            <h3>Configure &amp; Run</h3>
            <p>Set your gene ID type and analysis parameters in the sidebar
            (&#9664;), then click <strong>Run GeneWalk</strong>. Or upload a
            previous <code>genewalk_results.csv</code>.</p>
        </div>
        """, unsafe_allow_html=True)
    with col_c:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">3</div>
            <h3>Explore Results</h3>
            <p>Interactive volcano plots, gene-GO networks, heatmaps, and
            per-gene bar charts &mdash; all filterable in real time with
            CSV export.</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("")
    st.markdown(
        '<div class="info-tip">Select <strong>Use sample data (with demo results)</strong> '
        "in the sidebar to instantly explore the app with pre-computed MAPK/ERK pathway data.</div>",
        unsafe_allow_html=True,
    )

    if st.session_state.run_log:
        with st.expander("Run log"):
            st.markdown(st.session_state.run_log)

    # Footer
    st.markdown(
        '<div class="app-footer">Built with '
        '<a href="https://streamlit.io">Streamlit</a> &middot; '
        'Analysis by <a href="https://github.com/churchmanlab/genewalk">GeneWalk</a> &middot; '
        'Ietswaart et al., <em>Genome Biology</em> 2021</div>',
        unsafe_allow_html=True,
    )
    st.stop()

# =========================================================================
# Results loaded -- show analysis dashboard
# =========================================================================
df = st.session_state.results_df

if st.session_state.run_log:
    with st.expander("Run log", expanded=False):
        st.markdown(st.session_state.run_log)

# --- Filters in a styled container ---
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
        help="Gene-GO pairs with adjusted p-value below this threshold are considered significant.",
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

# ---- Tab: Overview --------------------------------------------------------
with tab_overview:
    col_left, col_right = st.columns(2)
    with col_left:
        st.plotly_chart(
            volcano_plot(df, padj_col=effective_padj_col, padj_threshold=padj_threshold),
            use_container_width=True,
        )
    with col_right:
        st.plotly_chart(
            go_domain_pie(df, padj_col=effective_padj_col, padj_threshold=padj_threshold),
            use_container_width=True,
        )

    if padj_col != "(none)":
        summary = get_gene_summary(df, padj_col=effective_padj_col, padj_threshold=padj_threshold)
        if not summary.empty:
            st.plotly_chart(summary_bar(summary), use_container_width=True)

    st.plotly_chart(
        pvalue_distribution(df, padj_col=effective_padj_col),
        use_container_width=True,
    )

# ---- Tab: Per-Gene Explorer -----------------------------------------------
with tab_gene:
    if available_genes:
        gene_col1, gene_col2 = st.columns([2, 1])
        selected_gene = gene_col1.selectbox("Select a gene", available_genes)
        top_n = gene_col2.slider("GO terms to show", 5, 50, 20)

        st.plotly_chart(
            gene_bar_chart(df, selected_gene, padj_col=effective_padj_col, top_n=top_n),
            use_container_width=True,
        )

        with st.expander(f"Raw data for {selected_gene}", expanded=False):
            gene_data = df[df["hgnc_symbol"] == selected_gene]
            if effective_padj_col in gene_data.columns:
                gene_data = gene_data.sort_values(effective_padj_col)
            st.dataframe(gene_data, use_container_width=True, height=300)
    else:
        st.info("No gene symbols found in results.")

# ---- Tab: Network ---------------------------------------------------------
with tab_network:
    if available_genes:
        st.markdown(
            '<div class="info-tip">Select genes and adjust the edge limit to explore how '
            "your genes connect to GO terms. Larger networks take longer to render.</div>",
            unsafe_allow_html=True,
        )
        net_col1, net_col2 = st.columns([3, 1])
        net_genes = net_col1.multiselect(
            "Select genes for network",
            available_genes,
            default=available_genes[:6] if len(available_genes) >= 6 else available_genes,
            key="net_genes",
        )
        max_edges = net_col2.slider("Max edges", 50, 500, 200, key="net_edges")

        st.plotly_chart(
            gene_go_network(
                df,
                genes=net_genes or None,
                padj_col=effective_padj_col,
                padj_threshold=padj_threshold,
                max_edges=max_edges,
            ),
            use_container_width=True,
        )
        st.caption(
            "Red = genes | Blue = biological process | Orange = molecular function | "
            "Green = cellular component. Node size scales with connections."
        )
    else:
        st.info("No gene symbols found in results.")

# ---- Tab: Heatmap ---------------------------------------------------------
with tab_heatmap:
    if available_genes:
        hm_col1, hm_col2 = st.columns([3, 1])
        heatmap_genes = hm_col1.multiselect(
            "Select genes for heatmap (leave empty for all significant)",
            available_genes,
            default=available_genes[:10] if len(available_genes) >= 10 else available_genes,
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
            use_container_width=True,
        )
    else:
        st.info("No gene symbols found in results.")

# ---- Tab: Data Table -------------------------------------------------------
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
    'Analysis by <a href="https://github.com/churchmanlab/genewalk">GeneWalk</a> &middot; '
    'Ietswaart et al., <em>Genome Biology</em> 2021</div>',
    unsafe_allow_html=True,
)
