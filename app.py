"""GeneWalk Analysis Pipeline -- Streamlit App.

Launch with:  streamlit run app.py
"""

import io
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from genewalk_app.runner import (
    DEFAULT_BASE_DIR,
    filter_results,
    find_results_csv,
    get_gene_summary,
    load_results,
    run_genewalk,
    save_gene_list,
)
from genewalk_app.visualizations import (
    gene_bar_chart,
    gene_similarity_heatmap,
    go_domain_pie,
    pvalue_distribution,
    summary_bar,
    volcano_plot,
)

SAMPLE_GENES_PATH = Path(__file__).parent / "sample_data" / "sample_genes.txt"

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GeneWalk Analysis Pipeline",
    page_icon=":dna:",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Sidebar -- input & parameters
# ---------------------------------------------------------------------------
st.sidebar.title("GeneWalk Analysis")
st.sidebar.markdown("Upload a gene list, configure parameters, and run GeneWalk.")

st.sidebar.markdown("---")
st.sidebar.subheader("1. Gene List Input")

input_method = st.sidebar.radio(
    "Choose input method",
    ["Upload file", "Paste genes", "Use sample data"],
)

genes: list[str] = []

if input_method == "Upload file":
    uploaded = st.sidebar.file_uploader(
        "Upload gene list (.txt, one gene per line)",
        type=["txt", "csv", "tsv"],
    )
    if uploaded:
        content = uploaded.getvalue().decode("utf-8")
        genes = [g.strip() for g in content.splitlines() if g.strip()]
elif input_method == "Paste genes":
    pasted = st.sidebar.text_area(
        "Paste genes (one per line)",
        height=200,
        placeholder="BRAF\nMAP2K1\nMAPK1\n...",
    )
    genes = [g.strip() for g in pasted.splitlines() if g.strip()]
else:
    if SAMPLE_GENES_PATH.exists():
        genes = [g.strip() for g in SAMPLE_GENES_PATH.read_text().splitlines() if g.strip()]
        st.sidebar.info(f"Loaded {len(genes)} sample genes (MAPK/ERK pathway)")

if genes:
    st.sidebar.success(f"{len(genes)} genes loaded")
    with st.sidebar.expander("Preview genes"):
        st.code("\n".join(genes[:50]) + ("\n..." if len(genes) > 50 else ""))

st.sidebar.markdown("---")
st.sidebar.subheader("2. Parameters")

id_type = st.sidebar.selectbox(
    "Gene ID type",
    ["hgnc_symbol", "hgnc_id", "ensembl_id", "mgi_id", "rgd_id", "entrez"],
    index=0,
)

col1, col2 = st.sidebar.columns(2)
nproc = col1.number_input("CPU cores", min_value=1, max_value=32, value=4)
nreps_graph = col2.number_input("Graph reps", min_value=1, max_value=20, value=3)
nreps_null = st.sidebar.number_input("Null reps", min_value=1, max_value=20, value=3)
alpha_fdr = st.sidebar.slider("FDR threshold (alpha)", 0.01, 1.0, 1.0, 0.01,
                               help="FDR threshold for GeneWalk output. Use 1.0 to keep all results.")
project_name = st.sidebar.text_input("Project name", value="genewalk_analysis")

st.sidebar.markdown("---")

# ---------------------------------------------------------------------------
# Run GeneWalk
# ---------------------------------------------------------------------------
run_clicked = st.sidebar.button("Run GeneWalk", type="primary", disabled=len(genes) == 0)

# Allow uploading previous results instead of running
st.sidebar.markdown("---")
st.sidebar.subheader("Or: Upload Existing Results")
uploaded_results = st.sidebar.file_uploader(
    "Upload genewalk_results.csv",
    type=["csv"],
    key="results_upload",
)

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
if "results_df" not in st.session_state:
    st.session_state.results_df = None
if "run_log" not in st.session_state:
    st.session_state.run_log = None

# Handle uploaded results
if uploaded_results is not None:
    st.session_state.results_df = pd.read_csv(uploaded_results)
    st.session_state.run_log = "Results loaded from uploaded CSV."

# Handle run button
if run_clicked and genes:
    with st.spinner("Running GeneWalk -- this may take a while (typically 1-4 hours)..."):
        tmp = Path(tempfile.mkdtemp())
        gene_file = save_gene_list(genes, tmp / "genes.txt")

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
            else:
                st.error("GeneWalk completed but results CSV not found.")
        else:
            st.error("GeneWalk exited with an error. Check the run log below.")

# ---------------------------------------------------------------------------
# Main content area
# ---------------------------------------------------------------------------
st.title("GeneWalk Analysis Pipeline")
st.markdown(
    "An interactive tool for running [GeneWalk](https://github.com/churchmanlab/genewalk) "
    "and exploring gene-function associations."
)

if st.session_state.results_df is None:
    # Landing / instructions
    st.markdown("---")
    col_a, col_b, col_c = st.columns(3)
    with col_a:
        st.markdown("### Step 1: Load Genes")
        st.markdown(
            "Upload a text file with one gene per line, paste them directly, "
            "or try the built-in sample data (MAPK/ERK pathway)."
        )
    with col_b:
        st.markdown("### Step 2: Configure & Run")
        st.markdown(
            "Set your gene ID type and analysis parameters in the sidebar, "
            "then click **Run GeneWalk**. Alternatively, upload a previous "
            "`genewalk_results.csv`."
        )
    with col_c:
        st.markdown("### Step 3: Explore Results")
        st.markdown(
            "Interactive volcano plots, per-gene bar charts, heatmaps, and "
            "summary statistics -- all filterable in real time."
        )

    if st.session_state.run_log:
        with st.expander("Run log"):
            st.markdown(st.session_state.run_log)
    st.stop()

# ---------------------------------------------------------------------------
# Results loaded -- show analysis
# ---------------------------------------------------------------------------
df = st.session_state.results_df

if st.session_state.run_log:
    with st.expander("Run log", expanded=False):
        st.markdown(st.session_state.run_log)

# --- Filters ---
st.markdown("---")
st.subheader("Filter Results")

filter_cols = st.columns(3)

padj_col = filter_cols[0].selectbox(
    "P-value column",
    [c for c in ["gene_padj", "global_padj"] if c in df.columns] or ["(none)"],
)

padj_threshold = filter_cols[1].slider(
    "Significance threshold",
    0.001, 1.0, 0.05, 0.001,
    format="%.3f",
)

domain_options = ["All"] + sorted(df["go_domain"].dropna().unique().tolist()) if "go_domain" in df.columns else ["All"]
domain_filter = filter_cols[2].selectbox("GO domain", domain_options)

filtered = filter_results(
    df,
    padj_col=padj_col if padj_col != "(none)" else "gene_padj",
    padj_threshold=padj_threshold,
    go_domain=domain_filter if domain_filter != "All" else None,
)

# --- Summary metrics ---
st.markdown("---")
m1, m2, m3, m4 = st.columns(4)
m1.metric("Total gene-GO pairs", f"{len(df):,}")
m2.metric("Significant pairs", f"{len(filtered):,}")
m3.metric("Unique genes", f"{filtered['hgnc_symbol'].nunique() if 'hgnc_symbol' in filtered.columns else 'N/A'}")
m4.metric("Unique GO terms", f"{filtered['go_name'].nunique() if 'go_name' in filtered.columns else 'N/A'}")

# --- Tabs for visualizations ---
st.markdown("---")
tab_overview, tab_gene, tab_heatmap, tab_table = st.tabs(
    ["Overview", "Per-Gene Explorer", "Heatmap", "Data Table"]
)

# ---- Tab: Overview ----
with tab_overview:
    col_left, col_right = st.columns(2)

    with col_left:
        st.plotly_chart(
            volcano_plot(df, padj_col=padj_col if padj_col != "(none)" else "gene_padj", padj_threshold=padj_threshold),
            use_container_width=True,
        )

    with col_right:
        st.plotly_chart(
            go_domain_pie(df, padj_col=padj_col if padj_col != "(none)" else "gene_padj", padj_threshold=padj_threshold),
            use_container_width=True,
        )

    # Gene summary bar chart
    if padj_col != "(none)":
        summary = get_gene_summary(df, padj_col=padj_col, padj_threshold=padj_threshold)
        if not summary.empty:
            st.plotly_chart(summary_bar(summary), use_container_width=True)

    # P-value distribution
    st.plotly_chart(
        pvalue_distribution(df, padj_col=padj_col if padj_col != "(none)" else "gene_padj"),
        use_container_width=True,
    )

# ---- Tab: Per-Gene Explorer ----
with tab_gene:
    available_genes = sorted(df["hgnc_symbol"].dropna().unique().tolist()) if "hgnc_symbol" in df.columns else []

    if available_genes:
        selected_gene = st.selectbox("Select a gene", available_genes)
        top_n = st.slider("Number of GO terms to show", 5, 50, 20)

        st.plotly_chart(
            gene_bar_chart(df, selected_gene, padj_col=padj_col if padj_col != "(none)" else "gene_padj", top_n=top_n),
            use_container_width=True,
        )

        # Show the raw data for this gene
        gene_data = df[df["hgnc_symbol"] == selected_gene]
        if padj_col in gene_data.columns:
            gene_data = gene_data.sort_values(padj_col)
        st.dataframe(gene_data, use_container_width=True, height=300)
    else:
        st.info("No gene symbols found in results.")

# ---- Tab: Heatmap ----
with tab_heatmap:
    if available_genes:
        heatmap_genes = st.multiselect(
            "Select genes for heatmap (leave empty for all significant genes)",
            available_genes,
            default=available_genes[:10] if len(available_genes) >= 10 else available_genes,
        )
        max_terms = st.slider("Max GO terms to display", 10, 100, 30)

        st.plotly_chart(
            gene_similarity_heatmap(
                df,
                genes=heatmap_genes or None,
                padj_col=padj_col if padj_col != "(none)" else "gene_padj",
                padj_threshold=padj_threshold,
                max_terms=max_terms,
            ),
            use_container_width=True,
        )
    else:
        st.info("No gene symbols found in results.")

# ---- Tab: Data Table ----
with tab_table:
    st.markdown(f"Showing **{len(filtered):,}** filtered rows")
    st.dataframe(filtered, use_container_width=True, height=500)

    # Download button
    csv_buf = io.StringIO()
    filtered.to_csv(csv_buf, index=False)
    st.download_button(
        "Download filtered results as CSV",
        csv_buf.getvalue(),
        file_name="genewalk_filtered_results.csv",
        mime="text/csv",
    )
