"""GeneWalk Results Visualizer -- Web App.

Online visualizer for exploring GeneWalk results. Users upload a
``genewalk_results.csv`` from a local GeneWalk run (or the desktop app)
and explore interactive charts, filters, and data tables.

Launch with:  streamlit run app.py
"""

from pathlib import Path

import pandas as pd
import streamlit as st

from genewalk_app.dashboard import render_dashboard
from genewalk_app.styles import get_custom_css

SAMPLE_RESULTS_PATH = Path(__file__).parent / "sample_data" / "sample_genewalk_results.csv"

# ---------------------------------------------------------------------------
# Page config & custom styles
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GeneWalk Results Visualizer",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(get_custom_css(), unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
if "results_df" not in st.session_state:
    st.session_state.results_df = None
if "run_log" not in st.session_state:
    st.session_state.run_log = None

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## GeneWalk Visualizer")
    st.caption("Upload GeneWalk results and explore interactively.")

    st.markdown('<p class="section-header">Load Results</p>',
                unsafe_allow_html=True)

    input_method = st.radio(
        "Choose input method",
        ["Upload results CSV", "Use sample data (MAPK/ERK demo)"],
        label_visibility="collapsed",
    )

    if input_method == "Upload results CSV":
        uploaded_results = st.file_uploader(
            "Upload genewalk_results.csv",
            type=["csv"],
            help="Upload a CSV from a previous GeneWalk run. "
                 "This is the ``genewalk_results.csv`` file produced by "
                 "GeneWalk or the GeneWalk Desktop app.",
        )
        if uploaded_results is not None:
            st.session_state.results_df = pd.read_csv(uploaded_results)
            st.session_state.run_log = "Results loaded from uploaded CSV."
    else:
        if SAMPLE_RESULTS_PATH.exists() and st.session_state.results_df is None:
            st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
            st.session_state.run_log = "Demo results loaded from sample data."

    st.markdown("---")
    st.markdown(
        '<div class="info-tip">'
        "<strong>Need to run GeneWalk?</strong><br>"
        'Download the <a href="https://github.com/SpaceSorcerer/'
        'GeneWalk-Analysis-Pipeline/releases">GeneWalk Desktop app</a> '
        "to run analyses on your own machine, then upload the results here."
        "</div>",
        unsafe_allow_html=True,
    )

# ---------------------------------------------------------------------------
# Main content area
# ---------------------------------------------------------------------------

# Hero banner
st.markdown("""
<div class="hero-banner">
    <div class="hero-badge">Results Visualizer</div>
    <h1>GeneWalk Analysis Pipeline</h1>
    <p>Explore GeneWalk results with interactive volcano plots, gene-GO
    networks, heatmaps, and more.
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
            <h3>Run GeneWalk</h3>
            <p>Download the <strong>GeneWalk Desktop app</strong> or run
            GeneWalk from the command line on your own machine. This
            produces a <code>genewalk_results.csv</code> file.</p>
        </div>
        """, unsafe_allow_html=True)
    with col_b:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">2</div>
            <h3>Upload Results</h3>
            <p>Upload your <code>genewalk_results.csv</code> using the
            sidebar on the left. Or try the built-in sample data to see
            the visualizer in action.</p>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Try Sample Data", key="load_sample_btn",
                      use_container_width=True, type="primary"):
            if SAMPLE_RESULTS_PATH.exists():
                st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
                st.session_state.run_log = "Demo results loaded from sample data."
                st.rerun()
            else:
                st.error("Sample data files not found.")
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
        '<div class="info-tip">Select <strong>Use sample data (MAPK/ERK demo)'
        "</strong> in the sidebar to instantly explore the app with "
        "pre-computed MAPK/ERK pathway data.</div>",
        unsafe_allow_html=True,
    )

    # ------------------------------------------------------------------
    # Getting Started Guide
    # ------------------------------------------------------------------
    st.markdown('<p class="guide-header">Getting Started Guide</p>',
                unsafe_allow_html=True)

    with st.expander("How to use this visualizer", expanded=False):
        st.markdown("""
<div class="guide-section">
<h4>Quick walkthrough</h4>
<ol>
<li><strong>Get your results</strong> &mdash; Run GeneWalk using the
<strong>Desktop app</strong> (download from
<a href="https://github.com/SpaceSorcerer/GeneWalk-Analysis-Pipeline/releases">
GitHub Releases</a>) or via the command line (<code>pip install genewalk</code>).
This produces a <code>genewalk_results.csv</code> file.</li>
<li><strong>Upload here</strong> &mdash; In the sidebar, select
<em>Upload results CSV</em> and drag in your file.</li>
<li><strong>Explore</strong> &mdash; The dashboard appears with five
interactive tabs:</li>
</ol>
<ul>
<li><strong>Overview</strong> &mdash; Volcano plot, GO domain breakdown,
per-gene summary bar chart, and p-value distribution.</li>
<li><strong>Per-Gene Explorer</strong> &mdash; Select any gene and see its top
GO term associations ranked by significance.</li>
<li><strong>Network</strong> &mdash; Interactive gene&ndash;GO bipartite
network. Select genes and adjust max edges to explore connections.</li>
<li><strong>Heatmap</strong> &mdash; Gene &times; GO term similarity heatmap
showing which genes share functional annotations.</li>
<li><strong>Data Table</strong> &mdash; Full filterable table with CSV
download.</li>
</ul>
<p>All visualizations respond to the <strong>Filters</strong> bar at the top of
the dashboard (p-value column, significance threshold, and GO domain).</p>
</div>
        """, unsafe_allow_html=True)

    with st.expander("What is GeneWalk?", expanded=False):
        st.markdown("""
<div class="guide-section">
<p><a href="https://github.com/churchmanlab/genewalk">GeneWalk</a> identifies
relevant gene functions for a biological context using network representation
learning. Given a list of genes (e.g. from a differential expression analysis),
GeneWalk:</p>
<ol>
<li>Builds a gene&ndash;GO term network from known annotations and
Pathway Commons interactions.</li>
<li>Uses DeepWalk random walks to learn vector representations of genes and
GO terms in the network.</li>
<li>Computes similarity scores and statistical significance for each
gene&ndash;GO term pair.</li>
</ol>
<p>The result is a ranked list of gene-function associations with p-values,
which this visualizer helps you explore.</p>
<p><strong>Citation:</strong> Ietswaart et al., <em>GeneWalk identifies relevant
gene functions for a biological context using network representation learning</em>,
Genome Biology 2021.</p>
</div>
        """, unsafe_allow_html=True)

    with st.expander("Expected CSV format", expanded=False):
        st.markdown("""
<div class="guide-section">
<p>This visualizer expects a <code>genewalk_results.csv</code> with these
columns:</p>
<table class="param-table">
<tr><th>Column</th><th>Description</th></tr>
<tr><td><code>hgnc_symbol</code></td><td>Gene symbol (e.g. BRAF)</td></tr>
<tr><td><code>go_name</code></td><td>GO term name</td></tr>
<tr><td><code>go_id</code></td><td>GO term ID (e.g. GO:0000165)</td></tr>
<tr><td><code>go_domain</code></td><td>biological_process, molecular_function,
or cellular_component</td></tr>
<tr><td><code>sim</code></td><td>Similarity score (0&ndash;1)</td></tr>
<tr><td><code>gene_padj</code></td><td>Gene-level adjusted p-value</td></tr>
<tr><td><code>global_padj</code></td><td>Global adjusted p-value</td></tr>
</table>
<p>This is the standard output format from GeneWalk. If you ran GeneWalk with
default settings, your CSV will already have the correct columns.</p>
</div>
        """, unsafe_allow_html=True)

    # Footer
    st.markdown(
        '<div class="app-footer">Built with '
        '<a href="https://streamlit.io">Streamlit</a> &middot; '
        'Analysis by <a href="https://github.com/churchmanlab/genewalk">'
        "GeneWalk</a> &middot; "
        "Ietswaart et al., <em>Genome Biology</em> 2021</div>",
        unsafe_allow_html=True,
    )
    st.stop()

# =========================================================================
# Results loaded -- show analysis dashboard
# =========================================================================
render_dashboard(st.session_state.results_df, st.session_state.run_log)
