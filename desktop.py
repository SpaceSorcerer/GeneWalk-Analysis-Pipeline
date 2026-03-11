"""GeneWalk App -- Desktop.

Full desktop application that can run GeneWalk analyses and visualize results.
This is the entry point bundled by PyInstaller into an executable.

Launch with:  streamlit run desktop.py
"""

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

# Signal to the launcher that a browser client has connected so it can
# skip opening a duplicate tab (handles browser session restore).
Path(tempfile.gettempdir(), ".genewalk_client_connected").write_text("")

from genewalk_app.dashboard import render_dashboard
from genewalk_app.runner import (
    _sanitize_project_name,
    find_results_csv,
    is_genewalk_available,
    load_results,
    run_genewalk,
    save_gene_list,
)
from genewalk_app.styles import get_custom_css

SAMPLE_GENES_PATH = Path(__file__).parent / "sample_data" / "sample_genes.txt"
SAMPLE_RESULTS_PATH = (
    Path(__file__).parent / "sample_data" / "sample_genewalk_results.csv"
)

# ---------------------------------------------------------------------------
# Page config & custom styles
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GeneWalk App",
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
    st.markdown("## GeneWalk Analysis")
    st.caption("Upload genes, configure, and explore results.")

    st.markdown('<p class="section-header">Gene List Input</p>',
                unsafe_allow_html=True)

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
            try:
                content = uploaded.getvalue().decode("utf-8")
            except UnicodeDecodeError:
                content = uploaded.getvalue().decode("latin-1")
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
            genes = [
                g.strip()
                for g in SAMPLE_GENES_PATH.read_text().splitlines()
                if g.strip()
            ]
        if SAMPLE_RESULTS_PATH.exists() and st.session_state.results_df is None:
            st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
            st.session_state.run_log = "Demo results loaded from sample data."

    if genes:
        st.success(f"{len(genes)} genes loaded")
        with st.expander("Preview gene list", expanded=False):
            st.code(
                "\n".join(genes[:50])
                + ("\n..." if len(genes) > 50 else ""),
                language=None,
            )

    st.markdown("---")
    st.markdown('<p class="section-header">Analysis Parameters</p>',
                unsafe_allow_html=True)

    id_type = st.selectbox(
        "Gene ID type",
        ["hgnc_symbol", "hgnc_id", "ensembl_id", "mgi_id", "rgd_id",
         "entrez_human", "entrez_mouse"],
        index=0,
        help="Must match your gene list. Use hgnc_symbol for human gene names "
             "(e.g. BRAF, KRAS), ensembl_id for Ensembl IDs, mgi_id for mouse, "
             "rgd_id for rat, entrez_human/entrez_mouse for Entrez Gene IDs.",
    )

    col1, col2 = st.columns(2)
    nproc = col1.number_input(
        "CPU cores", min_value=1, max_value=32, value=4,
        help="Parallel workers. Use 2-4 for laptops, up to 8-16 "
             "for workstations. Don't exceed your physical core count. "
             "More cores = faster but higher memory usage.",
    )
    nreps_graph = col2.number_input(
        "Graph reps", min_value=1, max_value=20, value=3,
        help="DeepWalk repetitions on the gene-GO network. "
             "Use 3 for quick testing, 10-15 for publication-quality "
             "results. Higher = more stable but slower.",
    )
    nreps_null = st.number_input(
        "Null reps", min_value=1, max_value=20, value=3,
        help="DeepWalk repetitions on randomized networks for "
             "statistical background. Should generally match Graph reps. "
             "Use 3 for testing, 10-15 for final analyses.",
    )
    alpha_fdr = st.slider(
        "FDR threshold", 0.01, 1.0, 1.0, 0.01,
        help="False discovery rate cutoff. Set to 1.0 (default) to keep "
             "all results and filter interactively in the dashboard. "
             "Set to 0.05 or 0.1 to pre-filter during the GeneWalk run.",
    )
    project_name = st.text_input(
        "Project name", value="genewalk_analysis",
        help="Names the output directory. Use something descriptive "
             "like 'mapk_screen' or 'patient_deg_analysis'.",
    )

    output_folder = st.text_input(
        "Output folder",
        value="",
        help="Where to save GeneWalk results. Leave blank to use the system "
             "temp directory. Examples: C:\\Users\\you\\genewalk_output or "
             "~/genewalk_output",
        placeholder="Leave blank for default (system temp)",
    )

    st.markdown("---")
    genewalk_installed = is_genewalk_available()
    if not genewalk_installed:
        st.warning(
            "GeneWalk is not installed (or not found on PATH). "
            "Install it with `pip install genewalk`. "
            "You can still upload a previous results CSV below.",
            icon="\u26a0\ufe0f",
        )
    run_clicked = st.button(
        "Run GeneWalk", type="primary",
        disabled=len(genes) == 0 or not genewalk_installed,
        use_container_width=True,
    )

    st.markdown("---")
    st.markdown('<p class="section-header">Or Upload Existing Results</p>',
                unsafe_allow_html=True)
    uploaded_results = st.file_uploader(
        "Upload genewalk_results.csv",
        type=["csv"],
        key="results_upload",
        help="Upload a CSV from a previous GeneWalk run to skip the "
             "analysis step.",
    )

# ---------------------------------------------------------------------------
# Handle uploaded / run results
# ---------------------------------------------------------------------------
if uploaded_results is not None:
    try:
        st.session_state.results_df = pd.read_csv(uploaded_results)
        st.session_state.run_log = "Results loaded from uploaded CSV."
    except Exception as exc:
        st.error(f"Could not read uploaded CSV: {exc}")

if run_clicked and genes:
    with st.status("Running GeneWalk analysis...", expanded=True) as status:
        st.write("Saving gene list...")
        tmp = Path(tempfile.mkdtemp())
        gene_file = save_gene_list(genes, tmp / "genes.txt")

        base_folder = (
            Path(output_folder).expanduser() if output_folder.strip() else None
        )
        st.write(
            f"Starting GeneWalk with {nproc} cores, {nreps_graph} graph reps, "
            f"{nreps_null} null reps..."
        )
        progress_placeholder = st.empty()

        def _update_progress(msg: str) -> None:
            progress_placeholder.write(msg)

        safe_project = _sanitize_project_name(project_name)
        result = run_genewalk(
            gene_file=gene_file,
            project=safe_project,
            id_type=id_type,
            nproc=nproc,
            nreps_graph=nreps_graph,
            nreps_null=nreps_null,
            alpha_fdr=alpha_fdr,
            base_folder=base_folder,
            on_progress=_update_progress,
        )

        st.session_state.run_log = (
            f"**Return code:** {result['return_code']}\n\n"
            f"**Output dir:** `{result['output_dir']}`\n\n"
            f"```\n{result['stderr'][-3000:] if result['stderr'] else '(no stderr)'}\n```"
        )

        csv_path = find_results_csv(result["output_dir"])
        if result["return_code"] == 0:
            if csv_path:
                try:
                    st.session_state.results_df = load_results(csv_path)
                except ValueError as exc:
                    st.error(str(exc))
                    status.update(label="Results CSV is invalid", state="error")
                    csv_path = None  # prevent duplicate handling below
                if csv_path:
                    status.update(label="Analysis complete!", state="complete")
            else:
                status.update(label="Results CSV not found", state="error")
        else:
            if csv_path:
                try:
                    st.session_state.results_df = load_results(csv_path)
                except ValueError as exc:
                    st.error(str(exc))
                    csv_path = None
                if csv_path:
                    st.session_state.run_log += (
                        "\n\n**Note:** GeneWalk exited with an error (likely "
                        "during HTML report generation), but the results CSV "
                        "was found and loaded successfully."
                    )
                    status.update(
                        label="Analysis complete (with warnings)", state="complete"
                    )
            if not csv_path and result["return_code"] != 0:
                status.update(
                    label="GeneWalk exited with an error", state="error"
                )

# ---------------------------------------------------------------------------
# Main content area
# ---------------------------------------------------------------------------

# Hero banner
st.markdown("""
<div class="hero-banner">
    <div class="hero-badge">Desktop App</div>
    <h1>GeneWalk App</h1>
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
        if st.button("Try Sample Data", key="load_sample_btn",
                      use_container_width=True, type="primary"):
            if SAMPLE_RESULTS_PATH.exists():
                st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
                st.session_state.run_log = (
                    "Demo results loaded from sample data."
                )
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
        '<div class="info-tip">Select <strong>Use sample data (with demo '
        "results)</strong> in the sidebar to instantly explore the app with "
        "pre-computed MAPK/ERK pathway data.</div>",
        unsafe_allow_html=True,
    )

    # ------------------------------------------------------------------
    # Getting Started Guide
    # ------------------------------------------------------------------
    st.markdown('<p class="guide-header">Getting Started Guide</p>',
                unsafe_allow_html=True)

    with st.expander("How to use this app", expanded=False):
        st.markdown("""
<div class="guide-section">
<h4>Quick walkthrough</h4>
<ol>
<li><strong>Load your genes</strong> &mdash; In the sidebar on the left, choose
an input method: upload a <code>.txt</code> file (one gene per line), paste gene
names directly, or select <em>Use sample data</em> to try the built-in
MAPK/ERK pathway demo.</li>
<li><strong>Configure parameters</strong> &mdash; Pick the correct gene ID type,
adjust CPU cores &amp; repetition counts, and set the FDR threshold (see the
<em>Understanding the parameters</em> section below for guidance).</li>
<li><strong>Run GeneWalk</strong> &mdash; Click the blue <em>Run GeneWalk</em>
button. The analysis runs locally on your machine and progress is shown in real
time. For a quick test, lower the repetition counts.</li>
<li><strong>Explore results</strong> &mdash; Once the run completes (or you
upload a previous <code>genewalk_results.csv</code>), the dashboard appears with
five interactive tabs:</li>
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

    with st.expander("Understanding the parameters", expanded=False):
        st.markdown("""
<div class="guide-section">
<table class="param-table">
<tr>
    <th>Parameter</th><th>What it does</th><th>Recommendation</th>
</tr>
<tr>
    <td><strong>Gene ID type</strong></td>
    <td>Tells GeneWalk which identifier system your gene list uses. Must match
    exactly or genes won't be recognized.</td>
    <td><span class="recommend-badge">hgnc_symbol</span> for human genes
    (e.g. <code>BRAF</code>, <code>KRAS</code>). Use <code>ensembl_id</code>
    for Ensembl IDs, <code>mgi_id</code> for mouse, <code>rgd_id</code>
    for rat.</td>
</tr>
<tr>
    <td><strong>CPU cores</strong></td>
    <td>Number of parallel workers for the DeepWalk random-walk step. More
    cores = faster analysis, but uses more CPU and memory.</td>
    <td><span class="recommend-badge">2&ndash;4</span> for laptops/small
    machines. Up to <strong>8&ndash;16</strong> if you have a workstation.
    Don't exceed your physical core count.</td>
</tr>
<tr>
    <td><strong>Graph reps</strong></td>
    <td>Number of DeepWalk repetitions on the condition-specific gene&ndash;GO
    network. More reps = more stable embeddings but longer runtime.</td>
    <td><span class="recommend-badge">3</span> for a quick exploratory run.
    <strong>10&ndash;15</strong> for publication-quality results.</td>
</tr>
<tr>
    <td><strong>Null reps</strong></td>
    <td>DeepWalk repetitions on randomized (null) networks, used to build the
    statistical background. More reps = tighter p-values.</td>
    <td><span class="recommend-badge">3</span> for quick testing.
    <strong>10&ndash;15</strong> for final analyses. Should generally match
    Graph reps.</td>
</tr>
<tr>
    <td><strong>FDR threshold</strong></td>
    <td>False discovery rate cutoff applied during the GeneWalk run.
    Gene&ndash;GO associations above this threshold may be discarded
    internally.</td>
    <td><span class="recommend-badge">1.0</span> (keep everything) so you can
    filter interactively in the dashboard. Set to <strong>0.05</strong> or
    <strong>0.1</strong> if you want GeneWalk to pre-filter for you.</td>
</tr>
<tr>
    <td><strong>Project name</strong></td>
    <td>Names the output directory. Useful when running multiple analyses to
    keep results organized.</td>
    <td>Something descriptive, e.g. <code>mapk_erk_screen</code> or
    <code>patient_cohort_deg</code>.</td>
</tr>
</table>

<h4>Performance &amp; resource allocation</h4>
<ul>
<li><strong>RAM:</strong> GeneWalk typically uses <strong>2&ndash;4 GB</strong>
for a list of ~20 genes. For 100+ genes, expect <strong>4&ndash;8 GB</strong>.
Make sure your machine has enough free memory before starting a run.</li>
<li><strong>Runtime:</strong> A 20-gene list with 3 graph/null reps on 4 cores
takes roughly <strong>5&ndash;15 minutes</strong>. With 10 reps each it can
take <strong>30&ndash;60+ minutes</strong>. Very large gene lists (200+) may
need hours.</li>
<li><strong>Disk:</strong> Each run produces ~5&ndash;50 MB of output files
depending on gene count and reps.</li>
<li><strong>Tip:</strong> Start with a small gene list and low reps (3/3) to
verify everything works, then increase for your real analysis.</li>
</ul>
</div>
        """, unsafe_allow_html=True)

    with st.expander("Example gene list & input format", expanded=False):
        st.markdown("""
<div class="guide-section">
<h4>Input format</h4>
<p>Your gene list should be a plain text file with <strong>one gene identifier
per line</strong>, no headers, no commas. Like this:</p>
<div class="gene-example-box">
BRAF<br>
MAP2K1<br>
MAP2K2<br>
MAPK1<br>
MAPK3<br>
KRAS<br>
NRAS<br>
RAF1<br>
MYC<br>
FOS
</div>
<p>This example contains 10 genes from the <strong>MAPK/ERK signaling
pathway</strong>. The built-in sample data uses 20 genes from this
pathway.</p>

<h4>Where should gene lists come from?</h4>
<ul>
<li><strong>Differential expression analysis</strong> &mdash; Genes that are
significantly up- or down-regulated in your RNA-seq / microarray experiment
(e.g. DESeq2, edgeR, limma output).</li>
<li><strong>Genetic screens</strong> &mdash; Hit genes from CRISPR screens,
siRNA screens, or mutagenesis studies.</li>
<li><strong>GWAS / variant studies</strong> &mdash; Genes harboring significant
variants from genome-wide association studies.</li>
<li><strong>Curated pathways</strong> &mdash; Genes from a known pathway of
interest (KEGG, Reactome, WikiPathways).</li>
<li><strong>Co-expression modules</strong> &mdash; Gene clusters from WGCNA or
similar co-expression network analyses.</li>
</ul>

<h4>Tips for good input</h4>
<ul>
<li><strong>10&ndash;200 genes</strong> is the sweet spot. Fewer than 5 may not
produce meaningful network context. More than 500 will be very slow.</li>
<li>Make sure identifiers match the <strong>Gene ID type</strong> you select
(e.g. <code>BRAF</code> for HGNC symbols, <code>ENSG00000157764</code> for
Ensembl).</li>
<li>Remove duplicates and blank lines beforehand.</li>
<li>Gene names are <strong>case-sensitive</strong> &mdash; use official symbols
(e.g. <code>BRAF</code>, not <code>braf</code> or <code>B-Raf</code>).</li>
</ul>
</div>
        """, unsafe_allow_html=True)

    if st.session_state.run_log:
        with st.expander("Run log"):
            st.markdown(st.session_state.run_log)

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
