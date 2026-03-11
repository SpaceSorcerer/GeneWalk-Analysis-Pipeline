"""GeneWalk App -- Web.

Web application for running GeneWalk (if installed) and interactively exploring
gene-function associations. Users can run a new analysis or upload a
``genewalk_results.csv`` from a local GeneWalk run (or the desktop app).

Launch with:  streamlit run app.py
"""

import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

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
SAMPLE_RESULTS_PATH = Path(__file__).parent / "sample_data" / "sample_genewalk_results.csv"

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
# Sidebar -- mode toggle + inputs
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## GeneWalk App")

    app_mode = st.radio(
        "Mode",
        ["Run Analysis", "View Results"],
        help="**Run Analysis** lets you execute GeneWalk with a gene list. "
             "**View Results** lets you upload or demo existing results.",
    )

    st.markdown("---")

    if app_mode == "Run Analysis":
        st.caption("Upload genes, configure, and run GeneWalk.")

        st.markdown('<p class="section-header">Gene List Input</p>',
                    unsafe_allow_html=True)

        input_method = st.radio(
            "Choose input method",
            ["Upload file", "Paste genes", "Use sample genes"],
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
                 "for workstations.",
        )
        nreps_graph = col2.number_input(
            "Graph reps", min_value=1, max_value=20, value=3,
            help="DeepWalk repetitions on the gene-GO network. "
                 "Use 3 for quick testing, 10-15 for publication.",
        )
        nreps_null = st.number_input(
            "Null reps", min_value=1, max_value=20, value=3,
            help="DeepWalk repetitions on randomized networks. "
                 "Should generally match Graph reps.",
        )
        alpha_fdr = st.slider(
            "FDR threshold", 0.01, 1.0, 1.0, 0.01,
            help="Set to 1.0 to keep all results and filter in the dashboard.",
        )
        project_name = st.text_input(
            "Project name", value="genewalk_analysis",
            help="Names the output directory.",
        )

        output_folder = st.text_input(
            "Output folder",
            value="",
            help="Where to save GeneWalk results. Leave blank for system temp.",
            placeholder="Leave blank for default (system temp)",
        )

        st.markdown("---")
        genewalk_installed = is_genewalk_available()
        if not genewalk_installed:
            st.warning(
                "GeneWalk is not installed (or not found on PATH). "
                "Install it with `pip install genewalk`. "
                "Switch to **View Results** mode to upload existing results.",
                icon="\u26a0\ufe0f",
            )
        run_clicked = st.button(
            "Run GeneWalk", type="primary",
            disabled=len(genes) == 0 or not genewalk_installed,
            use_container_width=True,
        )

    else:
        # View Results mode
        st.caption("Upload GeneWalk results and explore interactively.")

        st.markdown('<p class="section-header">Load Results</p>',
                    unsafe_allow_html=True)

        view_input = st.radio(
            "Choose input method",
            ["Upload results CSV", "Use sample data (MAPK/ERK demo)"],
            label_visibility="collapsed",
        )

        if view_input == "Upload results CSV":
            uploaded_results = st.file_uploader(
                "Upload genewalk_results.csv",
                type=["csv"],
                help="Upload a CSV from a previous GeneWalk run.",
            )
            if uploaded_results is not None:
                try:
                    st.session_state.results_df = pd.read_csv(uploaded_results)
                    st.session_state.run_log = "Results loaded from uploaded CSV."
                except Exception as exc:
                    st.error(f"Could not read uploaded CSV: {exc}")
        else:
            if SAMPLE_RESULTS_PATH.exists() and st.session_state.results_df is None:
                st.session_state.results_df = pd.read_csv(SAMPLE_RESULTS_PATH)
                st.session_state.run_log = "Demo results loaded from sample data."

        st.markdown("---")
        st.markdown(
            '<div class="info-tip">'
            "<strong>Need to run GeneWalk?</strong><br>"
            "Switch to <strong>Run Analysis</strong> mode (above) if you have "
            "GeneWalk installed, or download the "
            '<a href="https://github.com/SpaceSorcerer/'
            'GeneWalk-Analysis-Pipeline/releases">Desktop app</a>.'
            "</div>",
            unsafe_allow_html=True,
        )

        # Set defaults so Run Analysis variables exist when not in that mode
        run_clicked = False
        genes = []

# ---------------------------------------------------------------------------
# Handle Run Analysis execution
# ---------------------------------------------------------------------------
if app_mode == "Run Analysis" and run_clicked and genes:
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
                    csv_path = None
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
    <div class="hero-badge">GeneWalk App</div>
    <h1>GeneWalk App</h1>
    <p>Run GeneWalk analyses and explore gene-function associations with interactive
    volcano plots, gene-GO networks, heatmaps, and more.
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
            <p>Switch to <strong>Run Analysis</strong> mode and upload a gene
            list, or switch to <strong>View Results</strong> to upload existing
            GeneWalk output.</p>
        </div>
        """, unsafe_allow_html=True)
    with col_b:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">2</div>
            <h3>Run or Upload</h3>
            <p>Run GeneWalk directly from this app (requires
            <code>pip install genewalk</code>) or upload a previous
            <code>genewalk_results.csv</code>.</p>
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
        "</strong> in the sidebar (View Results mode) to instantly explore the "
        "app with pre-computed MAPK/ERK pathway data.</div>",
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
<li><strong>Run Analysis mode</strong> &mdash; Load a gene list (upload, paste,
or use sample genes), configure parameters, and click <strong>Run GeneWalk</strong>.
Requires GeneWalk to be installed (<code>pip install genewalk</code>).</li>
<li><strong>View Results mode</strong> &mdash; Upload a
<code>genewalk_results.csv</code> from a previous run, or try the built-in
MAPK/ERK demo data.</li>
<li><strong>Explore</strong> &mdash; The dashboard appears with five
interactive tabs: Overview, Per-Gene Explorer, Network, Heatmap, Data Table.</li>
</ol>
</div>
        """, unsafe_allow_html=True)

    with st.expander("What is GeneWalk?", expanded=False):
        st.markdown("""
<div class="guide-section">
<p><a href="https://github.com/churchmanlab/genewalk">GeneWalk</a> identifies
relevant gene functions for a biological context using network representation
learning. Given a list of genes, GeneWalk:</p>
<ol>
<li>Builds a gene&ndash;GO term network from known annotations and
Pathway Commons interactions.</li>
<li>Uses DeepWalk random walks to learn vector representations of genes and
GO terms in the network.</li>
<li>Computes similarity scores and statistical significance for each
gene&ndash;GO term pair.</li>
</ol>
<p><strong>Citation:</strong> Ietswaart et al., <em>GeneWalk identifies relevant
gene functions for a biological context using network representation learning</em>,
Genome Biology 2021.</p>
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
