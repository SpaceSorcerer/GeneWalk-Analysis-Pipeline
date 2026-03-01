"""GeneWalk + GSEA Comparison Pipeline.

Upload a differential expression table (gene, log2FC, padj) and this app
will:
  1. Split genes into up- and down-regulated lists
  2. Run GeneWalk on each list (context-specific functional annotation)
  3. Run GSEA prerank on the full ranked list (pathway enrichment)
  4. Run ORA on each list (over-representation analysis)
  5. Present a unified comparison dashboard

Launch with:  streamlit run comparison_app.py
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from genewalk_app.comparison import (
    cross_method_concordance,
    detect_deg_columns,
    make_ranked_list,
    parse_deg_table,
    shared_go_terms,
    split_deg_lists,
)
from genewalk_app.comparison_viz import (
    concordance_bar,
    deg_overview_volcano,
    direction_volcano,
    gene_set_summary_metrics,
    gsea_dot_plot,
    nes_bar_chart,
    ora_bar_chart,
    shared_terms_bar,
)
from genewalk_app.dashboard import render_dashboard
from genewalk_app.gsea_runner import DEFAULT_GENE_SETS, run_gsea_prerank, run_ora
from genewalk_app.runner import (
    _sanitize_project_name,
    find_results_csv,
    is_genewalk_available,
    load_results,
    run_genewalk,
    save_gene_list,
)
from genewalk_app.styles import get_custom_css

SAMPLE_DEG_PATH = Path(__file__).parent / "sample_data" / "sample_deg_table.csv"

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="GeneWalk + GSEA Comparison",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(get_custom_css(), unsafe_allow_html=True)

# ---------------------------------------------------------------------------
# Session state
# ---------------------------------------------------------------------------
_defaults = {
    "deg_table": None,
    "up_genes": [],
    "down_genes": [],
    "gw_results_up": None,
    "gw_results_down": None,
    "gsea_results": None,
    "ora_results_up": None,
    "ora_results_down": None,
    "comp_run_log": None,
    "analysis_complete": False,
}
for key, val in _defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## Comparison Pipeline")
    st.caption("Upload a DEG table to run GeneWalk + GSEA side by side.")

    # ---- Input ----
    st.markdown('<p class="section-header">Differential Expression Table</p>',
                unsafe_allow_html=True)

    input_method = st.radio(
        "Choose input method",
        ["Upload DEG table", "Use sample data"],
        label_visibility="collapsed",
    )

    deg_df: pd.DataFrame | None = None

    if input_method == "Upload DEG table":
        uploaded = st.file_uploader(
            "Upload DEG CSV/TSV",
            type=["csv", "tsv", "txt"],
            help="CSV or TSV with columns for gene name, log2 fold-change, "
                 "and adjusted p-value (from DESeq2, edgeR, limma, etc.).",
        )
        if uploaded:
            try:
                content = uploaded.getvalue().decode("utf-8")
            except UnicodeDecodeError:
                content = uploaded.getvalue().decode("latin-1")
            sep = "\t" if uploaded.name.endswith((".tsv", ".txt")) else ","
            try:
                deg_df = pd.read_csv(io.StringIO(content), sep=sep)
            except Exception as exc:
                st.error(f"Could not parse file: {exc}")
    else:
        if SAMPLE_DEG_PATH.exists():
            deg_df = pd.read_csv(SAMPLE_DEG_PATH)
            st.info("Sample MAPK/ERK DEG data loaded.")
        else:
            st.warning("Sample DEG data not found.")

    # ---- Column mapping ----
    gene_col = log2fc_col = padj_col = None
    if deg_df is not None:
        detected = detect_deg_columns(deg_df)
        all_cols = deg_df.columns.tolist()

        st.markdown('<p class="section-header">Column Mapping</p>',
                    unsafe_allow_html=True)

        gene_col = st.selectbox(
            "Gene column",
            all_cols,
            index=all_cols.index(detected["gene"]) if detected["gene"] else 0,
            help="Column containing gene identifiers.",
        )
        log2fc_col = st.selectbox(
            "log2FC column",
            all_cols,
            index=all_cols.index(detected["log2fc"]) if detected["log2fc"] else min(1, len(all_cols) - 1),
            help="Column containing log2 fold-change values.",
        )
        padj_col = st.selectbox(
            "Adjusted p-value column",
            all_cols,
            index=all_cols.index(detected["padj"]) if detected["padj"] else min(2, len(all_cols) - 1),
            help="Column containing adjusted p-values (FDR).",
        )

        with st.expander("Preview DEG table", expanded=False):
            st.dataframe(deg_df.head(20), height=250)

    # ---- Thresholds ----
    st.markdown("---")
    st.markdown('<p class="section-header">Thresholds</p>',
                unsafe_allow_html=True)

    fc_threshold = st.slider(
        "|log2FC| threshold", 0.0, 5.0, 1.0, 0.25,
        help="Genes with absolute log2FC above this are considered "
             "up- or down-regulated for GeneWalk and ORA.",
    )
    sig_threshold = st.slider(
        "Significance threshold (padj)", 0.001, 0.5, 0.05, 0.001,
        format="%.3f",
        help="Adjusted p-value cutoff for selecting differentially "
             "expressed genes.",
    )

    # ---- Analysis parameters ----
    st.markdown("---")
    st.markdown('<p class="section-header">Analysis Parameters</p>',
                unsafe_allow_html=True)

    id_type = st.selectbox(
        "Gene ID type",
        ["hgnc_symbol", "hgnc_id", "ensembl_id", "mgi_id", "rgd_id",
         "entrez_human", "entrez_mouse"],
        index=0,
    )

    col_a, col_b = st.columns(2)
    nproc = col_a.number_input("CPU cores", 1, 32, 4)
    nreps = col_b.number_input("DeepWalk reps", 1, 20, 3,
                                help="Used for both graph and null reps.")

    project_name = st.text_input("Project name", value="comparison_analysis")

    # ---- Gene set libraries ----
    st.markdown("---")
    st.markdown('<p class="section-header">Gene Set Libraries</p>',
                unsafe_allow_html=True)

    selected_gene_sets = st.multiselect(
        "Libraries for GSEA & ORA",
        DEFAULT_GENE_SETS,
        default=DEFAULT_GENE_SETS[:3],
        help="Gene set databases to test against. GO libraries match "
             "GeneWalk's ontology. KEGG/Reactome/Hallmarks add pathway "
             "context that GeneWalk doesn't cover.",
    )

    # ---- What to run ----
    st.markdown("---")
    st.markdown('<p class="section-header">Analyses to Run</p>',
                unsafe_allow_html=True)

    run_gw = st.checkbox("GeneWalk (up & down)", value=True,
                          help="Run GeneWalk separately on up- and down-regulated genes.")
    run_gsea = st.checkbox("GSEA prerank", value=True,
                            help="Run GSEA on the full ranked gene list (by log2FC).")
    run_ora_check = st.checkbox("ORA (up & down)", value=True,
                                 help="Run over-representation analysis on each gene list.")

    genewalk_installed = is_genewalk_available()
    if run_gw and not genewalk_installed:
        st.warning("GeneWalk is not installed. GeneWalk analysis will be skipped.")
        run_gw = False

    st.markdown("---")
    run_clicked = st.button(
        "Run Comparison Analysis",
        type="primary",
        disabled=deg_df is None or not selected_gene_sets,
        use_container_width=True,
    )

    # ---- Upload previous results ----
    st.markdown("---")
    st.markdown('<p class="section-header">Or Upload Previous Results</p>',
                unsafe_allow_html=True)
    st.caption(
        "Upload individual result CSVs from previous runs to skip analysis."
    )
    up_gw_upload = st.file_uploader("GeneWalk Up results CSV", type=["csv"], key="up_gw_csv")
    down_gw_upload = st.file_uploader("GeneWalk Down results CSV", type=["csv"], key="down_gw_csv")
    gsea_upload = st.file_uploader("GSEA results CSV", type=["csv"], key="gsea_csv")


# ---------------------------------------------------------------------------
# Handle uploads of previous results
# ---------------------------------------------------------------------------
for upload_key, state_key in [
    (up_gw_upload, "gw_results_up"),
    (down_gw_upload, "gw_results_down"),
    (gsea_upload, "gsea_results"),
]:
    if upload_key is not None:
        try:
            st.session_state[state_key] = pd.read_csv(upload_key)
        except Exception as exc:
            st.error(f"Could not read CSV: {exc}")

# ---------------------------------------------------------------------------
# Run analysis
# ---------------------------------------------------------------------------
if run_clicked and deg_df is not None and gene_col and log2fc_col and padj_col:
    try:
        parsed = parse_deg_table(deg_df, gene_col, log2fc_col, padj_col)
    except ValueError as exc:
        st.error(str(exc))
        st.stop()

    st.session_state.deg_table = parsed
    up_genes, down_genes = split_deg_lists(parsed, fc_threshold, sig_threshold)
    st.session_state.up_genes = up_genes
    st.session_state.down_genes = down_genes

    if not up_genes and not down_genes:
        st.error(
            f"No genes pass the thresholds (|log2FC| >= {fc_threshold}, "
            f"padj <= {sig_threshold}). Try relaxing the cutoffs."
        )
        st.stop()

    log_parts: list[str] = [
        f"**Up-regulated genes:** {len(up_genes)}",
        f"**Down-regulated genes:** {len(down_genes)}",
    ]

    with st.status("Running comparison analysis...", expanded=True) as status:
        progress = st.empty()
        safe_project = _sanitize_project_name(project_name)
        tmp = Path(tempfile.mkdtemp())

        # --- GeneWalk ---
        if run_gw:
            for label, genes, state_key in [
                ("up-regulated", up_genes, "gw_results_up"),
                ("down-regulated", down_genes, "gw_results_down"),
            ]:
                if not genes:
                    log_parts.append(f"Skipped GeneWalk ({label}): no genes.")
                    continue
                progress.write(f"Running GeneWalk on {len(genes)} {label} genes...")
                gene_file = save_gene_list(genes, tmp / f"genes_{label}.txt")
                gw_project = f"{safe_project}_{label.replace('-', '_')}"

                def _gw_progress(msg: str, _label=label) -> None:
                    progress.write(f"GeneWalk ({_label}): {msg}")

                result = run_genewalk(
                    gene_file=gene_file,
                    project=gw_project,
                    id_type=id_type,
                    nproc=nproc,
                    nreps_graph=nreps,
                    nreps_null=nreps,
                    alpha_fdr=1.0,
                    on_progress=_gw_progress,
                )
                csv_path = find_results_csv(result["output_dir"])
                if csv_path:
                    try:
                        st.session_state[state_key] = load_results(csv_path)
                        log_parts.append(f"GeneWalk ({label}): complete.")
                    except ValueError as exc:
                        log_parts.append(f"GeneWalk ({label}): CSV error - {exc}")
                else:
                    log_parts.append(
                        f"GeneWalk ({label}): finished (rc={result['return_code']}) "
                        f"but no results CSV found."
                    )

        # --- GSEA prerank ---
        if run_gsea:
            progress.write("Running GSEA prerank on full ranked list...")
            ranked = make_ranked_list(parsed)
            if len(ranked) < 15:
                log_parts.append("Skipped GSEA: fewer than 15 ranked genes.")
            else:
                def _gsea_progress(msg: str) -> None:
                    progress.write(msg)

                gsea_res = run_gsea_prerank(
                    ranked_genes=ranked,
                    gene_sets=selected_gene_sets,
                    permutation_num=1000,
                    min_size=15,
                    max_size=500,
                    threads=nproc,
                    on_progress=_gsea_progress,
                )
                st.session_state.gsea_results = gsea_res
                n_sig = len(gsea_res[gsea_res["fdr"] <= 0.25]) if "fdr" in gsea_res.columns else 0
                log_parts.append(f"GSEA prerank: {len(gsea_res)} terms tested, {n_sig} significant (FDR <= 0.25).")

        # --- ORA ---
        if run_ora_check:
            for label, genes, state_key in [
                ("up-regulated", up_genes, "ora_results_up"),
                ("down-regulated", down_genes, "ora_results_down"),
            ]:
                if not genes:
                    log_parts.append(f"Skipped ORA ({label}): no genes.")
                    continue
                progress.write(f"Running ORA on {len(genes)} {label} genes...")

                def _ora_progress(msg: str, _label=label) -> None:
                    progress.write(f"ORA ({_label}): {msg}")

                ora_res = run_ora(
                    gene_list=genes,
                    gene_sets=selected_gene_sets,
                    on_progress=_ora_progress,
                )
                st.session_state[state_key] = ora_res
                n_sig = len(ora_res[ora_res["fdr"] <= 0.05]) if "fdr" in ora_res.columns else 0
                log_parts.append(f"ORA ({label}): {len(ora_res)} terms, {n_sig} significant (FDR <= 0.05).")

        st.session_state.comp_run_log = "\n\n".join(log_parts)
        st.session_state.analysis_complete = True
        status.update(label="Analysis complete!", state="complete")


# ---------------------------------------------------------------------------
# Check if we have any results to show
# ---------------------------------------------------------------------------
has_results = any([
    st.session_state.gw_results_up is not None,
    st.session_state.gw_results_down is not None,
    st.session_state.gsea_results is not None,
    st.session_state.ora_results_up is not None,
    st.session_state.ora_results_down is not None,
])

# ---------------------------------------------------------------------------
# Hero banner
# ---------------------------------------------------------------------------
st.markdown("""
<div class="hero-banner">
    <div class="hero-badge">Comparison Pipeline</div>
    <h1>GeneWalk + GSEA Analysis</h1>
    <p>Upload a differential expression table to run context-specific
    functional analysis (GeneWalk) alongside pathway enrichment (GSEA/ORA)
    &mdash; all in one place.</p>
</div>
""", unsafe_allow_html=True)

if not has_results:
    # ---- Landing page ----
    col_a, col_b, col_c = st.columns(3, gap="medium")
    with col_a:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">1</div>
            <h3>Upload DEG Table</h3>
            <p>Upload a CSV/TSV from DESeq2, edgeR, or limma with gene names,
            log2 fold-changes, and adjusted p-values.</p>
        </div>
        """, unsafe_allow_html=True)
    with col_b:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">2</div>
            <h3>Configure &amp; Run</h3>
            <p>Map your columns, set thresholds, choose gene set libraries,
            and run GeneWalk + GSEA + ORA in a single click.</p>
        </div>
        """, unsafe_allow_html=True)
    with col_c:
        st.markdown("""
        <div class="step-card">
            <div class="step-number">3</div>
            <h3>Compare Results</h3>
            <p>Explore up vs down GeneWalk results, GSEA enrichment,
            ORA hits, and cross-method concordance.</p>
        </div>
        """, unsafe_allow_html=True)

    with st.expander("What this app does", expanded=False):
        st.markdown("""
<div class="guide-section">
<h4>Three complementary analyses in one pipeline</h4>
<ul>
<li><strong>GeneWalk</strong> &mdash; Runs separately on your up- and
down-regulated gene lists. Identifies which GO functions are specifically
relevant to each set of genes using network representation learning.
Answers: <em>"What functions are contextually relevant to these genes?"</em></li>
<li><strong>GSEA prerank</strong> &mdash; Uses the full ranked gene list
(sorted by log2FC) to find pathways enriched at either extreme. Unlike
GeneWalk, GSEA incorporates the <strong>magnitude</strong> of change.
Answers: <em>"Which pathways have coordinated up- or down-regulation?"</em></li>
<li><strong>ORA</strong> &mdash; Over-representation analysis tests whether
your up/down gene lists have more genes from a pathway than expected by
chance. A simpler, threshold-based complement to GSEA.
Answers: <em>"Are any pathways statistically over-represented in my gene list?"</em></li>
</ul>

<h4>Why run all three?</h4>
<p>Each method has blind spots. GeneWalk finds context-specific functions
but ignores fold-change magnitude. GSEA uses magnitude but misses
gene-specific functional nuance. ORA is simple and intuitive but loses
information by binarizing genes as significant/not. Running all three and
looking for <strong>concordance</strong> across methods gives the most
robust biological interpretation.</p>
</div>
        """, unsafe_allow_html=True)

    st.markdown(
        '<div class="app-footer">Built with '
        '<a href="https://streamlit.io">Streamlit</a> &middot; '
        'Analysis by <a href="https://github.com/churchmanlab/genewalk">'
        "GeneWalk</a> + "
        '<a href="https://github.com/zqfang/GSEApy">GSEApy</a></div>',
        unsafe_allow_html=True,
    )
    st.stop()

# =========================================================================
# Results dashboard
# =========================================================================
gw_up = st.session_state.gw_results_up
gw_down = st.session_state.gw_results_down
gsea_res = st.session_state.gsea_results
ora_up = st.session_state.ora_results_up
ora_down = st.session_state.ora_results_down
up_genes = st.session_state.up_genes
down_genes = st.session_state.down_genes
deg = st.session_state.deg_table

# ---- Run log ----
if st.session_state.comp_run_log:
    with st.expander("Run log", expanded=False):
        st.markdown(st.session_state.comp_run_log)

# ---- Filters ----
st.markdown('<p class="section-header">Dashboard Filters</p>', unsafe_allow_html=True)
fc1, fc2, fc3 = st.columns(3)

gw_padj_col = fc1.selectbox(
    "GeneWalk p-value column",
    ["global_padj", "gene_padj"],
    help="Which GeneWalk adjusted p-value to use for filtering.",
)
gw_padj_threshold = fc2.slider(
    "GeneWalk significance", 0.001, 1.0, 0.05, 0.001, format="%.3f",
)
gsea_fdr_threshold = fc3.slider(
    "GSEA FDR threshold", 0.01, 1.0, 0.25, 0.01, format="%.2f",
)

# ---- Summary metrics ----
metrics = gene_set_summary_metrics(
    up_genes, down_genes, gw_up, gw_down, gsea_res, ora_up, ora_down,
    padj_col=gw_padj_col, padj_threshold=gw_padj_threshold,
    gsea_fdr=gsea_fdr_threshold,
)
st.markdown("")
mc = st.columns(6)
mc[0].metric("Up Genes", metrics["n_up"])
mc[1].metric("Down Genes", metrics["n_down"])
mc[2].metric("GW Up Sig Terms", metrics["gw_up_sig_terms"])
mc[3].metric("GW Down Sig Terms", metrics["gw_down_sig_terms"])
mc[4].metric("GSEA Sig Up", metrics["gsea_sig_up"])
mc[5].metric("GSEA Sig Down", metrics["gsea_sig_down"])

# ---- Interpretive guide ----
with st.expander("How to interpret these results", expanded=False):
    st.markdown("""
<div class="guide-section">
<h4>Reading the comparison dashboard</h4>
<p>This dashboard shows results from three complementary methods. Here is
what to look for:</p>

<h4>GeneWalk results (Up vs Down tabs)</h4>
<p>GeneWalk identifies which GO functions are <strong>contextually
relevant</strong> to each gene set. A function significant in the
up-regulated set but not the down-regulated set suggests it is specifically
tied to the biology of activation/increase. Functions significant in
<strong>both</strong> directions may represent shared pathway architecture.</p>

<h4>GSEA results</h4>
<p><strong>NES > 0</strong> means genes in that pathway tend to be
up-regulated. <strong>NES < 0</strong> means down-regulated. Unlike
GeneWalk, GSEA incorporates the <em>magnitude</em> of fold-change, so
highly-changed genes contribute more to the enrichment score. FDR <= 0.25
is the conventional GSEA significance threshold.</p>

<h4>Cross-method concordance</h4>
<p>The strongest findings are those identified by <strong>multiple
methods</strong>. If a GO biological process is significant in GeneWalk
AND its parent pathway is enriched in GSEA AND over-represented in ORA,
that is high-confidence evidence of biological relevance.</p>

<h4>Key caveats</h4>
<ul>
<li>GeneWalk and ORA use only gene <em>membership</em> (in/out of your
list). GSEA uses the full ranked list with fold-change magnitudes.</li>
<li>GeneWalk uses GO terms directly. GSEA/ORA use curated pathway
databases (KEGG, Reactome, Hallmarks) which may group terms differently.</li>
<li>Term names differ across databases, so exact string matching between
methods will miss some overlaps.</li>
</ul>
</div>
    """, unsafe_allow_html=True)

# ---- Tabs ----
st.markdown("")
available_tabs = ["Overview"]
if gw_up is not None or gw_down is not None:
    available_tabs.append("GeneWalk Comparison")
if gsea_res is not None and not gsea_res.empty:
    available_tabs.append("GSEA Results")
if ora_up is not None or ora_down is not None:
    available_tabs.append("ORA Results")
available_tabs.append("Cross-Method")
available_tabs.append("Data Tables")

tabs = st.tabs(available_tabs)
tab_idx = 0

# ---- Tab: Overview ----
with tabs[tab_idx]:
    tab_idx += 1
    if deg is not None and not deg.empty:
        st.plotly_chart(
            deg_overview_volcano(deg, fc_threshold, sig_threshold),
            width="stretch",
        )

    if gw_up is not None and gw_down is not None:
        st.plotly_chart(
            direction_volcano(gw_up, gw_down, padj_col=gw_padj_col,
                              padj_threshold=gw_padj_threshold),
            width="stretch",
        )

# ---- Tab: GeneWalk Comparison ----
if "GeneWalk Comparison" in available_tabs:
    with tabs[tab_idx]:
        tab_idx += 1

        gw_sub = st.tabs(["Up-Regulated", "Down-Regulated", "Shared GO Terms"])

        with gw_sub[0]:
            if gw_up is not None and not gw_up.empty:
                st.markdown(f"**{len(up_genes)} up-regulated genes** analyzed by GeneWalk")
                render_dashboard(gw_up)
            else:
                st.info("No GeneWalk results for up-regulated genes.")

        with gw_sub[1]:
            if gw_down is not None and not gw_down.empty:
                st.markdown(f"**{len(down_genes)} down-regulated genes** analyzed by GeneWalk")
                render_dashboard(gw_down)
            else:
                st.info("No GeneWalk results for down-regulated genes.")

        with gw_sub[2]:
            if gw_up is not None and gw_down is not None:
                shared = shared_go_terms(
                    gw_up, gw_down,
                    padj_col=gw_padj_col,
                    padj_threshold=gw_padj_threshold,
                )
                if not shared.empty:
                    st.plotly_chart(shared_terms_bar(shared), width="stretch")
                    st.dataframe(shared, width="stretch", height=400)
                else:
                    st.info("No GO terms are significant in both up- and down-regulated gene sets at the current threshold.")
            else:
                st.info("Need both up and down GeneWalk results to compare.")

# ---- Tab: GSEA Results ----
if "GSEA Results" in available_tabs:
    with tabs[tab_idx]:
        tab_idx += 1

        gsea_sub = st.tabs(["NES Bar Chart", "Dot Plot", "Full Table"])

        with gsea_sub[0]:
            top_n_nes = st.slider("Top pathways", 10, 60, 30, key="nes_top")
            st.plotly_chart(
                nes_bar_chart(gsea_res, fdr_threshold=gsea_fdr_threshold,
                              top_n=top_n_nes),
                width="stretch",
            )

        with gsea_sub[1]:
            st.plotly_chart(
                gsea_dot_plot(gsea_res, fdr_threshold=gsea_fdr_threshold),
                width="stretch",
            )

        with gsea_sub[2]:
            if gsea_res is not None and not gsea_res.empty:
                display_gsea = gsea_res.copy()
                if "fdr" in display_gsea.columns:
                    display_gsea = display_gsea[display_gsea["fdr"] <= gsea_fdr_threshold]
                st.markdown(f"**{len(display_gsea)}** significant pathways (FDR <= {gsea_fdr_threshold})")
                st.dataframe(display_gsea, width="stretch", height=500)

# ---- Tab: ORA Results ----
if "ORA Results" in available_tabs:
    with tabs[tab_idx]:
        tab_idx += 1

        ora_sub = st.tabs(["Up-Regulated", "Down-Regulated"])

        with ora_sub[0]:
            if ora_up is not None and not ora_up.empty:
                st.plotly_chart(
                    ora_bar_chart(ora_up, label="Up-regulated"),
                    width="stretch",
                )
                with st.expander("Full ORA results (up)", expanded=False):
                    st.dataframe(ora_up, width="stretch", height=400)
            else:
                st.info("No ORA results for up-regulated genes.")

        with ora_sub[1]:
            if ora_down is not None and not ora_down.empty:
                st.plotly_chart(
                    ora_bar_chart(ora_down, label="Down-regulated"),
                    width="stretch",
                )
                with st.expander("Full ORA results (down)", expanded=False):
                    st.dataframe(ora_down, width="stretch", height=400)
            else:
                st.info("No ORA results for down-regulated genes.")

# ---- Tab: Cross-Method ----
with tabs[tab_idx]:
    tab_idx += 1

    st.markdown("""
    <div class="info-tip">Terms found by multiple methods (GeneWalk, GSEA, ORA)
    represent the highest-confidence functional findings. Terms matched here
    use exact name matching &mdash; some overlaps may be missed due to naming
    differences across databases.</div>
    """, unsafe_allow_html=True)

    # Merge GW results for concordance
    gw_combined = pd.DataFrame()
    if gw_up is not None and not gw_up.empty:
        gw_combined = gw_up
    if gw_down is not None and not gw_down.empty:
        gw_combined = pd.concat([gw_combined, gw_down], ignore_index=True) if not gw_combined.empty else gw_down

    ora_combined = pd.DataFrame()
    if ora_up is not None and not ora_up.empty:
        ora_combined = ora_up
    if ora_down is not None and not ora_down.empty:
        ora_combined = pd.concat([ora_combined, ora_down], ignore_index=True) if not ora_combined.empty else ora_down

    concordance = cross_method_concordance(
        gw_combined,
        gsea_res if gsea_res is not None else pd.DataFrame(),
        ora_combined,
        gw_padj_col=gw_padj_col,
        gw_padj_threshold=gw_padj_threshold,
        gsea_fdr_threshold=gsea_fdr_threshold,
    )

    if not concordance.empty:
        st.plotly_chart(concordance_bar(concordance), width="stretch")
        st.dataframe(concordance, width="stretch", height=400)
    else:
        st.info(
            "No terms found by multiple methods at the current thresholds. "
            "Try relaxing the FDR/p-value cutoffs, or note that term names "
            "differ across databases (GO vs KEGG vs Reactome)."
        )

# ---- Tab: Data Tables ----
with tabs[tab_idx]:
    dt_sub = st.tabs(["DEG Table", "GW Up", "GW Down", "GSEA", "ORA Up", "ORA Down"])

    datasets = [
        (deg, "deg_table.csv"),
        (gw_up, "genewalk_up_results.csv"),
        (gw_down, "genewalk_down_results.csv"),
        (gsea_res, "gsea_prerank_results.csv"),
        (ora_up, "ora_up_results.csv"),
        (ora_down, "ora_down_results.csv"),
    ]

    for sub_tab, (data, filename) in zip(dt_sub, datasets):
        with sub_tab:
            if data is not None and not data.empty:
                st.markdown(f"**{len(data):,}** rows")
                st.dataframe(data, width="stretch", height=500)
                csv_buf = io.StringIO()
                data.to_csv(csv_buf, index=False)
                st.download_button(
                    f"Download {filename}",
                    csv_buf.getvalue(),
                    file_name=filename,
                    mime="text/csv",
                    key=f"dl_{filename}",
                )
            else:
                st.info("No data available.")

# Footer
st.markdown(
    '<div class="app-footer">Built with '
    '<a href="https://streamlit.io">Streamlit</a> &middot; '
    'Analysis by <a href="https://github.com/churchmanlab/genewalk">'
    "GeneWalk</a> + "
    '<a href="https://github.com/zqfang/GSEApy">GSEApy</a> &middot; '
    "Ietswaart et al., <em>Genome Biology</em> 2021</div>",
    unsafe_allow_html=True,
)
