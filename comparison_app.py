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
    gsea_leading_edge_table,
    nes_bar_chart,
    ora_bar_chart,
    ora_dot_plot,
    shared_terms_bar,
)
from genewalk_app.deg_visualizations import (
    basemean_distribution,
    deg_category_pie,
    deg_summary_counts,
    lfc_se_scatter,
    lfc_vs_stat_scatter,
    ma_plot,
    padj_histogram,
    pvalue_histogram,
    top_genes_bar,
)
from genewalk_app.dashboard import render_dashboard
from genewalk_app.gene_investigator import render_gene_investigator
from genewalk_app.gsea_runner import DEFAULT_GENE_SETS, run_gsea_prerank, run_ora
from genewalk_app.runner import (
    _sanitize_project_name,
    find_results_csv,
    is_genewalk_available,
    load_results,
    run_genewalk,
    save_gene_list,
)
from genewalk_app.splicing import (
    filter_splicing_events,
    parse_rmats_file,
    parse_vasttools,
    splicing_summary,
)
from genewalk_app.styles import get_custom_css

SAMPLE_DEG_PATH = Path(__file__).parent / "sample_data" / "sample_deg_table.csv"


def run_comparison_ui() -> None:
    """Render the full GSEA + GeneWalk comparison UI.

    Can be called from the standalone ``comparison_app.py`` (after
    ``st.set_page_config``) **or** embedded inside ``desktop.py`` when the
    user selects the comparison mode.
    """
    # -------------------------------------------------------------------
    # Session state
    # -------------------------------------------------------------------
    _defaults = {
        "deg_table": None,
        "deg_raw": None,
        "deg_col_map": {},
        "up_genes": [],
        "down_genes": [],
        "gw_results_up": None,
        "gw_results_down": None,
        "gsea_results": None,
        "ora_results_up": None,
        "ora_results_down": None,
        "splicing_data": None,
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
                     "and adjusted p-value (from DESeq2, edgeR, limma, etc.). "
                     "R-exported CSVs (write.csv) are auto-detected.",
            )
            if uploaded:
                try:
                    content = uploaded.getvalue().decode("utf-8")
                except UnicodeDecodeError:
                    content = uploaded.getvalue().decode("latin-1")
                sep = "\t" if uploaded.name.endswith((".tsv", ".txt")) else ","
                try:
                    deg_df = pd.read_csv(io.StringIO(content), sep=sep)
                    # Handle R-exported CSVs: if first column is unnamed
                    # (Unnamed: 0) and contains string gene names, rename it
                    first_col = deg_df.columns[0]
                    if (
                        first_col.startswith("Unnamed")
                        or first_col.strip() == ""
                        or first_col == "X"
                        or first_col == "...1"
                    ):
                        # Check that the column contains string values
                        # (gene names), not numbers
                        sample = deg_df[first_col].dropna().head(20)
                        numeric_count = pd.to_numeric(
                            sample, errors="coerce"
                        ).notna().sum()
                        if numeric_count < len(sample) * 0.5:
                            deg_df = deg_df.rename(columns={first_col: "gene"})
                            st.info(
                                f"Detected R-exported CSV format. "
                                f"Renamed column '{first_col}' to 'gene'."
                            )
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
        basemean_col_name = lfcse_col_name = stat_col_name = pvalue_col_name = None
        if deg_df is not None:
            detected = detect_deg_columns(deg_df)
            all_cols = deg_df.columns.tolist()
            none_cols = ["(none)"] + all_cols

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

            with st.expander("Additional DESeq2 columns (optional)", expanded=False):
                basemean_col_name = st.selectbox(
                    "baseMean column",
                    none_cols,
                    index=none_cols.index(detected["basemean"]) if detected["basemean"] else 0,
                    help="Mean of normalized counts (DESeq2 baseMean).",
                )
                lfcse_col_name = st.selectbox(
                    "lfcSE column",
                    none_cols,
                    index=none_cols.index(detected["lfcse"]) if detected["lfcse"] else 0,
                    help="Standard error of the log2FC estimate.",
                )
                stat_col_name = st.selectbox(
                    "Wald statistic column",
                    none_cols,
                    index=none_cols.index(detected["stat"]) if detected["stat"] else 0,
                    help="Wald test statistic (DESeq2 stat column). "
                         "Used as the GSEA ranking metric when available — "
                         "produces better results than log2FC alone.",
                )
                pvalue_col_name = st.selectbox(
                    "Raw p-value column",
                    none_cols,
                    index=none_cols.index(detected["pvalue"]) if detected["pvalue"] else 0,
                    help="Unadjusted p-value (before BH correction).",
                )
                # Clean up "(none)" selections
                if basemean_col_name == "(none)":
                    basemean_col_name = None
                if lfcse_col_name == "(none)":
                    lfcse_col_name = None
                if stat_col_name == "(none)":
                    stat_col_name = None
                if pvalue_col_name == "(none)":
                    pvalue_col_name = None

            with st.expander("Preview DEG table", expanded=False):
                st.dataframe(deg_df.head(20), height=250)

        # ---- Thresholds ----
        st.markdown("---")
        st.markdown('<p class="section-header">Thresholds</p>',
                    unsafe_allow_html=True)

        fc_threshold = st.number_input(
            "|log2FC| threshold",
            min_value=0.0,
            max_value=20.0,
            value=1.0,
            step=0.1,
            format="%.2f",
            help="Genes with absolute log2FC above this are considered "
                 "up- or down-regulated for GeneWalk and ORA. "
                 "Type any value (e.g. 0.40, 0.58, 1.5).",
        )
        sig_threshold = st.number_input(
            "Significance threshold (padj)",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.001,
            format="%.4f",
            help="Adjusted p-value cutoff for selecting differentially "
                 "expressed genes. Type any value (e.g. 0.01, 0.05, 0.1).",
        )

        # baseMean threshold — only shown when the user mapped a baseMean column
        basemean_threshold = None
        if basemean_col_name and basemean_col_name != "(none)":
            basemean_threshold = st.number_input(
                "baseMean threshold (>=)",
                min_value=0.0,
                value=0.0,
                step=10.0,
                format="%.1f",
                help="Filter out lowly-expressed genes before splitting into "
                     "up/down lists. DESeq2 recommends baseMean >= 10–20 for "
                     "reliable results. Set to 0 to disable.",
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

        output_folder = st.text_input(
            "Output folder",
            value="",
            help="Where to save GeneWalk and GSEA results. Leave blank to "
                 "use the system temp directory. Examples: "
                 "C:\\Users\\you\\genewalk_output or ~/genewalk_output",
            placeholder="Leave blank for default (system temp)",
        )

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

        # ---- Splicing data (optional) ----
        st.markdown("---")
        st.markdown('<p class="section-header">Splicing Data (Optional)</p>',
                    unsafe_allow_html=True)
        st.caption(
            "Upload differential splicing results (vast-tools or rMATS) "
            "to integrate into the Gene Investigator."
        )
        splicing_source = st.radio(
            "Splicing source",
            ["vast-tools", "rMATS"],
            key="splicing_source_radio",
            label_visibility="collapsed",
        )
        splicing_upload = st.file_uploader(
            "Upload splicing results",
            type=["csv", "tsv", "txt", "tab"],
            key="splicing_upload",
        )
        rmats_event = "SE"
        if splicing_source == "rMATS":
            rmats_event = st.selectbox(
                "Event type",
                ["SE", "MXE", "A3SS", "A5SS", "RI"],
                key="rmats_event_comp",
            )


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

    # Handle splicing upload
    if splicing_upload is not None:
        try:
            try:
                spl_content = splicing_upload.getvalue().decode("utf-8")
            except UnicodeDecodeError:
                spl_content = splicing_upload.getvalue().decode("latin-1")
            spl_sep = "\t" if splicing_upload.name.endswith((".tsv", ".txt", ".tab")) else ","
            spl_raw = pd.read_csv(io.StringIO(spl_content), sep=spl_sep)
            if splicing_source == "vast-tools":
                st.session_state.splicing_data = parse_vasttools(spl_raw)
            else:
                st.session_state.splicing_data = parse_rmats_file(spl_raw, event_type=rmats_event)
        except Exception as exc:
            st.error(f"Could not parse splicing file: {exc}")

    # ---------------------------------------------------------------------------
    # Run analysis
    # ---------------------------------------------------------------------------
    if run_clicked and deg_df is not None and gene_col and log2fc_col and padj_col:
        try:
            parsed = parse_deg_table(deg_df, gene_col, log2fc_col, padj_col)
        except ValueError as exc:
            st.error(str(exc))
            st.stop()

        # Carry forward extra DESeq2 columns under standardized names
        col_map = {"gene": gene_col, "log2fc": log2fc_col, "padj": padj_col}
        if basemean_col_name and basemean_col_name in deg_df.columns and basemean_col_name not in (gene_col, log2fc_col, padj_col):
            if basemean_col_name in parsed.columns:
                parsed = parsed.rename(columns={basemean_col_name: "baseMean"})
            col_map["basemean"] = "baseMean"
        if lfcse_col_name and lfcse_col_name in deg_df.columns and lfcse_col_name not in (gene_col, log2fc_col, padj_col):
            if lfcse_col_name in parsed.columns:
                parsed = parsed.rename(columns={lfcse_col_name: "lfcSE"})
            col_map["lfcse"] = "lfcSE"
        if stat_col_name and stat_col_name in deg_df.columns and stat_col_name not in (gene_col, log2fc_col, padj_col):
            if stat_col_name in parsed.columns:
                parsed = parsed.rename(columns={stat_col_name: "stat"})
            col_map["stat"] = "stat"
        if pvalue_col_name and pvalue_col_name in deg_df.columns and pvalue_col_name not in (gene_col, log2fc_col, padj_col):
            if pvalue_col_name in parsed.columns:
                parsed = parsed.rename(columns={pvalue_col_name: "pvalue"})
            col_map["pvalue"] = "pvalue"

        # Convert extra columns to numeric
        for extra_col in ["baseMean", "lfcSE", "stat", "pvalue"]:
            if extra_col in parsed.columns:
                parsed[extra_col] = pd.to_numeric(parsed[extra_col], errors="coerce")

        st.session_state.deg_table = parsed
        st.session_state.deg_raw = deg_df
        st.session_state.deg_col_map = col_map

        # Apply baseMean pre-filter to remove unreliable low-count genes.
        # This filter applies to ALL analyses (GeneWalk, GSEA, ORA) so
        # that all methods analyze the same quality-filtered gene universe.
        filtered_deg = parsed
        if basemean_threshold and basemean_threshold > 0 and "baseMean" in parsed.columns:
            filtered_deg = parsed[parsed["baseMean"] >= basemean_threshold]

        up_genes, down_genes = split_deg_lists(filtered_deg, fc_threshold, sig_threshold)
        st.session_state.up_genes = up_genes
        st.session_state.down_genes = down_genes

        if not up_genes and not down_genes:
            threshold_msg = f"|log2FC| >= {fc_threshold}, padj <= {sig_threshold}"
            if basemean_threshold and basemean_threshold > 0:
                threshold_msg += f", baseMean >= {basemean_threshold}"
            st.error(
                f"No genes pass the thresholds ({threshold_msg}). "
                f"Try relaxing the cutoffs."
            )
            st.stop()

        log_parts: list[str] = [
            f"**Total measured genes:** {len(parsed)}",
        ]
        if len(filtered_deg) < len(parsed):
            log_parts.append(
                f"**After baseMean filter (>= {basemean_threshold:.0f}):** "
                f"{len(filtered_deg)} genes"
            )
        log_parts.extend([
            f"**Up-regulated genes** (|log2FC| >= {fc_threshold}, padj <= {sig_threshold}): {len(up_genes)}",
            f"**Down-regulated genes** (|log2FC| >= {fc_threshold}, padj <= {sig_threshold}): {len(down_genes)}",
        ])

        with st.status("Running comparison analysis...", expanded=True) as status:
            progress = st.empty()
            safe_project = _sanitize_project_name(project_name)

            # Resolve output directory
            base_folder: Path | None = None
            if output_folder.strip():
                base_folder = Path(output_folder).expanduser()
                base_folder.mkdir(parents=True, exist_ok=True)
                tmp = base_folder
            else:
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
                        base_folder=base_folder,
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
                # Use filtered_deg (baseMean-filtered) so GSEA analyzes
                # the same quality-filtered gene set as GeneWalk/ORA.
                # Note: GSEA gets ALL genes that pass the quality filter
                # (not just significant ones), because prerank needs the
                # full ranked list to compute enrichment scores properly.
                ranked = make_ranked_list(filtered_deg)

                # Report which ranking metric was selected
                metric_used = ranked.attrs.get("ranking_metric", "log2fc")
                metric_labels = {
                    "stat": "Wald statistic (DESeq2 stat column)",
                    "signed_pvalue": "-log10(p-value) \u00d7 sign(log2FC)",
                    "log2fc": "log2 fold-change",
                }
                log_parts.append(
                    f"GSEA ranking metric: **{metric_labels.get(metric_used, metric_used)}**"
                )
                if metric_used == "log2fc":
                    log_parts.append(
                        "\u26a0\ufe0f *Tip: Map the DESeq2 **stat** (Wald statistic) column "
                        "in the sidebar for better GSEA resolution and fewer tied ranks.*"
                    )

                # Diagnostic: check ranked list has both positive and negative scores
                if len(ranked) > 0 and "score" in ranked.columns:
                    scores = ranked["score"]
                    n_pos = (scores > 0).sum()
                    n_neg = (scores < 0).sum()
                    if n_pos == 0 or n_neg == 0:
                        log_parts.append(
                            "**Warning:** All ranking scores have the same sign. "
                            "GSEA expects scores with both positive "
                            "(up-regulated) and negative (down-regulated) values. "
                            "Check that the correct column was selected."
                        )

                if len(ranked) < 15:
                    log_parts.append("Skipped GSEA: fewer than 15 ranked genes.")
                else:
                    def _gsea_progress(msg: str) -> None:
                        progress.write(msg)

                    gsea_outdir = tmp / "gsea_prerank" if base_folder else None
                    gsea_res = run_gsea_prerank(
                        ranked_genes=ranked,
                        gene_sets=selected_gene_sets,
                        outdir=gsea_outdir,
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

                    ora_outdir = tmp / f"ora_{label.replace('-', '_')}" if base_folder else None
                    # Use quality-filtered genes as ORA background —
                    # this is the correct universe of reliably measured
                    # genes, not the entire genome.
                    all_measured_genes = filtered_deg["gene"].unique().tolist()
                    ora_res = run_ora(
                        gene_list=genes,
                        gene_sets=selected_gene_sets,
                        background=all_measured_genes,
                        outdir=ora_outdir,
                        on_progress=_ora_progress,
                    )
                    st.session_state[state_key] = ora_res
                    n_sig = len(ora_res[ora_res["fdr"] <= 0.05]) if "fdr" in ora_res.columns else 0
                    log_parts.append(f"ORA ({label}): {len(ora_res)} terms, {n_sig} significant (FDR <= 0.05).")

            if base_folder:
                log_parts.append(f"**Output folder:** `{tmp}`")
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
    <li>Term names differ across databases. The concordance analysis uses
    normalized matching (stripping prefixes and formatting), but some overlaps
    may still be missed when databases use very different naming conventions.</li>
    </ul>
    </div>
        """, unsafe_allow_html=True)

    splicing = st.session_state.splicing_data

    # ---- Tabs ----
    st.markdown("")
    available_tabs = ["Overview", "DEG Explorer"]
    if gw_up is not None or gw_down is not None:
        available_tabs.append("GeneWalk Comparison")
    if gsea_res is not None and not gsea_res.empty:
        available_tabs.append("GSEA Results")
    if ora_up is not None or ora_down is not None:
        available_tabs.append("ORA Results")
    if splicing is not None and not splicing.empty:
        available_tabs.append("Splicing")
    available_tabs.append("Cross-Method")
    available_tabs.append("Gene Investigator")
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
                key="plt_1",
            )

        if gw_up is not None and gw_down is not None:
            st.plotly_chart(
                direction_volcano(gw_up, gw_down, padj_col=gw_padj_col,
                                  padj_threshold=gw_padj_threshold),
                width="stretch",
                key="plt_2",
            )

    # ---- Tab: DEG Explorer ----
    with tabs[tab_idx]:
        tab_idx += 1

        if deg is not None and not deg.empty:
            # ---- Dynamic column filters ----
            st.markdown(
                '<div class="info-tip">Use these filters to explore your DESeq2 '
                "results. All plots and the table below update in real time. "
                "Filtered genes are used to compute summary statistics.</div>",
                unsafe_allow_html=True,
            )

            # Build filter controls dynamically for all numeric columns
            numeric_cols = deg.select_dtypes(include="number").columns.tolist()
            # Separate known columns from generic numeric ones
            known_filters = []
            if "baseMean" in numeric_cols:
                known_filters.append("baseMean")
            if "log2fc" in numeric_cols:
                known_filters.append("log2fc")
            if "lfcSE" in numeric_cols:
                known_filters.append("lfcSE")
            if "stat" in numeric_cols:
                known_filters.append("stat")
            if "pvalue" in numeric_cols:
                known_filters.append("pvalue")
            if "padj" in numeric_cols:
                known_filters.append("padj")
            other_numeric = [c for c in numeric_cols if c not in known_filters]

            st.markdown('<p class="section-header">Filters</p>',
                        unsafe_allow_html=True)

            deg_filtered = deg.copy()

            # Known DESeq2 column filters with sensible defaults
            fcols = st.columns(min(len(known_filters) + len(other_numeric), 4) or 1)
            fi = 0

            if "baseMean" in known_filters:
                col_min = float(deg["baseMean"].min()) if not deg["baseMean"].isna().all() else 0.0
                col_max = float(deg["baseMean"].max()) if not deg["baseMean"].isna().all() else 100000.0
                bm_min = fcols[fi % len(fcols)].number_input(
                    "baseMean >=",
                    min_value=0.0,
                    max_value=col_max,
                    value=0.0,
                    step=10.0,
                    help="Filter out lowly-expressed genes. DESeq2 recommends "
                         "baseMean >= 10-20 for reliable results.",
                )
                deg_filtered = deg_filtered[deg_filtered["baseMean"] >= bm_min]
                fi += 1

            if "padj" in known_filters:
                padj_max = fcols[fi % len(fcols)].number_input(
                    "padj <=",
                    min_value=0.0,
                    max_value=1.0,
                    value=1.0,
                    step=0.01,
                    format="%.3f",
                    help="Show only genes below this adjusted p-value.",
                )
                deg_filtered = deg_filtered[
                    deg_filtered["padj"].isna() | (deg_filtered["padj"] <= padj_max)
                ]
                fi += 1

            if "log2fc" in known_filters:
                lfc_abs_min = fcols[fi % len(fcols)].number_input(
                    "|log2FC| >=",
                    min_value=0.0,
                    max_value=20.0,
                    value=0.0,
                    step=0.1,
                    format="%.2f",
                    help="Show only genes with absolute fold-change above this.",
                )
                deg_filtered = deg_filtered[deg_filtered["log2fc"].abs() >= lfc_abs_min]
                fi += 1

            if "pvalue" in known_filters:
                pv_max = fcols[fi % len(fcols)].number_input(
                    "pvalue <=",
                    min_value=0.0,
                    max_value=1.0,
                    value=1.0,
                    step=0.01,
                    format="%.4f",
                    help="Filter by raw (unadjusted) p-value.",
                )
                deg_filtered = deg_filtered[
                    deg_filtered["pvalue"].isna() | (deg_filtered["pvalue"] <= pv_max)
                ]
                fi += 1

            # Extra numeric column filters in an expander
            if other_numeric:
                with st.expander("Additional column filters", expanded=False):
                    extra_cols = st.columns(min(len(other_numeric), 3))
                    for j, col_name in enumerate(other_numeric):
                        col_vals = deg_filtered[col_name].dropna()
                        if col_vals.empty:
                            continue
                        cmin, cmax = float(col_vals.min()), float(col_vals.max())
                        if cmin == cmax:
                            continue
                        rng = extra_cols[j % len(extra_cols)].slider(
                            col_name,
                            min_value=cmin,
                            max_value=cmax,
                            value=(cmin, cmax),
                            key=f"deg_filter_{col_name}",
                        )
                        deg_filtered = deg_filtered[
                            deg_filtered[col_name].between(rng[0], rng[1])
                            | deg_filtered[col_name].isna()
                        ]

            # Gene search
            gene_search = st.text_input(
                "Search genes",
                placeholder="e.g. BRAF, KRAS, MYC",
                help="Comma-separated list. Shows only matching genes in the table.",
            )
            if gene_search.strip():
                search_terms = [g.strip().upper() for g in gene_search.split(",") if g.strip()]
                if "gene" in deg_filtered.columns:
                    deg_filtered = deg_filtered[
                        deg_filtered["gene"].str.upper().isin(search_terms)
                    ]

            # ---- Summary metrics ----
            st.markdown("")
            deg_counts = deg_summary_counts(
                deg_filtered, "log2fc", "padj", fc_threshold, sig_threshold,
            )
            dm = st.columns(6)
            dm[0].metric("Total Genes", f"{deg_counts['total_genes']:,}")
            dm[1].metric("Significant", f"{deg_counts['significant']:,}")
            dm[2].metric("Up-regulated", f"{deg_counts['up_regulated']:,}")
            dm[3].metric("Down-regulated", f"{deg_counts['down_regulated']:,}")
            dm[4].metric("Not Significant", f"{deg_counts['not_significant']:,}")
            dm[5].metric("% Significant", f"{deg_counts['pct_significant']}%")

            # ---- Sub-tabs for plots ----
            deg_sub = st.tabs([
                "Volcano", "MA Plot", "Top Genes",
                "P-value Dist.", "Category", "Diagnostics",
            ])

            with deg_sub[0]:
                st.plotly_chart(
                    deg_overview_volcano(deg_filtered, fc_threshold, sig_threshold),
                    width="stretch",
                    key="plt_3",
                )

            with deg_sub[1]:
                if "baseMean" in deg_filtered.columns:
                    st.plotly_chart(
                        ma_plot(deg_filtered, "baseMean", "log2fc", "padj",
                                fc_threshold, sig_threshold),
                        width="stretch",
                        key="plt_4",
                    )
                else:
                    st.info("MA plot requires a baseMean column. Map it in the sidebar under 'Additional DESeq2 columns'.")

            with deg_sub[2]:
                top_n_deg = st.slider("Number of top genes", 10, 100, 30, key="deg_top_n")
                st.plotly_chart(
                    top_genes_bar(deg_filtered, "log2fc", "padj",
                                  sig_threshold, top_n_deg),
                    width="stretch",
                    key="plt_5",
                )

            with deg_sub[3]:
                pv_col1, pv_col2 = st.columns(2)
                with pv_col1:
                    if "pvalue" in deg_filtered.columns:
                        st.plotly_chart(
                            pvalue_histogram(deg_filtered, "pvalue"),
                            width="stretch",
                            key="plt_6",
                        )
                    else:
                        st.info("Raw p-value column not mapped.")
                with pv_col2:
                    st.plotly_chart(
                        padj_histogram(deg_filtered, "padj"),
                        width="stretch",
                        key="plt_7",
                    )

            with deg_sub[4]:
                cat_col1, cat_col2 = st.columns(2)
                with cat_col1:
                    st.plotly_chart(
                        deg_category_pie(deg_filtered, "log2fc", "padj",
                                         fc_threshold, sig_threshold),
                        width="stretch",
                        key="plt_8",
                    )
                with cat_col2:
                    if "baseMean" in deg_filtered.columns:
                        st.plotly_chart(
                            basemean_distribution(deg_filtered, "baseMean"),
                            width="stretch",
                            key="plt_9",
                        )

            with deg_sub[5]:
                diag_col1, diag_col2 = st.columns(2)
                with diag_col1:
                    if "stat" in deg_filtered.columns:
                        st.plotly_chart(
                            lfc_vs_stat_scatter(deg_filtered, "log2fc", "stat",
                                                "padj", sig_threshold),
                            width="stretch",
                            key="plt_10",
                        )
                    else:
                        st.info("Wald statistic column not mapped.")
                with diag_col2:
                    if "baseMean" in deg_filtered.columns and "lfcSE" in deg_filtered.columns:
                        st.plotly_chart(
                            lfc_se_scatter(deg_filtered, "baseMean", "lfcSE"),
                            width="stretch",
                            key="plt_11",
                        )
                    else:
                        st.info("Need both baseMean and lfcSE for this plot.")

            # ---- Filtered data table ----
            st.markdown(f"**{len(deg_filtered):,}** genes after filtering")
            st.dataframe(
                deg_filtered.sort_values("padj", ascending=True)
                if "padj" in deg_filtered.columns
                else deg_filtered,
                width="stretch",
                height=400,
            )
            csv_buf_deg = io.StringIO()
            deg_filtered.to_csv(csv_buf_deg, index=False)
            st.download_button(
                "Download filtered DEG table",
                csv_buf_deg.getvalue(),
                file_name="filtered_deg_table.csv",
                mime="text/csv",
                key="dl_filtered_deg",
            )
        else:
            st.info("No DEG table loaded. Upload one in the sidebar or run the analysis first.")

    # ---- Tab: GeneWalk Comparison ----
    if "GeneWalk Comparison" in available_tabs:
        with tabs[tab_idx]:
            tab_idx += 1

            gw_sub = st.tabs(["Up-Regulated", "Down-Regulated", "Shared GO Terms"])

            with gw_sub[0]:
                if gw_up is not None and not gw_up.empty:
                    st.markdown(f"**{len(up_genes)} up-regulated genes** analyzed by GeneWalk")
                    render_dashboard(gw_up, key_prefix="gw_up_")
                else:
                    st.info("No GeneWalk results for up-regulated genes.")

            with gw_sub[1]:
                if gw_down is not None and not gw_down.empty:
                    st.markdown(f"**{len(down_genes)} down-regulated genes** analyzed by GeneWalk")
                    render_dashboard(gw_down, key_prefix="gw_down_")
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
                        st.plotly_chart(shared_terms_bar(shared), width="stretch", key="plt_12")
                        st.dataframe(shared, width="stretch", height=400)
                    else:
                        st.info("No GO terms are significant in both up- and down-regulated gene sets at the current threshold.")
                else:
                    st.info("Need both up and down GeneWalk results to compare.")

    # ---- Tab: GSEA Results ----
    if "GSEA Results" in available_tabs:
        with tabs[tab_idx]:
            tab_idx += 1

            gsea_sub = st.tabs(["NES Bar Chart", "Dot Plot", "Leading Edge Genes", "Full Table"])

            with gsea_sub[0]:
                top_n_nes = st.slider("Top pathways", 10, 60, 30, key="nes_top")
                st.plotly_chart(
                    nes_bar_chart(gsea_res, fdr_threshold=gsea_fdr_threshold,
                                  top_n=top_n_nes),
                    width="stretch",
                    key="plt_13",
                )

            with gsea_sub[1]:
                st.plotly_chart(
                    gsea_dot_plot(gsea_res, fdr_threshold=gsea_fdr_threshold),
                    width="stretch",
                    key="plt_14",
                )

            with gsea_sub[2]:
                st.markdown(
                    '<div class="info-tip">Leading edge genes are the core subset '
                    "driving each pathway's enrichment score. These are the genes "
                    "that appear before the running enrichment score reaches its "
                    "maximum (for up-regulated) or minimum (for down-regulated) "
                    "peak. They represent the most biologically relevant genes "
                    "for each pathway.</div>",
                    unsafe_allow_html=True,
                )
                le_table = gsea_leading_edge_table(
                    gsea_res, fdr_threshold=gsea_fdr_threshold, top_n=30,
                )
                if not le_table.empty:
                    st.dataframe(le_table, width="stretch", height=500)
                    csv_le = io.StringIO()
                    le_table.to_csv(csv_le, index=False)
                    st.download_button(
                        "Download leading edge genes",
                        csv_le.getvalue(),
                        file_name="gsea_leading_edge_genes.csv",
                        mime="text/csv",
                        key="dl_leading_edge",
                    )
                else:
                    st.info("No significant pathways at current FDR threshold.")

            with gsea_sub[3]:
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
                    ora_up_viz = st.tabs(["Bar Chart", "Dot Plot", "Data"])
                    with ora_up_viz[0]:
                        st.plotly_chart(
                            ora_bar_chart(ora_up, label="Up-regulated"),
                            width="stretch",
                            key="plt_15",
                        )
                    with ora_up_viz[1]:
                        st.plotly_chart(
                            ora_dot_plot(ora_up, label="Up-regulated"),
                            width="stretch",
                            key="plt_16",
                        )
                    with ora_up_viz[2]:
                        st.dataframe(ora_up, width="stretch", height=400)
                        csv_ora_up = io.StringIO()
                        ora_up.to_csv(csv_ora_up, index=False)
                        st.download_button(
                            "Download ORA up results",
                            csv_ora_up.getvalue(),
                            file_name="ora_up_results.csv",
                            mime="text/csv",
                            key="dl_ora_up_tab",
                        )
                else:
                    st.info("No ORA results for up-regulated genes.")

            with ora_sub[1]:
                if ora_down is not None and not ora_down.empty:
                    ora_down_viz = st.tabs(["Bar Chart", "Dot Plot", "Data"])
                    with ora_down_viz[0]:
                        st.plotly_chart(
                            ora_bar_chart(ora_down, label="Down-regulated"),
                            width="stretch",
                            key="plt_17",
                        )
                    with ora_down_viz[1]:
                        st.plotly_chart(
                            ora_dot_plot(ora_down, label="Down-regulated"),
                            width="stretch",
                            key="plt_18",
                        )
                    with ora_down_viz[2]:
                        st.dataframe(ora_down, width="stretch", height=400)
                        csv_ora_down = io.StringIO()
                        ora_down.to_csv(csv_ora_down, index=False)
                        st.download_button(
                            "Download ORA down results",
                            csv_ora_down.getvalue(),
                            file_name="ora_down_results.csv",
                            mime="text/csv",
                            key="dl_ora_down_tab",
                        )
                else:
                    st.info("No ORA results for down-regulated genes.")

    # ---- Tab: Splicing ----
    if "Splicing" in available_tabs:
        with tabs[tab_idx]:
            tab_idx += 1

            from genewalk_app.splicing_viz import (
                dpsi_volcano as spl_volcano,
                event_type_summary as spl_event_summary,
                top_splicing_events_bar as spl_top_bar,
                dpsi_distribution as spl_dpsi_dist,
                gene_splicing_detail as spl_gene_detail,
                genes_by_event_count as spl_gene_counts,
            )

            spl_fc1, spl_fc2 = st.columns(2)
            spl_dpsi_thresh = spl_fc1.number_input(
                "|\u0394PSI| threshold", 0.0, 1.0, 0.1, 0.01,
                format="%.2f", key="spl_dpsi_thresh",
            )
            spl_fdr_thresh = spl_fc2.number_input(
                "Splicing FDR", 0.0, 1.0, 0.05, 0.01,
                format="%.3f", key="spl_fdr_thresh",
            )

            spl_filtered = filter_splicing_events(
                splicing, dpsi_threshold=spl_dpsi_thresh, fdr_threshold=spl_fdr_thresh,
            )
            spl_summary = splicing_summary(spl_filtered)

            sm = st.columns(4)
            sm[0].metric("Sig. Splicing Events", f"{spl_summary['total_events']:,}")
            sm[1].metric("Affected Genes", f"{spl_summary['unique_genes']:,}")
            sm[2].metric("Mean |\u0394PSI|", f"{spl_summary['mean_abs_dpsi']:.3f}")
            et_str = ", ".join(f"{k}: {v}" for k, v in spl_summary.get("event_types", {}).items())
            sm[3].metric("Event Types", et_str if et_str else "N/A")

            spl_tabs = st.tabs(["Volcano", "Top Events", "Distributions", "Per-Gene"])

            with spl_tabs[0]:
                st.plotly_chart(
                    spl_volcano(splicing, spl_dpsi_thresh, spl_fdr_thresh),
                    width="stretch",
                    key="plt_19",
                )

            with spl_tabs[1]:
                st.plotly_chart(
                    spl_top_bar(spl_filtered, top_n=30),
                    width="stretch",
                    key="plt_20",
                )

            with spl_tabs[2]:
                dist_c1, dist_c2 = st.columns(2)
                with dist_c1:
                    st.plotly_chart(spl_dpsi_dist(splicing), width="stretch", key="plt_21")
                with dist_c2:
                    st.plotly_chart(spl_event_summary(spl_filtered), width="stretch", key="plt_22")
                st.plotly_chart(spl_gene_counts(spl_filtered), width="stretch", key="plt_23")

            with spl_tabs[3]:
                spl_genes = sorted(spl_filtered["gene"].dropna().unique().tolist())
                if spl_genes:
                    spl_gene = st.selectbox("Gene", spl_genes, key="spl_gene_comp")
                    st.plotly_chart(spl_gene_detail(splicing, spl_gene), width="stretch", key="plt_24")
                else:
                    st.info("No significant splicing events at current thresholds.")

    # ---- Tab: Cross-Method ----
    with tabs[tab_idx]:
        tab_idx += 1

        st.markdown("""
        <div class="info-tip">Terms found by multiple methods (GeneWalk, GSEA, ORA)
        represent the highest-confidence functional findings. Terms are matched
        using normalized names (stripping database prefixes, lowercasing, and
        removing formatting differences), so GeneWalk GO terms can match GSEA
        pathway names even when databases use different naming conventions.</div>
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
            st.plotly_chart(concordance_bar(concordance), width="stretch", key="plt_25")
            st.dataframe(concordance, width="stretch", height=400)
        else:
            st.info(
                "No terms found by multiple methods at the current thresholds. "
                "Try relaxing the FDR/p-value cutoffs, or note that term names "
                "differ across databases (GO vs KEGG vs Reactome)."
            )

    # ---- Tab: Gene Investigator ----
    with tabs[tab_idx]:
        tab_idx += 1

        render_gene_investigator(
            deg=deg,
            gw_up=gw_up,
            gw_down=gw_down,
            gsea_res=gsea_res,
            ora_up=ora_up,
            ora_down=ora_down,
            splicing=splicing,
            gw_padj_col=gw_padj_col,
            gw_padj_threshold=gw_padj_threshold,
            fc_threshold=fc_threshold,
            padj_threshold=sig_threshold,
        )

    # ---- Tab: Data Tables ----
    with tabs[tab_idx]:
        dt_sub_names = ["DEG Table", "GW Up", "GW Down", "GSEA", "ORA Up", "ORA Down"]
        dt_datasets = [
            (deg, "deg_table.csv"),
            (gw_up, "genewalk_up_results.csv"),
            (gw_down, "genewalk_down_results.csv"),
            (gsea_res, "gsea_prerank_results.csv"),
            (ora_up, "ora_up_results.csv"),
            (ora_down, "ora_down_results.csv"),
        ]
        if splicing is not None and not splicing.empty:
            dt_sub_names.append("Splicing")
            dt_datasets.append((splicing, "splicing_results.csv"))
        dt_sub = st.tabs(dt_sub_names)

        for sub_tab, (data, filename) in zip(dt_sub, dt_datasets):
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


# ---------------------------------------------------------------------------
# Standalone entry point — runs when executed via 'streamlit run comparison_app.py'
# When imported by desktop.py, only run_comparison_ui() is called (page config
# is already set by desktop.py).
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    st.set_page_config(
        page_title="GeneWalk + GSEA Comparison",
        page_icon=":dna:",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.markdown(get_custom_css(), unsafe_allow_html=True)
    run_comparison_ui()
