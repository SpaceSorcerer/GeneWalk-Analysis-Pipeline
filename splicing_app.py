"""Standalone RNA Splicing Analysis App.

Upload differential splicing results from vast-tools or rMATS and
explore splicing events with interactive visualizations.

Launch with:  streamlit run splicing_app.py
"""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import streamlit as st

from genewalk_app.splicing import (
    auto_detect_format,
    filter_splicing_events,
    get_gene_splicing_events,
    parse_rmats_file,
    parse_vasttools,
    splicing_summary,
)
from genewalk_app.splicing_viz import (
    dpsi_distribution,
    dpsi_volcano,
    event_type_pie,
    event_type_summary,
    gene_splicing_detail,
    genes_by_event_count,
    psi_distribution,
    top_splicing_events_bar,
)
from genewalk_app.styles import get_custom_css


def run_splicing_ui() -> None:
    """Render the full splicing analysis UI.

    Can be called standalone or embedded inside desktop.py.
    """
    # ---- Session state ----
    if "splicing_data" not in st.session_state:
        st.session_state.splicing_data = None

    # ---- Sidebar ----
    with st.sidebar:
        st.markdown("## Splicing Analysis")
        st.caption("Upload vast-tools or rMATS results to explore differential splicing.")

        st.markdown('<p class="section-header">Input</p>', unsafe_allow_html=True)

        source_tool = st.radio(
            "Source tool",
            ["vast-tools", "rMATS"],
            help="Select which tool produced your splicing results.",
        )

        uploaded = st.file_uploader(
            "Upload splicing results",
            type=["csv", "tsv", "txt", "tab"],
            help=(
                "For vast-tools: upload the diff/compare output table.\n"
                "For rMATS: upload an event-type file (e.g. SE.MATS.JC.txt)."
            ),
        )

        rmats_event_type = "SE"
        if source_tool == "rMATS":
            rmats_event_type = st.selectbox(
                "rMATS event type",
                ["SE", "MXE", "A3SS", "A5SS", "RI"],
                help="Which alternative splicing event type is in this file.",
            )

            st.markdown(
                '<div class="info-tip">For a complete rMATS analysis, upload '
                "each event type file separately and combine results.</div>",
                unsafe_allow_html=True,
            )

        # Additional rMATS file uploads for combined analysis
        additional_rmats = []
        if source_tool == "rMATS":
            with st.expander("Upload additional event types", expanded=False):
                for et in ["MXE", "A3SS", "A5SS", "RI"]:
                    if et != rmats_event_type:
                        extra = st.file_uploader(
                            f"{et} file",
                            type=["csv", "tsv", "txt", "tab"],
                            key=f"rmats_{et}",
                        )
                        if extra:
                            additional_rmats.append((et, extra))

        # ---- Thresholds ----
        st.markdown("---")
        st.markdown('<p class="section-header">Thresholds</p>', unsafe_allow_html=True)

        dpsi_threshold = st.number_input(
            "|\u0394PSI| threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.1,
            step=0.01,
            format="%.2f",
            help="Minimum absolute change in PSI to consider significant. "
                 "0.1 (10%) is a common cutoff; 0.05 is more permissive.",
        )
        fdr_threshold = st.number_input(
            "FDR threshold",
            min_value=0.0,
            max_value=1.0,
            value=0.05,
            step=0.01,
            format="%.3f",
            help="False discovery rate cutoff for significance.",
        )

        st.markdown("---")
        parse_clicked = st.button(
            "Load & Analyze",
            type="primary",
            disabled=uploaded is None,
            use_container_width=True,
        )

    # ---- Parse data ----
    if parse_clicked and uploaded is not None:
        try:
            content = uploaded.getvalue().decode("utf-8")
        except UnicodeDecodeError:
            content = uploaded.getvalue().decode("latin-1")

        sep = "\t" if uploaded.name.endswith((".tsv", ".txt", ".tab")) else ","

        try:
            raw_df = pd.read_csv(io.StringIO(content), sep=sep)

            if source_tool == "vast-tools":
                parsed = parse_vasttools(raw_df)
            else:
                parsed = parse_rmats_file(raw_df, event_type=rmats_event_type)

                # Parse additional event type files
                for et, extra_file in additional_rmats:
                    try:
                        extra_content = extra_file.getvalue().decode("utf-8")
                    except UnicodeDecodeError:
                        extra_content = extra_file.getvalue().decode("latin-1")
                    extra_sep = "\t" if extra_file.name.endswith((".tsv", ".txt", ".tab")) else ","
                    extra_df = pd.read_csv(io.StringIO(extra_content), sep=extra_sep)
                    extra_parsed = parse_rmats_file(extra_df, event_type=et)
                    parsed = pd.concat([parsed, extra_parsed], ignore_index=True)

            st.session_state.splicing_data = parsed
        except (ValueError, pd.errors.EmptyDataError) as exc:
            st.error(f"Could not parse file: {exc}")

    # ---- Hero banner ----
    st.markdown("""
    <div class="hero-banner">
        <div class="hero-badge">Splicing Analysis</div>
        <h1>RNA Splicing Explorer</h1>
        <p>Explore differential splicing events from vast-tools or rMATS.
        Identify alternative splicing changes, classify event types, and
        investigate individual genes.</p>
    </div>
    """, unsafe_allow_html=True)

    # ---- Landing page ----
    if st.session_state.splicing_data is None:
        col_a, col_b, col_c = st.columns(3, gap="medium")
        with col_a:
            st.markdown("""
            <div class="step-card">
                <div class="step-number">1</div>
                <h3>Upload Results</h3>
                <p>Upload a vast-tools diff/compare table or an rMATS
                event-type file (SE, RI, A3SS, A5SS, MXE).</p>
            </div>
            """, unsafe_allow_html=True)
        with col_b:
            st.markdown("""
            <div class="step-card">
                <div class="step-number">2</div>
                <h3>Set Thresholds</h3>
                <p>Configure &Delta;PSI and FDR cutoffs to define
                significant splicing changes.</p>
            </div>
            """, unsafe_allow_html=True)
        with col_c:
            st.markdown("""
            <div class="step-card">
                <div class="step-number">3</div>
                <h3>Explore Events</h3>
                <p>Volcano plots, event type distributions,
                per-gene splicing views, and data export.</p>
            </div>
            """, unsafe_allow_html=True)

        with st.expander("About this app", expanded=False):
            st.markdown("""
<div class="guide-section">
<h4>What is differential splicing?</h4>
<p>While differential expression (DESeq2) measures changes in total RNA
abundance, differential splicing detects changes in <strong>how</strong>
pre-mRNA is processed. A gene can be alternatively spliced without any
change in total expression, or vice versa.</p>

<h4>Event types</h4>
<ul>
<li><strong>SE</strong> &mdash; Skipped Exon (cassette exon): an exon
is either included or excluded from the mature mRNA.</li>
<li><strong>RI</strong> &mdash; Retained Intron: an intron is retained
in the mature mRNA instead of being spliced out.</li>
<li><strong>A3SS</strong> &mdash; Alternative 3' Splice Site: the
3' end of an exon is extended or shortened.</li>
<li><strong>A5SS</strong> &mdash; Alternative 5' Splice Site: the
5' end of an exon is extended or shortened.</li>
<li><strong>MXE</strong> &mdash; Mutually Exclusive Exons: one of two
exons is included, but never both.</li>
<li><strong>MIC</strong> &mdash; Microexon (vast-tools only): very
small exons (3-27 nt) often in neuronal genes.</li>
</ul>

<h4>Key metric: &Delta;PSI</h4>
<p><strong>PSI</strong> (Percent Spliced In) ranges from 0 to 1 and
represents what fraction of transcripts include a particular exon or
splicing event. <strong>&Delta;PSI</strong> is the change between
conditions. A &Delta;PSI of +0.2 means the exon is included in 20%
more transcripts in the treatment vs control.</p>

<h4>Supported tools</h4>
<ul>
<li><strong>vast-tools</strong> &mdash; Upload the output from
<code>vast-tools compare</code> or <code>vast-tools diff</code>.</li>
<li><strong>rMATS</strong> &mdash; Upload individual event-type files
(e.g. <code>SE.MATS.JC.txt</code>). For comprehensive analysis,
upload all five event types.</li>
</ul>
</div>
            """, unsafe_allow_html=True)

        st.markdown(
            '<div class="app-footer">Built with '
            '<a href="https://streamlit.io">Streamlit</a></div>',
            unsafe_allow_html=True,
        )
        return

    # =========================================================================
    # Results dashboard
    # =========================================================================
    data = st.session_state.splicing_data

    # Apply filters
    filtered = filter_splicing_events(
        data, dpsi_threshold=dpsi_threshold, fdr_threshold=fdr_threshold,
    )

    # ---- Summary metrics ----
    summary_all = splicing_summary(data)
    summary_sig = splicing_summary(filtered)

    st.markdown("")
    mc = st.columns(6)
    mc[0].metric("Total Events", f"{summary_all['total_events']:,}")
    mc[1].metric("Significant Events", f"{summary_sig['total_events']:,}")
    mc[2].metric("Total Genes", f"{summary_all['unique_genes']:,}")
    mc[3].metric("Sig. Genes", f"{summary_sig['unique_genes']:,}")
    mc[4].metric("Mean |\u0394PSI| (sig)", f"{summary_sig['mean_abs_dpsi']:.3f}")
    mc[5].metric("Source", data["source"].iloc[0] if "source" in data.columns else "?")

    # ---- Tabs ----
    tabs = st.tabs([
        "Volcano", "Event Types", "Top Events",
        "Distributions", "Per-Gene", "Data Table",
    ])

    # ---- Volcano ----
    with tabs[0]:
        st.plotly_chart(
            dpsi_volcano(data, dpsi_threshold=dpsi_threshold, fdr_threshold=fdr_threshold),
            width="stretch",
            key="spl_volcano",
        )

    # ---- Event Types ----
    with tabs[1]:
        evt_col1, evt_col2 = st.columns(2)
        with evt_col1:
            st.plotly_chart(event_type_summary(filtered), width="stretch", key="spl_evt_bar")
        with evt_col2:
            st.plotly_chart(event_type_pie(filtered), width="stretch", key="spl_evt_pie")

        st.markdown("**All events (before filtering):**")
        all_col1, all_col2 = st.columns(2)
        with all_col1:
            st.plotly_chart(event_type_summary(data), width="stretch", key="spl_evt_bar_all")

    # ---- Top Events ----
    with tabs[2]:
        top_n_spl = st.slider("Number of top events", 10, 60, 30, key="spl_top_n")
        sort_metric = st.selectbox(
            "Sort by", ["dpsi", "fdr"],
            help="Sort events by absolute delta-PSI or statistical significance.",
        )
        st.plotly_chart(
            top_splicing_events_bar(filtered, top_n=top_n_spl, sort_by=sort_metric),
            width="stretch",
            key="spl_top_bar",
        )

    # ---- Distributions ----
    with tabs[3]:
        dist_col1, dist_col2 = st.columns(2)
        with dist_col1:
            st.plotly_chart(psi_distribution(data), width="stretch", key="spl_psi_dist")
        with dist_col2:
            st.plotly_chart(dpsi_distribution(data), width="stretch", key="spl_dpsi_dist")

        st.plotly_chart(genes_by_event_count(filtered, top_n=30), width="stretch", key="spl_gene_count")

    # ---- Per-Gene ----
    with tabs[4]:
        sig_genes = sorted(filtered["gene"].dropna().unique().tolist())
        all_genes = sorted(data["gene"].dropna().unique().tolist())
        gene_pool = sig_genes if sig_genes else all_genes

        if gene_pool:
            selected_gene = st.selectbox(
                "Select gene", gene_pool, key="spl_gene_select",
            )
            if selected_gene:
                gene_events = get_gene_splicing_events(data, selected_gene)
                st.plotly_chart(
                    gene_splicing_detail(data, selected_gene),
                    width="stretch",
                    key="spl_gene_detail",
                )
                st.dataframe(gene_events, width="stretch", height=300)
        else:
            st.info("No genes found in the data.")

    # ---- Data Table ----
    with tabs[5]:
        st.markdown(f"**{len(filtered):,}** significant events (|\u0394PSI| >= {dpsi_threshold}, FDR <= {fdr_threshold})")
        st.dataframe(filtered, width="stretch", height=500)

        csv_buf = io.StringIO()
        filtered.to_csv(csv_buf, index=False)
        st.download_button(
            "Download filtered splicing results",
            csv_buf.getvalue(),
            file_name="splicing_filtered_results.csv",
            mime="text/csv",
            key="dl_spl_filtered",
        )

    # Footer
    st.markdown(
        '<div class="app-footer">Built with '
        '<a href="https://streamlit.io">Streamlit</a></div>',
        unsafe_allow_html=True,
    )


# ── Standalone entry point ──
if __name__ == "__main__":
    st.set_page_config(
        page_title="RNA Splicing Explorer",
        page_icon=":scissors:",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    st.markdown(get_custom_css(), unsafe_allow_html=True)
    run_splicing_ui()
