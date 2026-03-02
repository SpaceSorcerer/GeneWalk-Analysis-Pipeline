"""Gene Investigator: cross-analysis per-gene deep dive.

Collects all evidence for a single gene across DEG, GeneWalk, GSEA,
ORA, and splicing analyses and renders a unified Streamlit view.
"""

from __future__ import annotations

import io

import pandas as pd
import plotly.graph_objects as go
import streamlit as st


# ───────────────────────────────────────────────────────────────────────
# Evidence collection
# ───────────────────────────────────────────────────────────────────────

def collect_all_genes(
    deg: pd.DataFrame | None = None,
    gw_up: pd.DataFrame | None = None,
    gw_down: pd.DataFrame | None = None,
    gsea_res: pd.DataFrame | None = None,
    ora_up: pd.DataFrame | None = None,
    ora_down: pd.DataFrame | None = None,
    splicing: pd.DataFrame | None = None,
) -> list[str]:
    """Collect a sorted, deduplicated list of all gene names across analyses."""
    genes: set[str] = set()

    if deg is not None and "gene" in deg.columns:
        genes.update(deg["gene"].dropna().unique())

    for gw in (gw_up, gw_down):
        if gw is not None and "hgnc_symbol" in gw.columns:
            genes.update(gw["hgnc_symbol"].dropna().unique())

    # GSEA / ORA — extract genes from leading edge / gene lists
    for enrich_df in (gsea_res, ora_up, ora_down):
        if enrich_df is not None:
            for col in ("lead_genes", "Lead_genes", "genes", "Genes"):
                if col in enrich_df.columns:
                    for val in enrich_df[col].dropna():
                        for g in str(val).replace(",", ";").split(";"):
                            g = g.strip()
                            if g and g.lower() != "nan":
                                genes.add(g)

    if splicing is not None and "gene" in splicing.columns:
        genes.update(splicing["gene"].dropna().unique())

    return sorted(genes)


def _get_deg_evidence(gene: str, deg: pd.DataFrame | None) -> dict | None:
    """Get DEG evidence for a gene."""
    if deg is None or "gene" not in deg.columns:
        return None
    row = deg[deg["gene"].str.upper() == gene.upper()]
    if row.empty:
        return None
    best = row.sort_values("padj", ascending=True).iloc[0] if "padj" in row.columns else row.iloc[0]
    result = {"gene": gene}
    for col in ("log2fc", "padj", "pvalue", "baseMean", "lfcSE", "stat"):
        if col in best.index:
            result[col] = best[col]
    return result


def _get_gw_evidence(
    gene: str,
    gw_df: pd.DataFrame | None,
    padj_col: str = "global_padj",
    padj_threshold: float = 0.05,
) -> pd.DataFrame:
    """Get GeneWalk GO terms for a gene."""
    if gw_df is None or "hgnc_symbol" not in gw_df.columns:
        return pd.DataFrame()
    gdf = gw_df[gw_df["hgnc_symbol"].str.upper() == gene.upper()].copy()
    if gdf.empty:
        return gdf
    # Mark significance
    if padj_col in gdf.columns:
        gdf["significant"] = gdf[padj_col] <= padj_threshold
        gdf = gdf.sort_values(padj_col, ascending=True)
    return gdf


def _get_gsea_pathways(gene: str, gsea_res: pd.DataFrame | None) -> pd.DataFrame:
    """Find GSEA pathways where this gene appears in the leading edge."""
    if gsea_res is None or gsea_res.empty:
        return pd.DataFrame()

    lead_col = None
    for candidate in ("lead_genes", "Lead_genes", "genes", "ledge_genes"):
        if candidate in gsea_res.columns:
            lead_col = candidate
            break

    if lead_col is None:
        return pd.DataFrame()

    gene_upper = gene.upper()
    mask = gsea_res[lead_col].apply(
        lambda x: gene_upper in [
            g.strip().upper()
            for g in str(x).replace(",", ";").split(";")
        ] if pd.notna(x) else False
    )

    result = gsea_res[mask].copy()
    if not result.empty:
        cols = [c for c in ["term", "nes", "fdr", "gene_set_library"] if c in result.columns]
        if cols:
            result = result[cols]
        if "nes" in result.columns:
            result = result.sort_values("nes", key=abs, ascending=False)
    return result


def _get_ora_terms(gene: str, ora_df: pd.DataFrame | None) -> pd.DataFrame:
    """Find ORA terms where this gene is in the overlap."""
    if ora_df is None or ora_df.empty:
        return pd.DataFrame()

    gene_col = None
    for candidate in ("genes", "Genes"):
        if candidate in ora_df.columns:
            gene_col = candidate
            break

    if gene_col is None:
        return pd.DataFrame()

    gene_upper = gene.upper()
    mask = ora_df[gene_col].apply(
        lambda x: gene_upper in [
            g.strip().upper()
            for g in str(x).replace(",", ";").split(";")
        ] if pd.notna(x) else False
    )

    result = ora_df[mask].copy()
    if not result.empty:
        cols = [c for c in ["term", "fdr", "overlap", "combined_score"] if c in result.columns]
        if cols:
            result = result[cols]
        if "fdr" in result.columns:
            result = result.sort_values("fdr", ascending=True)
    return result


def _get_splicing_evidence(gene: str, splicing: pd.DataFrame | None) -> pd.DataFrame:
    """Get splicing events for a gene."""
    if splicing is None or splicing.empty or "gene" not in splicing.columns:
        return pd.DataFrame()
    return splicing[splicing["gene"].str.upper() == gene.upper()].copy()


# ───────────────────────────────────────────────────────────────────────
# Mini volcano highlight
# ───────────────────────────────────────────────────────────────────────

def _mini_volcano(
    deg: pd.DataFrame,
    gene: str,
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> go.Figure:
    """Small volcano plot with the selected gene highlighted."""
    import numpy as np

    plot = deg.dropna(subset=["log2fc"]).copy()
    if "padj" not in plot.columns:
        return go.Figure()

    plot["neg_log10_padj"] = -np.log10(plot["padj"].clip(lower=1e-300))

    fig = go.Figure()

    # Background genes
    fig.add_trace(go.Scatter(
        x=plot["log2fc"],
        y=plot["neg_log10_padj"],
        mode="markers",
        marker=dict(size=4, color="#B0BEC5", opacity=0.4),
        name="Other genes",
        hoverinfo="skip",
    ))

    # Highlight the selected gene
    gene_rows = plot[plot["gene"].str.upper() == gene.upper()]
    if not gene_rows.empty:
        fig.add_trace(go.Scatter(
            x=gene_rows["log2fc"],
            y=gene_rows["neg_log10_padj"],
            mode="markers+text",
            marker=dict(size=14, color="#D9534F", line=dict(color="white", width=2)),
            text=[gene],
            textposition="top center",
            name=gene,
        ))

    fig.add_vline(x=fc_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.add_vline(x=-fc_threshold, line_dash="dash", line_color="#888", line_width=0.8)
    fig.add_hline(
        y=-np.log10(padj_threshold), line_dash="dash", line_color="#888", line_width=0.8,
    )
    fig.update_layout(
        title=f"{gene} on Volcano Plot",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(padj)",
        template="plotly_white",
        height=350,
        showlegend=False,
        margin=dict(t=40, b=40),
    )
    return fig


# ───────────────────────────────────────────────────────────────────────
# Main render function
# ───────────────────────────────────────────────────────────────────────

def render_gene_investigator(
    deg: pd.DataFrame | None = None,
    gw_up: pd.DataFrame | None = None,
    gw_down: pd.DataFrame | None = None,
    gsea_res: pd.DataFrame | None = None,
    ora_up: pd.DataFrame | None = None,
    ora_down: pd.DataFrame | None = None,
    splicing: pd.DataFrame | None = None,
    gw_padj_col: str = "global_padj",
    gw_padj_threshold: float = 0.05,
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> None:
    """Render the Gene Investigator panel in Streamlit.

    Allows the user to select any gene and see all evidence for it
    across every loaded analysis.
    """
    all_genes = collect_all_genes(deg, gw_up, gw_down, gsea_res, ora_up, ora_down, splicing)

    if not all_genes:
        st.info("No gene data loaded. Run an analysis first.")
        return

    st.markdown(
        '<div class="info-tip">Select any gene to view all evidence across '
        "DEG, GeneWalk, GSEA, ORA, and splicing analyses. Genes appear here "
        "if they are found in <strong>any</strong> loaded analysis.</div>",
        unsafe_allow_html=True,
    )

    col_sel, col_info = st.columns([2, 1])
    with col_sel:
        gene = st.selectbox(
            "Select a gene to investigate",
            all_genes,
            key="gene_investigator_select",
            help="All genes found across any loaded analysis.",
        )
    with col_info:
        st.metric("Total genes available", f"{len(all_genes):,}")

    if not gene:
        return

    # ---- Collect evidence ----
    deg_ev = _get_deg_evidence(gene, deg)
    gw_up_ev = _get_gw_evidence(gene, gw_up, gw_padj_col, gw_padj_threshold)
    gw_down_ev = _get_gw_evidence(gene, gw_down, gw_padj_col, gw_padj_threshold)
    gsea_pathways = _get_gsea_pathways(gene, gsea_res)
    ora_up_terms = _get_ora_terms(gene, ora_up)
    ora_down_terms = _get_ora_terms(gene, ora_down)
    splice_ev = _get_splicing_evidence(gene, splicing)

    # Count evidence sources
    sources = []
    if deg_ev:
        sources.append("DEG")
    if not gw_up_ev.empty:
        sources.append("GW-Up")
    if not gw_down_ev.empty:
        sources.append("GW-Down")
    if not gsea_pathways.empty:
        sources.append("GSEA")
    if not ora_up_terms.empty:
        sources.append("ORA-Up")
    if not ora_down_terms.empty:
        sources.append("ORA-Down")
    if not splice_ev.empty:
        sources.append("Splicing")

    st.markdown(f"### {gene}")
    if sources:
        st.markdown(f"Found in: **{', '.join(sources)}**")
    else:
        st.warning(f"No evidence found for {gene} in any loaded analysis.")
        return

    # ---- DEG Panel ----
    if deg_ev:
        st.markdown("#### Differential Expression")
        deg_cols = st.columns(4)
        if "log2fc" in deg_ev and pd.notna(deg_ev["log2fc"]):
            lfc = deg_ev["log2fc"]
            direction = "Up" if lfc > 0 else "Down"
            deg_cols[0].metric("log2FC", f"{lfc:.3f}", delta=direction)
        if "padj" in deg_ev and pd.notna(deg_ev["padj"]):
            deg_cols[1].metric("padj", f"{deg_ev['padj']:.2e}")
        elif "padj" in deg_ev:
            deg_cols[1].metric("padj", "N/A")
        if "baseMean" in deg_ev and pd.notna(deg_ev["baseMean"]):
            deg_cols[2].metric("baseMean", f"{deg_ev['baseMean']:.1f}")
        if "stat" in deg_ev and pd.notna(deg_ev["stat"]):
            deg_cols[3].metric("Wald stat", f"{deg_ev['stat']:.2f}")

        if deg is not None and "log2fc" in deg.columns and "padj" in deg.columns:
            st.plotly_chart(
                _mini_volcano(deg, gene, fc_threshold, padj_threshold),
                width="stretch",
            )

    # ---- GeneWalk Panel ----
    gw_has = not gw_up_ev.empty or not gw_down_ev.empty
    if gw_has:
        st.markdown("#### GeneWalk: Context-Specific GO Terms")

        for label, gw_ev in [("Up-regulated set", gw_up_ev), ("Down-regulated set", gw_down_ev)]:
            if gw_ev.empty:
                continue
            sig_count = gw_ev["significant"].sum() if "significant" in gw_ev.columns else "?"
            st.markdown(f"**{label}**: {sig_count} significant GO terms")
            display_cols = [c for c in ["go_name", "go_domain", "sim", gw_padj_col]
                           if c in gw_ev.columns]
            if display_cols:
                st.dataframe(
                    gw_ev[display_cols].head(15),
                    width="stretch",
                    height=min(300, len(gw_ev[display_cols].head(15)) * 35 + 40),
                )

    # ---- GSEA Panel ----
    if not gsea_pathways.empty:
        st.markdown("#### GSEA: Enriched Pathways Containing This Gene")
        st.markdown(
            f"**{len(gsea_pathways)}** enriched pathways include {gene} "
            "in their leading edge"
        )
        st.dataframe(
            gsea_pathways.head(15),
            width="stretch",
            height=min(300, len(gsea_pathways.head(15)) * 35 + 40),
        )

    # ---- ORA Panel ----
    ora_has = not ora_up_terms.empty or not ora_down_terms.empty
    if ora_has:
        st.markdown("#### ORA: Over-Represented Terms Containing This Gene")
        for label, ora_terms in [("Up-regulated ORA", ora_up_terms),
                                  ("Down-regulated ORA", ora_down_terms)]:
            if ora_terms.empty:
                continue
            st.markdown(f"**{label}**: {len(ora_terms)} terms contain {gene}")
            st.dataframe(
                ora_terms.head(15),
                width="stretch",
                height=min(300, len(ora_terms.head(15)) * 35 + 40),
            )

    # ---- Splicing Panel ----
    if not splice_ev.empty:
        st.markdown("#### Splicing: Differential Splicing Events")

        spl_cols = st.columns(3)
        spl_cols[0].metric("Splicing Events", f"{len(splice_ev)}")
        if "event_type" in splice_ev.columns:
            types = splice_ev["event_type"].value_counts().to_dict()
            type_str = ", ".join(f"{k}: {v}" for k, v in types.items())
            spl_cols[1].metric("Event Types", type_str)
        if "dpsi" in splice_ev.columns and splice_ev["dpsi"].notna().any():
            max_dpsi = splice_ev["dpsi"].abs().max()
            spl_cols[2].metric("Max |\u0394PSI|", f"{max_dpsi:.3f}")

        # Show per-gene splicing detail plot
        from genewalk_app.splicing_viz import gene_splicing_detail
        st.plotly_chart(gene_splicing_detail(splice_ev, gene), width="stretch")

        display_spl_cols = [c for c in ["event_id", "event_type", "dpsi", "psi_1", "psi_2", "fdr", "coord"]
                           if c in splice_ev.columns]
        st.dataframe(
            splice_ev[display_spl_cols] if display_spl_cols else splice_ev,
            width="stretch",
            height=min(300, len(splice_ev) * 35 + 40),
        )

    # ---- Export ----
    with st.expander(f"Export all evidence for {gene}", expanded=False):
        export_parts = []
        if deg_ev:
            export_parts.append(f"=== Differential Expression ===\n{pd.DataFrame([deg_ev]).to_csv(index=False)}")
        for label, gw_ev in [("GeneWalk Up", gw_up_ev), ("GeneWalk Down", gw_down_ev)]:
            if not gw_ev.empty:
                export_parts.append(f"=== {label} GO Terms ===\n{gw_ev.to_csv(index=False)}")
        if not gsea_pathways.empty:
            export_parts.append(f"=== GSEA Pathways ===\n{gsea_pathways.to_csv(index=False)}")
        for label, ora_ev in [("ORA Up", ora_up_terms), ("ORA Down", ora_down_terms)]:
            if not ora_ev.empty:
                export_parts.append(f"=== {label} Terms ===\n{ora_ev.to_csv(index=False)}")
        if not splice_ev.empty:
            export_parts.append(f"=== Splicing Events ===\n{splice_ev.to_csv(index=False)}")

        if export_parts:
            full_export = "\n\n".join(export_parts)
            st.download_button(
                f"Download all {gene} evidence",
                full_export,
                file_name=f"{gene}_investigation.txt",
                mime="text/plain",
                key=f"dl_investigate_{gene}",
            )
