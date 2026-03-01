"""DEG table parsing, gene list splitting, and cross-analysis comparison."""

from __future__ import annotations

import pandas as pd


def detect_deg_columns(df: pd.DataFrame) -> dict[str, str | None]:
    """Auto-detect gene, log2FC, p-value, and DESeq2-specific columns.

    Returns a dict with keys 'gene', 'log2fc', 'padj', 'basemean',
    'lfcse', 'stat', 'pvalue' mapped to the best matching column name,
    or *None* if no match is found.
    """
    cols = {c.lower(): c for c in df.columns}

    gene_candidates = [
        "gene", "gene_symbol", "hgnc_symbol", "symbol", "gene_name",
        "geneid", "gene_id", "name", "external_gene_name",
    ]
    log2fc_candidates = [
        "log2foldchange", "log2fc", "logfc", "log2_fold_change",
        "lfc", "fc", "foldchange", "fold_change",
    ]
    padj_candidates = [
        "padj", "adj.p.val", "adj_pval", "fdr", "q_value", "qvalue",
        "p_adj", "adjusted_pvalue", "adjpvalue", "bh",
    ]
    basemean_candidates = [
        "basemean", "base_mean", "meanexpression", "averageexpression",
        "aveexpr",
    ]
    lfcse_candidates = [
        "lfcse", "lfc_se", "logfc_se", "se", "standarderror",
    ]
    stat_candidates = [
        "stat", "statistic", "wald_statistic", "t_statistic",
        "test_statistic",
    ]
    pvalue_candidates = [
        "pvalue", "p_value", "pval", "p-value", "p.value",
        "nominal_pvalue", "rawpvalue",
    ]

    def _find(candidates: list[str]) -> str | None:
        for c in candidates:
            if c in cols:
                return cols[c]
        return None

    return {
        "gene": _find(gene_candidates),
        "log2fc": _find(log2fc_candidates),
        "padj": _find(padj_candidates),
        "basemean": _find(basemean_candidates),
        "lfcse": _find(lfcse_candidates),
        "stat": _find(stat_candidates),
        "pvalue": _find(pvalue_candidates),
    }


def parse_deg_table(
    df: pd.DataFrame,
    gene_col: str,
    log2fc_col: str,
    padj_col: str,
) -> pd.DataFrame:
    """Validate and normalize a DEG table to standard columns.

    Returns a DataFrame with columns: gene, log2fc, padj (plus any
    additional columns from the original table).

    Genes with NA padj are **preserved** because:
    - DESeq2 sets padj=NA for low-count genes, outliers, and genes
      removed by independent filtering.
    - GSEA prerank needs the full ranked list (all genes with a log2FC),
      not just genes with valid adjusted p-values.
    - Only genes missing a gene name or log2FC are dropped.
    """
    missing = [
        name for name, col in [("gene", gene_col), ("log2fc", log2fc_col), ("padj", padj_col)]
        if col not in df.columns
    ]
    if missing:
        raise ValueError(f"Missing columns: {', '.join(missing)}")

    out = df.copy()
    out = out.rename(columns={
        gene_col: "gene",
        log2fc_col: "log2fc",
        padj_col: "padj",
    })
    out["log2fc"] = pd.to_numeric(out["log2fc"], errors="coerce")
    out["padj"] = pd.to_numeric(out["padj"], errors="coerce")
    # Only drop rows missing gene name or log2FC — keep NA padj rows
    # so GSEA prerank gets the full ranked list.
    out = out.dropna(subset=["gene", "log2fc"])
    out["gene"] = out["gene"].astype(str).str.strip()
    out = out[out["gene"] != ""]
    return out.reset_index(drop=True)


def split_deg_lists(
    deg: pd.DataFrame,
    fc_threshold: float = 1.0,
    padj_threshold: float = 0.05,
) -> tuple[list[str], list[str]]:
    """Split a normalized DEG table into up- and down-regulated gene lists.

    Parameters
    ----------
    deg : pd.DataFrame
        Must have columns 'gene', 'log2fc', 'padj' (from *parse_deg_table*).
    fc_threshold : float
        Absolute log2 fold-change threshold.
    padj_threshold : float
        Adjusted p-value threshold.

    Returns
    -------
    (up_genes, down_genes) : tuple[list[str], list[str]]
    """
    # NA padj values are treated as not significant (they represent
    # genes DESeq2 could not test, e.g. low counts or outliers).
    sig = deg[deg["padj"].notna() & (deg["padj"] <= padj_threshold)]
    up = sig[sig["log2fc"] >= fc_threshold]["gene"].unique().tolist()
    down = sig[sig["log2fc"] <= -fc_threshold]["gene"].unique().tolist()
    return up, down


def make_ranked_list(deg: pd.DataFrame) -> pd.DataFrame:
    """Create a ranked gene list for GSEA prerank from a DEG table.

    Returns a two-column DataFrame (gene, score) sorted by score descending.

    Deduplication strategy: when a gene appears multiple times, keep the
    entry with the smallest padj (most significant).  If padj is not
    available, fall back to keeping the entry with the largest |log2FC|.
    This avoids inflating ranking metrics with outlier fold-change values.
    """
    cols = ["gene", "log2fc"]
    has_padj = "padj" in deg.columns
    if has_padj:
        cols.append("padj")

    rnk = deg[cols].copy()
    rnk = rnk.dropna(subset=["gene", "log2fc"])

    if has_padj:
        # Prefer the most statistically significant entry per gene
        rnk = rnk.sort_values("padj", ascending=True, na_position="last")
        rnk = rnk.drop_duplicates(subset="gene", keep="first")
        rnk = rnk.drop(columns="padj")
    else:
        rnk["abs_score"] = rnk["log2fc"].abs()
        rnk = rnk.sort_values("abs_score", ascending=False).drop_duplicates(
            subset="gene", keep="first"
        )
        rnk = rnk.drop(columns="abs_score")

    rnk = rnk.sort_values("log2fc", ascending=False)
    return rnk.reset_index(drop=True)


def shared_go_terms(
    gw_up: pd.DataFrame,
    gw_down: pd.DataFrame,
    padj_col: str = "global_padj",
    padj_threshold: float = 0.05,
) -> pd.DataFrame:
    """Find GO terms significant in both the up and down GeneWalk results.

    Returns a DataFrame with columns: go_name, go_id, go_domain,
    best_padj_up, best_padj_down.
    """
    def _sig_terms(df: pd.DataFrame) -> pd.DataFrame:
        if padj_col not in df.columns or "go_name" not in df.columns:
            return pd.DataFrame(columns=["go_name", "go_id", "go_domain", "best_padj"])
        sig = df[df[padj_col] <= padj_threshold]
        if sig.empty:
            return pd.DataFrame(columns=["go_name", "go_id", "go_domain", "best_padj"])
        agg_cols = {"best_padj": (padj_col, "min")}
        group_cols = ["go_name"]
        if "go_id" in sig.columns:
            group_cols.append("go_id")
        if "go_domain" in sig.columns:
            group_cols.append("go_domain")
        return (
            sig.groupby(group_cols, dropna=False)
            .agg(**agg_cols)
            .reset_index()
        )

    up_terms = _sig_terms(gw_up)
    down_terms = _sig_terms(gw_down)

    if up_terms.empty or down_terms.empty:
        return pd.DataFrame()

    merged = up_terms.merge(
        down_terms,
        on=[c for c in ["go_name", "go_id", "go_domain"]
            if c in up_terms.columns and c in down_terms.columns],
        suffixes=("_up", "_down"),
    )
    return merged.sort_values("best_padj_up", ascending=True).reset_index(drop=True)


def _normalize_term(term: str) -> str:
    """Normalize a term name for fuzzy matching across databases.

    Strips common prefixes (GO_, KEGG_, REACTOME_, HALLMARK_, WP_),
    converts to lowercase, and removes underscores/hyphens for
    comparison.
    """
    import re
    t = term.lower()
    # Strip common database prefixes
    t = re.sub(
        r"^(go_biological_process_|go_molecular_function_|"
        r"go_cellular_component_|kegg_|reactome_|hallmark_|wp_|"
        r"go:|hsa\d+_)",
        "", t,
    )
    # Replace underscores, hyphens, extra spaces with single space
    t = re.sub(r"[_\-]+", " ", t).strip()
    # Remove parenthetical qualifiers like "(GO:0001234)"
    t = re.sub(r"\s*\(go:\d+\)\s*", "", t)
    return t


def cross_method_concordance(
    gw_results: pd.DataFrame,
    gsea_results: pd.DataFrame,
    ora_results: pd.DataFrame,
    gw_padj_col: str = "global_padj",
    gw_padj_threshold: float = 0.05,
    gsea_fdr_threshold: float = 0.25,
    ora_fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Build a concordance table of terms found by multiple methods.

    Uses normalized term matching: strips database prefixes, lowercases,
    and compares cleaned names so that e.g. GeneWalk's
    "cell proliferation" can match GSEA's
    "GO_BIOLOGICAL_PROCESS_CELL_PROLIFERATION" or ORA's
    "Cell Proliferation (GO:0008283)".

    Returns a DataFrame with columns: term, found_in (list of methods),
    n_methods, gsea_fdr, ora_fdr.
    """
    records: list[dict] = []

    # Collect significant terms: normalized_name -> (display_name, score)
    gw_map: dict[str, str] = {}  # normalized -> display name
    if not gw_results.empty and "go_name" in gw_results.columns and gw_padj_col in gw_results.columns:
        gw_sig = gw_results[gw_results[gw_padj_col] <= gw_padj_threshold]
        for name in gw_sig["go_name"].dropna().unique():
            gw_map[_normalize_term(name)] = name

    gsea_map: dict[str, tuple[str, float]] = {}  # normalized -> (display, fdr)
    if not gsea_results.empty and "term" in gsea_results.columns and "fdr" in gsea_results.columns:
        gsea_sig = gsea_results[gsea_results["fdr"] <= gsea_fdr_threshold]
        for _, row in gsea_sig.iterrows():
            norm = _normalize_term(row["term"])
            if norm not in gsea_map or row["fdr"] < gsea_map[norm][1]:
                gsea_map[norm] = (row["term"], row["fdr"])

    ora_map: dict[str, tuple[str, float]] = {}
    if not ora_results.empty and "term" in ora_results.columns and "fdr" in ora_results.columns:
        ora_sig = ora_results[ora_results["fdr"] <= ora_fdr_threshold]
        for _, row in ora_sig.iterrows():
            norm = _normalize_term(row["term"])
            if norm not in ora_map or row["fdr"] < ora_map[norm][1]:
                ora_map[norm] = (row["term"], row["fdr"])

    # Match across all normalized terms
    all_norms = set(gw_map.keys()) | set(gsea_map.keys()) | set(ora_map.keys())
    for norm in sorted(all_norms):
        methods = []
        display_name = norm  # fallback
        if norm in gw_map:
            methods.append("GeneWalk")
            display_name = gw_map[norm]
        if norm in gsea_map:
            methods.append("GSEA")
            if display_name == norm:
                display_name = gsea_map[norm][0]
        if norm in ora_map:
            methods.append("ORA")
            if display_name == norm:
                display_name = ora_map[norm][0]

        if len(methods) >= 2:
            records.append({
                "term": display_name,
                "found_in": ", ".join(methods),
                "n_methods": len(methods),
                "gsea_fdr": gsea_map[norm][1] if norm in gsea_map else None,
                "ora_fdr": ora_map[norm][1] if norm in ora_map else None,
            })

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values("n_methods", ascending=False).reset_index(drop=True)
