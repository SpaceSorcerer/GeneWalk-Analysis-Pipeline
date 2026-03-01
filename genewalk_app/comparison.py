"""DEG table parsing, gene list splitting, and cross-analysis comparison."""

from __future__ import annotations

import pandas as pd


def detect_deg_columns(df: pd.DataFrame) -> dict[str, str | None]:
    """Auto-detect gene, log2FC, and p-value columns in a DEG table.

    Returns a dict with keys 'gene', 'log2fc', 'padj' mapped to the best
    matching column name, or *None* if no match is found.
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
        "pvalue", "p_value", "pval", "p-value",
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
    out = out.dropna(subset=["gene", "log2fc", "padj"])
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
    sig = deg[deg["padj"] <= padj_threshold]
    up = sig[sig["log2fc"] >= fc_threshold]["gene"].unique().tolist()
    down = sig[sig["log2fc"] <= -fc_threshold]["gene"].unique().tolist()
    return up, down


def make_ranked_list(deg: pd.DataFrame) -> pd.DataFrame:
    """Create a ranked gene list for GSEA prerank from a DEG table.

    Returns a two-column DataFrame (gene, score) sorted by score descending.
    Genes are deduplicated by keeping the entry with the largest absolute score.
    """
    rnk = deg[["gene", "log2fc"]].copy()
    rnk = rnk.dropna()
    # Deduplicate: keep the row with the largest absolute log2fc per gene
    rnk["abs_score"] = rnk["log2fc"].abs()
    rnk = rnk.sort_values("abs_score", ascending=False).drop_duplicates(
        subset="gene", keep="first"
    )
    rnk = rnk.drop(columns="abs_score").sort_values("log2fc", ascending=False)
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

    Matches terms by substring overlap since GeneWalk uses GO term names
    while GSEA/ORA use pathway database names.

    Returns a DataFrame with columns: term, found_in (list of methods),
    gw_padj, gsea_fdr, ora_fdr.
    """
    records: list[dict] = []

    # Collect significant terms from each method
    gw_terms: set[str] = set()
    if not gw_results.empty and "go_name" in gw_results.columns and gw_padj_col in gw_results.columns:
        gw_sig = gw_results[gw_results[gw_padj_col] <= gw_padj_threshold]
        gw_terms = set(gw_sig["go_name"].dropna().unique())

    gsea_terms: dict[str, float] = {}
    if not gsea_results.empty and "term" in gsea_results.columns and "fdr" in gsea_results.columns:
        gsea_sig = gsea_results[gsea_results["fdr"] <= gsea_fdr_threshold]
        for _, row in gsea_sig.iterrows():
            gsea_terms[row["term"]] = row["fdr"]

    ora_terms: dict[str, float] = {}
    if not ora_results.empty and "term" in ora_results.columns and "fdr" in ora_results.columns:
        ora_sig = ora_results[ora_results["fdr"] <= ora_fdr_threshold]
        for _, row in ora_sig.iterrows():
            ora_terms[row["term"]] = row["fdr"]

    # Combine all terms
    all_terms = gw_terms | set(gsea_terms.keys()) | set(ora_terms.keys())
    for term in sorted(all_terms):
        methods = []
        if term in gw_terms:
            methods.append("GeneWalk")
        if term in gsea_terms:
            methods.append("GSEA")
        if term in ora_terms:
            methods.append("ORA")
        if len(methods) >= 2:
            records.append({
                "term": term,
                "found_in": ", ".join(methods),
                "n_methods": len(methods),
                "gsea_fdr": gsea_terms.get(term),
                "ora_fdr": ora_terms.get(term),
            })

    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values("n_methods", ascending=False).reset_index(drop=True)
