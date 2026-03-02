"""GSEA and Over-Representation Analysis runner using GSEApy."""

from __future__ import annotations

import tempfile
from collections.abc import Callable
from pathlib import Path

import pandas as pd

# Default gene set libraries offered in the UI.
DEFAULT_GENE_SETS: list[str] = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "WikiPathway_2023_Human",
]


def run_gsea_prerank(
    ranked_genes: pd.DataFrame,
    gene_sets: list[str],
    outdir: Path | None = None,
    permutation_num: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
    threads: int = 4,
    seed: int = 123,
    on_progress: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    """Run GSEA prerank on a ranked gene list.

    Parameters
    ----------
    ranked_genes : pd.DataFrame
        Two-column DataFrame: gene identifier and numeric score (e.g. log2FC).
        The first column is the gene name, the second is the ranking metric.
    gene_sets : list[str]
        Gene set library names (from Enrichr/MSigDB).
    outdir : Path, optional
        Directory to write GSEA output.  Uses a temp dir if *None*.
    permutation_num : int
        Number of permutations for significance testing.
    min_size, max_size : int
        Min/max gene set sizes to consider.
    threads : int
        Parallel workers.
    seed : int
        Random seed for reproducibility.
    on_progress : callable, optional
        Called with status messages.

    Returns
    -------
    pd.DataFrame
        Columns: Term, ES, NES, NOM p-val, FDR q-val, FWER p-val,
        Tag %, Gene %, Lead_genes, gene_set_library.
    """
    import gseapy

    out = str(outdir) if outdir else tempfile.mkdtemp(prefix="gsea_prerank_")

    all_results: list[pd.DataFrame] = []
    for i, gs in enumerate(gene_sets, 1):
        if on_progress:
            on_progress(f"GSEA prerank: {gs} ({i}/{len(gene_sets)})...")
        try:
            pre = gseapy.prerank(
                rnk=ranked_genes,
                gene_sets=gs,
                outdir=out,
                permutation_num=permutation_num,
                min_size=min_size,
                max_size=max_size,
                threads=threads,
                seed=seed,
                no_plot=True,
                verbose=False,
            )
            res = pre.res2d.copy()
            res["gene_set_library"] = gs
            all_results.append(res)
        except Exception as exc:
            if on_progress:
                on_progress(f"Warning: {gs} failed ({exc})")

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    return _clean_gsea_results(combined)


def run_ora(
    gene_list: list[str],
    gene_sets: list[str],
    background: list[str] | int | None = None,
    outdir: Path | None = None,
    on_progress: Callable[[str], None] | None = None,
) -> pd.DataFrame:
    """Run over-representation analysis (ORA) on a gene list.

    Parameters
    ----------
    gene_list : list[str]
        Gene identifiers to test for enrichment.
    gene_sets : list[str]
        Gene set library names.
    background : list[str] | int | None
        Background gene set or genome size.  *None* uses the library default.
    outdir : Path, optional
        Output directory.  Uses a temp dir if *None*.
    on_progress : callable, optional
        Called with status messages.

    Returns
    -------
    pd.DataFrame
        Columns: Term, Overlap, P-value, Adjusted P-value, Odds Ratio,
        Combined Score, Genes, gene_set_library.
    """
    import gseapy

    out = str(outdir) if outdir else tempfile.mkdtemp(prefix="ora_")

    all_results: list[pd.DataFrame] = []
    for i, gs in enumerate(gene_sets, 1):
        if on_progress:
            on_progress(f"ORA: {gs} ({i}/{len(gene_sets)})...")
        try:
            enr = gseapy.enrich(
                gene_list=gene_list,
                gene_sets=gs,
                background=background,
                outdir=out,
                no_plot=True,
                verbose=False,
            )
            res = enr.results.copy()
            res["gene_set_library"] = gs
            all_results.append(res)
        except Exception as exc:
            if on_progress:
                on_progress(f"Warning: {gs} failed ({exc})")

    if not all_results:
        return pd.DataFrame()

    combined = pd.concat(all_results, ignore_index=True)
    return _clean_ora_results(combined)


def _clean_gsea_results(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize GSEA prerank result column names.

    GSEApy >= 1.1 uses lowercase column names directly (term, es, nes,
    fdr, etc.).  Older versions used Title Case (Term, ES, NES, etc.).
    This handles both.
    """
    rename = {
        # Older GSEApy versions
        "Term": "term",
        "Name": "term",
        "ES": "es",
        "NES": "nes",
        "NOM p-val": "pval",
        "FDR q-val": "fdr",
        "FWER p-val": "fwerp",
        "Tag %": "tag_pct",
        "Gene %": "gene_pct",
        "Lead_genes": "lead_genes",
        # GSEApy >= 1.1 sometimes uses these
        "NOM p-value": "pval",
        "FDR q-value": "fdr",
        "FWER p-value": "fwerp",
    }
    # Only rename columns that exist and don't conflict with existing lowercase
    for old, new in rename.items():
        if old in df.columns and new not in df.columns:
            df = df.rename(columns={old: new})
    # Ensure key numeric columns exist
    for col in ("nes", "es", "fdr", "pval"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _clean_ora_results(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize ORA result column names and parse overlap ratios."""
    rename = {
        "Term": "term",
        "Overlap": "overlap",
        "P-value": "pval",
        "Adjusted P-value": "fdr",
        "Odds Ratio": "odds_ratio",
        "Combined Score": "combined_score",
        "Genes": "genes",
    }
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    for col in ("pval", "fdr", "odds_ratio", "combined_score"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Parse "3/50" overlap strings into numeric columns for visualization
    if "overlap" in df.columns:
        parts = df["overlap"].astype(str).str.split("/", expand=True)
        if parts.shape[1] == 2:
            df["overlap_count"] = pd.to_numeric(parts[0], errors="coerce")
            df["overlap_size"] = pd.to_numeric(parts[1], errors="coerce")
            df["overlap_ratio"] = df["overlap_count"] / df["overlap_size"]
    return df
