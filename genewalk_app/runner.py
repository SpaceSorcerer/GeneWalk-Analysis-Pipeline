"""Backend module for running GeneWalk and parsing results."""

import subprocess
import tempfile
from pathlib import Path

import pandas as pd


DEFAULT_BASE_DIR = Path(tempfile.gettempdir()) / "genewalk_runs"


def save_gene_list(genes: list[str], dest: Path) -> Path:
    """Write a list of gene identifiers to a text file."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(g.strip() for g in genes if g.strip()) + "\n")
    return dest


def run_genewalk(
    gene_file: Path,
    project: str,
    id_type: str = "hgnc_symbol",
    nproc: int = 4,
    nreps_graph: int = 3,
    nreps_null: int = 3,
    alpha_fdr: float = 1.0,
    base_folder: Path | None = None,
) -> dict:
    """Run GeneWalk CLI and return paths to output files.

    Returns a dict with keys: 'return_code', 'stdout', 'stderr', 'output_dir'.
    """
    base = base_folder or DEFAULT_BASE_DIR
    base.mkdir(parents=True, exist_ok=True)

    cmd = [
        "genewalk",
        "--project", project,
        "--genes", str(gene_file),
        "--id_type", id_type,
        "--nproc", str(nproc),
        "--nreps_graph", str(nreps_graph),
        "--nreps_null", str(nreps_null),
        "--alpha_fdr", str(alpha_fdr),
        "--base_folder", str(base),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=14400,  # 4 hour timeout
    )

    output_dir = base / project

    return {
        "return_code": result.returncode,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "output_dir": output_dir,
    }


def find_results_csv(output_dir: Path) -> Path | None:
    """Locate the genewalk_results.csv file in the output directory."""
    candidate = output_dir / "genewalk_results.csv"
    if candidate.exists():
        return candidate
    # Search recursively as a fallback
    for p in output_dir.rglob("genewalk_results.csv"):
        return p
    return None


def load_results(csv_path: Path) -> pd.DataFrame:
    """Load and clean GeneWalk results CSV."""
    df = pd.read_csv(csv_path)
    return df


def filter_results(
    df: pd.DataFrame,
    padj_col: str = "gene_padj",
    padj_threshold: float = 0.05,
    go_domain: str | None = None,
) -> pd.DataFrame:
    """Filter results by adjusted p-value and optional GO domain."""
    filtered = df.copy()
    if padj_col in filtered.columns:
        filtered = filtered[filtered[padj_col] <= padj_threshold]
    if go_domain and "go_domain" in filtered.columns:
        filtered = filtered[filtered["go_domain"] == go_domain]
    return filtered.sort_values(padj_col, ascending=True) if padj_col in filtered.columns else filtered


def get_gene_summary(df: pd.DataFrame, padj_col: str = "gene_padj", padj_threshold: float = 0.05) -> pd.DataFrame:
    """Summarize results per gene: count of significant GO terms."""
    if padj_col not in df.columns:
        return pd.DataFrame()
    sig = df[df[padj_col] <= padj_threshold]
    summary = (
        sig.groupby("hgnc_symbol")
        .agg(
            significant_go_terms=("go_name", "count"),
            top_go_term=("go_name", "first"),
            best_padj=(padj_col, "min"),
            mean_similarity=("sim", "mean"),
        )
        .sort_values("significant_go_terms", ascending=False)
        .reset_index()
    )
    return summary
