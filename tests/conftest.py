"""Shared fixtures for GeneWalk App tests."""

import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Columns matching real GeneWalk output
RESULTS_COLUMNS = [
    "hgnc_symbol",
    "hgnc_id",
    "go_name",
    "go_id",
    "go_domain",
    "sim",
    "sem_sim",
    "cilow",
    "ciupp",
    "global_padj",
    "gene_padj",
]


@pytest.fixture()
def sample_results_df() -> pd.DataFrame:
    """DataFrame mimicking genewalk_results.csv with ~20 rows."""
    data = [
        # BRAF -- several GO terms, mix of significant and not
        ("BRAF", "HGNC:1097", "protein kinase activity", "GO:0004672",
         "molecular_function", 0.85, 0.02, 0.81, 0.89, 0.001, 0.002),
        ("BRAF", "HGNC:1097", "MAPK cascade", "GO:0000165",
         "biological_process", 0.92, 0.01, 0.90, 0.94, 0.0005, 0.001),
        ("BRAF", "HGNC:1097", "cytoplasm", "GO:0005737",
         "cellular_component", 0.45, 0.05, 0.35, 0.55, 0.6, 0.5),
        ("BRAF", "HGNC:1097", "protein phosphorylation", "GO:0006468",
         "biological_process", 0.78, 0.03, 0.72, 0.84, 0.01, 0.02),
        # MAP2K1 -- significant terms
        ("MAP2K1", "HGNC:6840", "MAPK cascade", "GO:0000165",
         "biological_process", 0.90, 0.01, 0.88, 0.92, 0.0003, 0.0005),
        ("MAP2K1", "HGNC:6840", "protein kinase activity", "GO:0004672",
         "molecular_function", 0.82, 0.02, 0.78, 0.86, 0.005, 0.008),
        ("MAP2K1", "HGNC:6840", "nucleus", "GO:0005634",
         "cellular_component", 0.38, 0.06, 0.26, 0.50, 0.8, 0.7),
        # MAPK1 -- mix
        ("MAPK1", "HGNC:6871", "MAPK cascade", "GO:0000165",
         "biological_process", 0.88, 0.02, 0.84, 0.92, 0.002, 0.003),
        ("MAPK1", "HGNC:6871", "protein binding", "GO:0005515",
         "molecular_function", 0.55, 0.04, 0.47, 0.63, 0.3, 0.25),
        ("MAPK1", "HGNC:6871", "cell proliferation", "GO:0008283",
         "biological_process", 0.72, 0.03, 0.66, 0.78, 0.02, 0.03),
        # KRAS -- mostly significant
        ("KRAS", "HGNC:6407", "GTPase activity", "GO:0003924",
         "molecular_function", 0.91, 0.01, 0.89, 0.93, 0.0001, 0.0002),
        ("KRAS", "HGNC:6407", "signal transduction", "GO:0007165",
         "biological_process", 0.84, 0.02, 0.80, 0.88, 0.003, 0.005),
        ("KRAS", "HGNC:6407", "plasma membrane", "GO:0005886",
         "cellular_component", 0.68, 0.03, 0.62, 0.74, 0.04, 0.06),
        # RAF1 -- not significant
        ("RAF1", "HGNC:9829", "protein kinase activity", "GO:0004672",
         "molecular_function", 0.35, 0.08, 0.19, 0.51, 0.5, 0.4),
        ("RAF1", "HGNC:9829", "MAPK cascade", "GO:0000165",
         "biological_process", 0.40, 0.07, 0.26, 0.54, 0.45, 0.35),
        # MYC -- one significant
        ("MYC", "HGNC:7553", "transcription factor activity", "GO:0003700",
         "molecular_function", 0.80, 0.02, 0.76, 0.84, 0.008, 0.01),
        ("MYC", "HGNC:7553", "cell cycle", "GO:0007049",
         "biological_process", 0.42, 0.06, 0.30, 0.54, 0.55, 0.45),
        # FOS
        ("FOS", "HGNC:3796", "transcription factor activity", "GO:0003700",
         "molecular_function", 0.76, 0.03, 0.70, 0.82, 0.012, 0.015),
        ("FOS", "HGNC:3796", "AP-1 complex", "GO:0035976",
         "cellular_component", 0.65, 0.04, 0.57, 0.73, 0.08, 0.1),
        # NRAS
        ("NRAS", "HGNC:7989", "GTPase activity", "GO:0003924",
         "molecular_function", 0.87, 0.02, 0.83, 0.91, 0.002, 0.003),
    ]
    return pd.DataFrame(data, columns=RESULTS_COLUMNS)


@pytest.fixture()
def sample_gene_list() -> list[str]:
    """List of 10 gene symbols."""
    return ["BRAF", "MAP2K1", "MAPK1", "KRAS", "RAF1", "MYC", "FOS", "NRAS", "MAP2K2", "MAPK3"]


@pytest.fixture()
def empty_results_df() -> pd.DataFrame:
    """Empty DataFrame with correct GeneWalk results columns."""
    return pd.DataFrame(columns=RESULTS_COLUMNS)


@pytest.fixture()
def tmp_results_csv(sample_results_df: pd.DataFrame, tmp_path: Path) -> Path:
    """Write sample_results_df to a temporary CSV file and return the path."""
    csv_path = tmp_path / "genewalk_results.csv"
    sample_results_df.to_csv(csv_path, index=False)
    return csv_path
