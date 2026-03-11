"""Tests for genewalk_app.runner — backend functions for running and parsing GeneWalk."""

import csv
from pathlib import Path

import pandas as pd
import pytest

from genewalk_app.runner import (
    _parse_progress,
    _sanitize_project_name,
    filter_results,
    find_results_csv,
    get_gene_summary,
    load_results,
    save_gene_list,
)


# ---------------------------------------------------------------------------
# save_gene_list
# ---------------------------------------------------------------------------

class TestSaveGeneList:
    def test_save_gene_list(self, sample_gene_list: list[str], tmp_path: Path):
        """Verify gene list is written correctly to a file."""
        dest = tmp_path / "genes.txt"
        result = save_gene_list(sample_gene_list, dest)
        assert result == dest
        assert dest.exists()
        lines = dest.read_text(encoding="utf-8").strip().splitlines()
        assert lines == sample_gene_list

    def test_save_gene_list_strips_whitespace(self, tmp_path: Path):
        """Whitespace-padded genes should be cleaned."""
        genes = ["  BRAF  ", " KRAS\t", "  ", "MAPK1"]
        dest = tmp_path / "genes.txt"
        save_gene_list(genes, dest)
        lines = dest.read_text(encoding="utf-8").strip().splitlines()
        assert lines == ["BRAF", "KRAS", "MAPK1"]

    def test_save_gene_list_empty_raises(self, tmp_path: Path):
        """An empty gene list should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            save_gene_list([], tmp_path / "genes.txt")

    def test_save_gene_list_all_blank_raises(self, tmp_path: Path):
        """A list of only blank strings should raise ValueError."""
        with pytest.raises(ValueError, match="empty"):
            save_gene_list(["", "  ", "\t"], tmp_path / "genes.txt")

    def test_save_gene_list_creates_parent_dirs(self, tmp_path: Path):
        """Parent directories should be created if they don't exist."""
        dest = tmp_path / "nested" / "deep" / "genes.txt"
        save_gene_list(["BRAF"], dest)
        assert dest.exists()


# ---------------------------------------------------------------------------
# load_results / find_results_csv
# ---------------------------------------------------------------------------

class TestLoadResults:
    def test_parse_results_csv(self, tmp_results_csv: Path):
        """Loading a valid results CSV should return a DataFrame with correct shape."""
        df = load_results(tmp_results_csv)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 20
        assert "hgnc_symbol" in df.columns
        assert "go_name" in df.columns
        assert "sim" in df.columns

    def test_empty_csv_raises(self, tmp_path: Path):
        """An empty CSV should raise ValueError."""
        empty_csv = tmp_path / "empty.csv"
        empty_csv.write_text("", encoding="utf-8")
        with pytest.raises(ValueError):
            load_results(empty_csv)

    def test_invalid_results_format(self, tmp_path: Path):
        """A CSV missing required columns should raise ValueError."""
        bad_csv = tmp_path / "bad.csv"
        bad_csv.write_text("col_a,col_b\n1,2\n3,4\n", encoding="utf-8")
        with pytest.raises(ValueError, match="missing required columns"):
            load_results(bad_csv)

    def test_malformed_csv_raises(self, tmp_path: Path):
        """A completely malformed file should raise ValueError."""
        bad_file = tmp_path / "garbage.csv"
        bad_file.write_bytes(b"\x00\x01\x02\x03")
        with pytest.raises((ValueError, Exception)):
            load_results(bad_file)

    def test_find_results_csv_direct(self, tmp_path: Path):
        """find_results_csv should locate genewalk_results.csv at top level."""
        csv_path = tmp_path / "genewalk_results.csv"
        csv_path.write_text("hgnc_symbol,go_name,sim\nBRAF,kinase,0.8\n")
        result = find_results_csv(tmp_path)
        assert result == csv_path

    def test_find_results_csv_nested(self, tmp_path: Path):
        """find_results_csv should find the file in subdirectories."""
        nested = tmp_path / "subdir"
        nested.mkdir()
        csv_path = nested / "genewalk_results.csv"
        csv_path.write_text("hgnc_symbol,go_name,sim\nBRAF,kinase,0.8\n")
        result = find_results_csv(tmp_path)
        assert result is not None
        assert result.name == "genewalk_results.csv"

    def test_find_results_csv_missing(self, tmp_path: Path):
        """find_results_csv should return None when no CSV exists."""
        assert find_results_csv(tmp_path) is None


# ---------------------------------------------------------------------------
# filter_results
# ---------------------------------------------------------------------------

class TestFilterResults:
    def test_filter_results_by_padj(self, sample_results_df: pd.DataFrame):
        """Filtering by gene_padj <= 0.05 should keep only significant rows."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert len(filtered) > 0
        assert all(filtered["gene_padj"] <= 0.05)
        # Some rows should have been removed
        assert len(filtered) < len(sample_results_df)

    def test_filter_results_by_global_padj(self, sample_results_df: pd.DataFrame):
        """Filtering by global_padj should use that column."""
        filtered = filter_results(sample_results_df, padj_col="global_padj", padj_threshold=0.01)
        assert all(filtered["global_padj"] <= 0.01)

    def test_filter_results_by_domain(self, sample_results_df: pd.DataFrame):
        """Filtering by GO domain should keep only matching rows."""
        filtered = filter_results(
            sample_results_df,
            padj_col="gene_padj",
            padj_threshold=1.0,  # keep all by p-value
            go_domain="biological_process",
        )
        assert all(filtered["go_domain"] == "biological_process")

    def test_filter_results_combined(self, sample_results_df: pd.DataFrame):
        """Filtering by both p-value and domain simultaneously."""
        filtered = filter_results(
            sample_results_df,
            padj_col="gene_padj",
            padj_threshold=0.05,
            go_domain="molecular_function",
        )
        assert all(filtered["gene_padj"] <= 0.05)
        assert all(filtered["go_domain"] == "molecular_function")

    def test_filter_results_no_padj_col(self, sample_results_df: pd.DataFrame):
        """Filtering with a nonexistent padj column should return all rows."""
        filtered = filter_results(
            sample_results_df, padj_col="nonexistent_col", padj_threshold=0.05,
        )
        assert len(filtered) == len(sample_results_df)

    def test_filter_empty_results(self, empty_results_df: pd.DataFrame):
        """Filtering an empty DataFrame should return an empty DataFrame."""
        filtered = filter_results(empty_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert isinstance(filtered, pd.DataFrame)
        assert len(filtered) == 0

    def test_filter_strict_threshold(self, sample_results_df: pd.DataFrame):
        """A very strict threshold should return very few rows."""
        filtered = filter_results(
            sample_results_df, padj_col="gene_padj", padj_threshold=0.0001,
        )
        # Only the single KRAS GTPase entry has gene_padj=0.0002; nothing <= 0.0001
        assert len(filtered) == 0


# ---------------------------------------------------------------------------
# get_gene_summary
# ---------------------------------------------------------------------------

class TestGetGeneSummary:
    def test_summarize_gene_results(self, sample_results_df: pd.DataFrame):
        """Per-gene summary should have one row per gene with correct columns."""
        summary = get_gene_summary(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert not summary.empty
        expected_cols = {"hgnc_symbol", "significant_go_terms", "top_go_term", "best_padj", "mean_similarity"}
        assert expected_cols.issubset(set(summary.columns))
        # Each row should be a unique gene
        assert summary["hgnc_symbol"].is_unique

    def test_summary_counts_correct(self, sample_results_df: pd.DataFrame):
        """Verify the significant GO term counts are accurate."""
        summary = get_gene_summary(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        # BRAF has gene_padj values: 0.002, 0.001, 0.5, 0.02 => 3 significant at 0.05
        braf_row = summary[summary["hgnc_symbol"] == "BRAF"]
        if not braf_row.empty:
            assert braf_row.iloc[0]["significant_go_terms"] == 3

    def test_summary_empty_when_nothing_significant(self, sample_results_df: pd.DataFrame):
        """A very strict threshold should produce an empty summary."""
        summary = get_gene_summary(sample_results_df, padj_col="gene_padj", padj_threshold=0.00001)
        assert summary.empty

    def test_summary_missing_padj_col(self, sample_results_df: pd.DataFrame):
        """Missing padj column should return an empty DataFrame."""
        summary = get_gene_summary(sample_results_df, padj_col="nonexistent")
        assert summary.empty

    def test_empty_results_handling(self, empty_results_df: pd.DataFrame):
        """Empty DataFrame should produce an empty summary."""
        summary = get_gene_summary(empty_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert summary.empty


# ---------------------------------------------------------------------------
# _sanitize_project_name
# ---------------------------------------------------------------------------

class TestSanitizeProjectName:
    def test_normal_name(self):
        assert _sanitize_project_name("my_project") == "my_project"

    def test_strips_path_separators(self):
        assert "/" not in _sanitize_project_name("../../etc/passwd")
        assert "\\" not in _sanitize_project_name("..\\..\\windows\\system32")

    def test_replaces_special_chars(self):
        result = _sanitize_project_name('my:project*name?"<>|')
        assert all(c not in result for c in ':*?"<>|')

    def test_empty_string_fallback(self):
        assert _sanitize_project_name("") == "genewalk_analysis"

    def test_whitespace_only_fallback(self):
        assert _sanitize_project_name("   ") == "genewalk_analysis"


# ---------------------------------------------------------------------------
# _parse_progress
# ---------------------------------------------------------------------------

class TestParseProgress:
    def test_downloading_pattern(self):
        assert _parse_progress("Downloading GO OBO file...") == "Downloading resources..."

    def test_deepwalk_graph_pattern(self):
        result = _parse_progress("DeepWalk on gene-GO graph...")
        assert result == "Running DeepWalk on gene network..."

    def test_no_match_returns_none(self):
        assert _parse_progress("some random log line") is None

    def test_statistics_pattern(self):
        result = _parse_progress("Performing statistical tests...")
        assert result == "Computing statistics..."
