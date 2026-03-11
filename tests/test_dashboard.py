"""Tests for genewalk_app.dashboard — pure logic functions.

These tests exercise the filtering and metric computation logic used by the
dashboard without importing Streamlit.  The actual Streamlit rendering
(render_dashboard) requires a running Streamlit session and is not tested here.
"""

import io

import pandas as pd
import pytest

# Import the pure-logic functions from runner (used by dashboard)
from genewalk_app.runner import filter_results, get_gene_summary


# ---------------------------------------------------------------------------
# Filter by p-value column
# ---------------------------------------------------------------------------

class TestFilterByPadjColumn:
    def test_filter_by_gene_padj(self, sample_results_df: pd.DataFrame):
        """Filtering by gene_padj should use that column."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert all(filtered["gene_padj"] <= 0.05)

    def test_filter_by_global_padj(self, sample_results_df: pd.DataFrame):
        """Filtering by global_padj should use that column."""
        filtered = filter_results(sample_results_df, padj_col="global_padj", padj_threshold=0.05)
        assert all(filtered["global_padj"] <= 0.05)

    def test_filter_nonexistent_padj_col(self, sample_results_df: pd.DataFrame):
        """A nonexistent padj column should return all rows unfiltered."""
        filtered = filter_results(sample_results_df, padj_col="nonexistent", padj_threshold=0.05)
        assert len(filtered) == len(sample_results_df)


# ---------------------------------------------------------------------------
# Filter by threshold
# ---------------------------------------------------------------------------

class TestFilterByThreshold:
    def test_threshold_0_05(self, sample_results_df: pd.DataFrame):
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        count_05 = len(filtered)
        assert count_05 > 0

    def test_threshold_0_01(self, sample_results_df: pd.DataFrame):
        filtered_01 = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.01)
        filtered_05 = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        # Stricter threshold should yield fewer or equal rows
        assert len(filtered_01) <= len(filtered_05)

    def test_threshold_1_0_keeps_all(self, sample_results_df: pd.DataFrame):
        """Threshold of 1.0 should keep all rows."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=1.0)
        assert len(filtered) == len(sample_results_df)


# ---------------------------------------------------------------------------
# Filter by GO domain
# ---------------------------------------------------------------------------

class TestFilterByGoDomain:
    def test_biological_process(self, sample_results_df: pd.DataFrame):
        filtered = filter_results(
            sample_results_df, padj_col="gene_padj", padj_threshold=1.0,
            go_domain="biological_process",
        )
        assert len(filtered) > 0
        assert all(filtered["go_domain"] == "biological_process")

    def test_molecular_function(self, sample_results_df: pd.DataFrame):
        filtered = filter_results(
            sample_results_df, padj_col="gene_padj", padj_threshold=1.0,
            go_domain="molecular_function",
        )
        assert len(filtered) > 0
        assert all(filtered["go_domain"] == "molecular_function")

    def test_cellular_component(self, sample_results_df: pd.DataFrame):
        filtered = filter_results(
            sample_results_df, padj_col="gene_padj", padj_threshold=1.0,
            go_domain="cellular_component",
        )
        assert len(filtered) > 0
        assert all(filtered["go_domain"] == "cellular_component")

    def test_none_domain_keeps_all(self, sample_results_df: pd.DataFrame):
        """go_domain=None should not filter by domain."""
        filtered = filter_results(
            sample_results_df, padj_col="gene_padj", padj_threshold=1.0,
            go_domain=None,
        )
        assert len(filtered) == len(sample_results_df)


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

class TestMetricComputation:
    def test_total_pairs(self, sample_results_df: pd.DataFrame):
        """Total pairs = number of rows in unfiltered DataFrame."""
        assert len(sample_results_df) == 20

    def test_significant_count(self, sample_results_df: pd.DataFrame):
        """Significant count should be rows passing the filter."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        assert len(filtered) > 0
        assert len(filtered) < len(sample_results_df)

    def test_unique_genes(self, sample_results_df: pd.DataFrame):
        """Unique gene count in filtered results."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        unique_genes = filtered["hgnc_symbol"].nunique()
        assert unique_genes > 0
        # We know from fixture: BRAF(3 sig), MAP2K1(2 sig), MAPK1(2 sig),
        # KRAS(3 sig), MYC(1 sig), FOS(1 sig), NRAS(1 sig) = 7 genes with
        # at least one significant term. RAF1 has padj 0.4 and 0.35 => both sig at 0.05
        # Actually let me count: RAF1 gene_padj values are 0.4, 0.35 => both <= 0.05? No!
        # 0.4 > 0.05 and 0.35 > 0.05, so RAF1 has 0 significant.
        # So 7 unique genes.
        assert unique_genes >= 5

    def test_unique_go_terms(self, sample_results_df: pd.DataFrame):
        """Unique GO term count in filtered results."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        unique_terms = filtered["go_name"].nunique()
        assert unique_terms > 0


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

class TestCsvExport:
    def test_csv_export(self, sample_results_df: pd.DataFrame):
        """Verify filtered data can be exported to CSV and re-read."""
        filtered = filter_results(sample_results_df, padj_col="gene_padj", padj_threshold=0.05)
        buf = io.StringIO()
        filtered.to_csv(buf, index=False)

        # Re-read and verify
        buf.seek(0)
        reloaded = pd.read_csv(buf)
        assert len(reloaded) == len(filtered)
        assert list(reloaded.columns) == list(filtered.columns)

    def test_csv_export_empty(self, empty_results_df: pd.DataFrame):
        """Exporting an empty DataFrame should produce a valid CSV with headers."""
        buf = io.StringIO()
        empty_results_df.to_csv(buf, index=False)
        buf.seek(0)
        reloaded = pd.read_csv(buf)
        assert len(reloaded) == 0
        assert list(reloaded.columns) == list(empty_results_df.columns)
