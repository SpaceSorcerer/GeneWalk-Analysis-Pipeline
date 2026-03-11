"""Tests for genewalk_app.visualizations — Plotly chart functions."""

import pandas as pd
import plotly.graph_objects as go
import pytest

from genewalk_app.visualizations import (
    gene_bar_chart,
    gene_go_network,
    gene_similarity_heatmap,
    go_domain_pie,
    pvalue_distribution,
    summary_bar,
    volcano_plot,
)


def _is_figure(obj) -> bool:
    """Check if obj is a Plotly Figure."""
    return isinstance(obj, go.Figure)


# ---------------------------------------------------------------------------
# volcano_plot
# ---------------------------------------------------------------------------

class TestVolcanoPlot:
    def test_volcano_plot_renders(self, sample_results_df: pd.DataFrame):
        fig = volcano_plot(sample_results_df)
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_volcano_plot_custom_padj(self, sample_results_df: pd.DataFrame):
        fig = volcano_plot(sample_results_df, padj_col="global_padj", padj_threshold=0.01)
        assert _is_figure(fig)

    def test_volcano_plot_empty_data(self, empty_results_df: pd.DataFrame):
        fig = volcano_plot(empty_results_df)
        assert _is_figure(fig)

    def test_volcano_plot_missing_columns(self):
        df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
        fig = volcano_plot(df)
        assert _is_figure(fig)
        # Should have a title indicating missing columns
        assert "Missing" in fig.layout.title.text or "missing" in fig.layout.title.text.lower()


# ---------------------------------------------------------------------------
# gene_bar_chart
# ---------------------------------------------------------------------------

class TestGeneBarChart:
    def test_gene_bar_chart_renders(self, sample_results_df: pd.DataFrame):
        fig = gene_bar_chart(sample_results_df, "BRAF")
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_gene_bar_chart_nonexistent_gene(self, sample_results_df: pd.DataFrame):
        fig = gene_bar_chart(sample_results_df, "NONEXISTENT_GENE")
        assert _is_figure(fig)
        # Title should mention no results
        assert "No results" in fig.layout.title.text

    def test_gene_bar_chart_top_n(self, sample_results_df: pd.DataFrame):
        fig = gene_bar_chart(sample_results_df, "BRAF", top_n=2)
        assert _is_figure(fig)

    def test_gene_bar_chart_empty_data(self, empty_results_df: pd.DataFrame):
        fig = gene_bar_chart(empty_results_df, "BRAF")
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# gene_go_network
# ---------------------------------------------------------------------------

class TestNetworkGraph:
    def test_network_graph_renders(self, sample_results_df: pd.DataFrame):
        fig = gene_go_network(sample_results_df)
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_network_graph_with_gene_filter(self, sample_results_df: pd.DataFrame):
        fig = gene_go_network(sample_results_df, genes=["BRAF", "KRAS"])
        assert _is_figure(fig)

    def test_network_graph_max_edges(self, sample_results_df: pd.DataFrame):
        fig = gene_go_network(sample_results_df, max_edges=5)
        assert _is_figure(fig)

    def test_network_graph_empty_data(self, empty_results_df: pd.DataFrame):
        fig = gene_go_network(empty_results_df)
        assert _is_figure(fig)

    def test_network_graph_strict_threshold(self, sample_results_df: pd.DataFrame):
        """Very strict threshold may yield no significant results."""
        fig = gene_go_network(sample_results_df, padj_threshold=0.00001)
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# gene_similarity_heatmap
# ---------------------------------------------------------------------------

class TestHeatmap:
    def test_heatmap_renders(self, sample_results_df: pd.DataFrame):
        fig = gene_similarity_heatmap(sample_results_df)
        assert _is_figure(fig)

    def test_heatmap_with_gene_filter(self, sample_results_df: pd.DataFrame):
        fig = gene_similarity_heatmap(sample_results_df, genes=["BRAF", "KRAS"])
        assert _is_figure(fig)

    def test_heatmap_max_terms(self, sample_results_df: pd.DataFrame):
        fig = gene_similarity_heatmap(sample_results_df, max_terms=3)
        assert _is_figure(fig)

    def test_heatmap_empty_data(self, empty_results_df: pd.DataFrame):
        fig = gene_similarity_heatmap(empty_results_df)
        assert _is_figure(fig)

    def test_heatmap_strict_threshold(self, sample_results_df: pd.DataFrame):
        """Very strict threshold may yield no significant results."""
        fig = gene_similarity_heatmap(sample_results_df, padj_threshold=0.00001)
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# go_domain_pie
# ---------------------------------------------------------------------------

class TestDomainPie:
    def test_domain_pie_chart_renders(self, sample_results_df: pd.DataFrame):
        fig = go_domain_pie(sample_results_df)
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_domain_pie_empty_data(self, empty_results_df: pd.DataFrame):
        fig = go_domain_pie(empty_results_df)
        assert _is_figure(fig)

    def test_domain_pie_no_domain_column(self):
        df = pd.DataFrame({"hgnc_symbol": ["BRAF"], "sim": [0.8]})
        fig = go_domain_pie(df)
        assert _is_figure(fig)

    def test_domain_pie_strict_threshold(self, sample_results_df: pd.DataFrame):
        fig = go_domain_pie(sample_results_df, padj_threshold=0.00001)
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# pvalue_distribution
# ---------------------------------------------------------------------------

class TestPvalueHistogram:
    def test_pvalue_histogram_renders(self, sample_results_df: pd.DataFrame):
        fig = pvalue_distribution(sample_results_df)
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_pvalue_histogram_global_padj(self, sample_results_df: pd.DataFrame):
        fig = pvalue_distribution(sample_results_df, padj_col="global_padj")
        assert _is_figure(fig)

    def test_pvalue_histogram_missing_col(self):
        df = pd.DataFrame({"hgnc_symbol": ["BRAF"], "sim": [0.8]})
        fig = pvalue_distribution(df, padj_col="gene_padj")
        assert _is_figure(fig)

    def test_pvalue_histogram_empty_data(self, empty_results_df: pd.DataFrame):
        fig = pvalue_distribution(empty_results_df)
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# summary_bar
# ---------------------------------------------------------------------------

class TestSummaryBar:
    def test_gene_summary_bar_renders(self, sample_results_df: pd.DataFrame):
        from genewalk_app.runner import get_gene_summary
        summary = get_gene_summary(sample_results_df)
        fig = summary_bar(summary)
        assert _is_figure(fig)
        assert len(fig.data) > 0

    def test_summary_bar_empty(self):
        fig = summary_bar(pd.DataFrame())
        assert _is_figure(fig)

    def test_summary_bar_missing_columns(self):
        """DataFrame without expected columns should return gracefully."""
        df = pd.DataFrame({"x": [1, 2], "y": [3, 4]})
        fig = summary_bar(df)
        assert _is_figure(fig)


# ---------------------------------------------------------------------------
# Edge cases: single gene / single GO term
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_single_gene_data(self):
        """All charts should handle a DataFrame with only one gene."""
        df = pd.DataFrame({
            "hgnc_symbol": ["BRAF", "BRAF"],
            "hgnc_id": ["HGNC:1097", "HGNC:1097"],
            "go_name": ["kinase activity", "MAPK cascade"],
            "go_id": ["GO:0004672", "GO:0000165"],
            "go_domain": ["molecular_function", "biological_process"],
            "sim": [0.85, 0.92],
            "sem_sim": [0.02, 0.01],
            "cilow": [0.81, 0.90],
            "ciupp": [0.89, 0.94],
            "global_padj": [0.001, 0.0005],
            "gene_padj": [0.002, 0.001],
        })
        assert _is_figure(volcano_plot(df))
        assert _is_figure(gene_bar_chart(df, "BRAF"))
        assert _is_figure(gene_go_network(df))
        assert _is_figure(gene_similarity_heatmap(df))
        assert _is_figure(go_domain_pie(df))
        assert _is_figure(pvalue_distribution(df))

    def test_single_goterm_data(self):
        """All charts should handle a DataFrame with only one GO term."""
        df = pd.DataFrame({
            "hgnc_symbol": ["BRAF", "KRAS"],
            "hgnc_id": ["HGNC:1097", "HGNC:6407"],
            "go_name": ["kinase activity", "kinase activity"],
            "go_id": ["GO:0004672", "GO:0004672"],
            "go_domain": ["molecular_function", "molecular_function"],
            "sim": [0.85, 0.78],
            "sem_sim": [0.02, 0.03],
            "cilow": [0.81, 0.72],
            "ciupp": [0.89, 0.84],
            "global_padj": [0.001, 0.005],
            "gene_padj": [0.002, 0.008],
        })
        assert _is_figure(volcano_plot(df))
        assert _is_figure(gene_bar_chart(df, "BRAF"))
        assert _is_figure(gene_go_network(df))
        assert _is_figure(gene_similarity_heatmap(df))
        assert _is_figure(go_domain_pie(df))
        assert _is_figure(pvalue_distribution(df))

    def test_empty_data_all_charts(self, empty_results_df: pd.DataFrame):
        """All chart functions should handle empty DataFrames gracefully."""
        assert _is_figure(volcano_plot(empty_results_df))
        assert _is_figure(gene_bar_chart(empty_results_df, "BRAF"))
        assert _is_figure(gene_go_network(empty_results_df))
        assert _is_figure(gene_similarity_heatmap(empty_results_df))
        assert _is_figure(go_domain_pie(empty_results_df))
        assert _is_figure(pvalue_distribution(empty_results_df))
        assert _is_figure(summary_bar(empty_results_df))
