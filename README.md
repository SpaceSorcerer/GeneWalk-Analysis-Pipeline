# GeneWalk App

> Run GeneWalk and explore gene function results interactively.

## What It Does

[GeneWalk](https://github.com/churchmanlab/genewalk) uses network representation learning on gene-Gene Ontology (GO) networks to identify which biological functions are specifically relevant to your set of genes, rather than generically associated.

- **Context-specific gene function analysis** -- GeneWalk builds a representation of your gene list within the GO network and uses DeepWalk embeddings to score gene-function associations by similarity and statistical significance.
- **Interactive dashboards** -- This app wraps GeneWalk with a Streamlit UI featuring eight Plotly visualizations: volcano plot, per-gene explorer, network graph, similarity heatmap, GO domain distribution, p-value histogram, gene summary ranking, and a filterable data table.
- **Easy input, flexible output** -- Paste a gene list, upload a file, or load demo data. Filter results by p-value and GO domain. Export filtered tables to CSV.

## Installation

### Option 1: Download Executable (Recommended)

Download a pre-built executable from [Releases](https://github.com/SpaceSorcerer/GeneWalk-Analysis-Pipeline/releases). No Python installation required.

Available for Windows, macOS, and Linux.

### Option 2: pip install

```bash
# Install the app with desktop analysis support (includes GeneWalk):
pip install genewalk-app[desktop]

# Or install just the web viewer (no GeneWalk, view pre-computed results only):
pip install genewalk-app
```

Alternatively, install from source:

```bash
git clone https://github.com/SpaceSorcerer/GeneWalk-Analysis-Pipeline.git
cd GeneWalk-Analysis-Pipeline
pip install -r requirements-desktop.txt   # full analysis
# or
pip install -r requirements.txt           # web viewer only
```

### Option 3: Docker (visualization only)

```bash
docker build -t genewalk-app .
docker run -p 8501:8501 genewalk-app
```

The Docker image runs the web viewer at `http://localhost:8501`. Upload a `genewalk_results.csv` to visualize, or use the built-in demo data.

## Quick Start

### Run GeneWalk

1. Launch the app (`streamlit run app.py` or open the desktop executable)
2. Select **Run Analysis** in the sidebar
3. Paste or upload your gene list (one gene per line, HGNC symbols or other supported ID types)
4. Select species and configure parameters (CPU cores, graph reps, null reps)
5. Click **Run GeneWalk**
6. When the run completes, results load automatically into the dashboard

### View Existing Results

1. Select **View Results** in the sidebar
2. Upload a `genewalk_results.csv` from a previous GeneWalk run
3. Or click **Load Sample Data** to explore the built-in MAPK/ERK pathway demo

## Dashboard Features

- **Volcano Plot** -- Similarity score vs. statistical significance, with adjustable threshold line
- **Per-Gene Explorer** -- Select any gene to see its top GO term associations ranked by similarity
- **Network Graph** -- Interactive gene-GO term bipartite network with adjustable edge limits
- **Similarity Heatmap** -- Gene x GO term matrix with hierarchical clustering
- **GO Domain Distribution** -- Pie chart breakdown by Biological Process, Molecular Function, and Cellular Component
- **P-value Distribution** -- Histogram of adjusted p-values as a quality check
- **Gene Summary** -- Rank genes by number of significant GO term associations
- **Data Table** -- Filter, sort, and export the full results table to CSV

All visualizations respond to the global p-value threshold slider and GO domain filter in the sidebar.

## Configuration

### Gene ID Types

| ID Type | Example | Source |
|---------|---------|--------|
| `hgnc_symbol` | BRAF | HUGO Gene Nomenclature |
| `hgnc_id` | HGNC:1097 | HUGO Gene Nomenclature |
| `ensembl_id` | ENSG00000157764 | Ensembl |
| `mgi_id` | MGI:88190 | Mouse Genome Informatics |
| `rgd_id` | RGD:2261 | Rat Genome Database |
| `entrez_human` | 673 | NCBI Entrez Gene (human) |
| `entrez_mouse` | 109880 | NCBI Entrez Gene (mouse) |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| CPU cores | 4 | Parallel workers for DeepWalk. Use 2--4 for laptops, 8--16 for workstations. |
| Graph reps | 3 | DeepWalk repetitions on the gene-GO network. Use 3 for testing, 10--15 for publication. |
| Null reps | 3 | DeepWalk repetitions on null networks. Should match graph reps. |
| FDR threshold | 1.0 | Set to 1.0 to keep all results; set to 0.05 for pre-filtering. |

## Sample Data

The built-in demo includes 20 MAPK/ERK signaling pathway genes (BRAF, MAP2K1, MAPK1, KRAS, RAF1, etc.) with ~400 pre-computed gene-GO term associations. Click **Load Sample Data** in the app to explore it instantly.

To regenerate sample data from scratch:

```bash
python generate_sample_data.py
```

## Citation

If you use this software, please cite:

> Ietswaart R, Gyori BM, Bachman JA, Sorger PK, Bhatt DL, Churchman LS.
> GeneWalk identifies relevant gene functions for a biological context using
> network representation learning. *Genome Biology* 22, 55 (2021).
> https://doi.org/10.1186/s13059-021-02264-8

See [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

## License

MIT License. See [LICENSE](LICENSE) for details.
