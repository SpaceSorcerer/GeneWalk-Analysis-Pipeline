# GeneWalk App

A focused desktop and web app for running [GeneWalk](https://github.com/churchmanlab/genewalk) and interactively exploring gene-function associations.

## What It Does

GeneWalk identifies which Gene Ontology (GO) functions are specifically relevant to your biological context using network representation learning. This app provides:

- **Run GeneWalk** directly from the browser with configurable parameters
- **Upload previous results** (`genewalk_results.csv`) to skip the analysis step
- **Interactive visualizations** powered by Plotly:
  - Volcano plot (similarity vs. significance)
  - Per-gene GO term bar charts
  - Gene-GO bipartite network graph
  - GO domain distribution pie chart
  - Gene-GO term similarity heatmap
  - P-value distribution histogram
  - Gene summary ranking
- **Real-time filtering** by p-value threshold and GO domain
- **CSV export** of filtered results
- **Instant demo mode** with pre-computed MAPK/ERK pathway data

## Installation

### Option 1: pip install (recommended)

```bash
pip install -r requirements.txt

# To run GeneWalk analyses (not just view results):
pip install genewalk>=1.6
```

### Option 2: Download executable

Download a pre-built executable from [GitHub Releases](https://github.com/SpaceSorcerer/GeneWalk-Analysis-Pipeline/releases). No Python installation required.

### Option 3: Docker

```bash
docker build -t genewalk-app .
docker run -p 8501:8501 genewalk-app
```

## Quick Start

### Web app (view results + run analysis if GeneWalk installed)

```bash
streamlit run app.py
```

### Desktop app (full analysis + results)

```bash
streamlit run desktop.py
```

The app opens in your browser at `http://localhost:8501`.

**Run a new analysis:**
1. Load genes via the sidebar (upload, paste, or use sample data)
2. Set parameters (ID type, CPU cores, replication counts)
3. Click **Run GeneWalk**

**Explore existing results:**
1. Upload a `genewalk_results.csv` from a previous GeneWalk run
2. Results load instantly and all visualizations become available

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
| CPU cores | 4 | Parallel workers for DeepWalk. 2-4 for laptops, 8-16 for workstations. |
| Graph reps | 3 | DeepWalk repetitions on gene-GO network. 3 for testing, 10-15 for publication. |
| Null reps | 3 | DeepWalk repetitions on null networks. Should match graph reps. |
| FDR threshold | 1.0 | Set to 1.0 to keep all results, 0.05 for pre-filtering. |

## Sample Data

The `sample_data/` directory includes:

- `sample_genes.txt` -- 20 MAPK/ERK pathway genes
- `sample_genewalk_results.csv` -- Pre-computed GeneWalk results for instant demo

To regenerate sample data: `python generate_sample_data.py`

## Project Structure

```
GeneWalk-Analysis-Pipeline/
  app.py                      # Web app (Streamlit)
  desktop.py                  # Desktop app (Streamlit)
  desktop_launcher.py         # PyInstaller entry point
  build_desktop.py            # Build script for executable
  Dockerfile                  # Container deployment
  requirements.txt            # Web dependencies
  requirements-desktop.txt    # Desktop dependencies (includes genewalk)
  genewalk_app/
    __init__.py
    runner.py                 # GeneWalk execution & result parsing
    _gw_wrapper.py            # GeneWalk CLI wrapper (User-Agent + patches)
    dashboard.py              # Shared results dashboard
    visualizations.py         # Plotly charts (7 chart types)
    styles.py                 # Custom CSS theme
  sample_data/
    sample_genes.txt          # Demo gene list
    sample_genewalk_results.csv  # Demo results
  .streamlit/
    config.toml               # Streamlit theme & server config
```

## Citation

If you use this software, please cite:

```
Ietswaart R, Gyori BM, Bachman JA, Sorger PK, Bhatt DL, Churchman LS.
GeneWalk identifies relevant gene functions for a biological context using
network representation learning. Genome Biology 22, 55 (2021).
https://doi.org/10.1186/s13059-021-02264-8
```

See `CITATION.cff` for machine-readable citation metadata.

## License

MIT License. See [LICENSE](LICENSE) for details.
