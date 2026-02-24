# GeneWalk Analysis Pipeline

A Streamlit web application for running [GeneWalk](https://github.com/churchmanlab/genewalk) and interactively exploring gene-function associations.

## Features

- **Upload or paste gene lists** -- supports text files, CSV, or direct pasting
- **Configurable parameters** -- gene ID type, CPU cores, replication counts, FDR threshold
- **Run GeneWalk** directly from the browser or **upload previous results**
- **Interactive visualizations** powered by Plotly:
  - Volcano plot (similarity vs. significance)
  - Per-gene GO term bar charts
  - Interactive gene-GO term network graph
  - GO domain distribution pie chart
  - Gene-GO term similarity heatmap
  - P-value distribution histogram
  - Gene summary ranking
- **Real-time filtering** by p-value threshold and GO domain
- **CSV export** of filtered results
- **Instant demo mode** -- sample MAPK/ERK pathway genes with pre-computed results load immediately

## Quick Start

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

> GeneWalk requires Python 3.8+ and at least 16 GB RAM for analysis runs.

### 2. Launch the app

```bash
streamlit run app.py
```

The app opens in your browser at `http://localhost:8501`.

### 3. Use the app

**Option A -- Run a new analysis:**
1. Load genes via the sidebar (upload, paste, or use sample data)
2. Set parameters (ID type, CPU cores, replication counts)
3. Click **Run GeneWalk**

**Option B -- Explore existing results:**
1. Upload a `genewalk_results.csv` file from a previous GeneWalk run
2. Results load instantly and all visualizations become available

## Project Structure

```
GeneWalk-Analysis-Pipeline/
  app.py                          # Main Streamlit application
  Dockerfile                      # Container deployment
  requirements.txt                # Python dependencies
  generate_sample_data.py         # Script to regenerate demo results
  genewalk_app/
    __init__.py
    runner.py                     # GeneWalk execution & result parsing
    visualizations.py             # Plotly visualization functions (7 chart types)
  sample_data/
    sample_genes.txt              # Example gene list (MAPK/ERK pathway)
    sample_genewalk_results.csv   # Pre-computed results for instant demo
  .streamlit/
    config.toml                   # Streamlit theme & server config
```

## Deployment

### Local (recommended for analysis runs)

```bash
streamlit run app.py
```

### Streamlit Community Cloud (for sharing results)

1. Push this repo to GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your repo and deploy

> Note: Running GeneWalk analysis requires significant compute. For cloud deployment, the "upload existing results" workflow is recommended.

### Docker

```bash
docker build -t genewalk-app .
docker run -p 8501:8501 genewalk-app
```

## Gene ID Types

| ID Type | Example | Source |
|---------|---------|--------|
| `hgnc_symbol` | BRAF | HUGO Gene Nomenclature |
| `hgnc_id` | HGNC:1097 | HUGO Gene Nomenclature |
| `ensembl_id` | ENSG00000157764 | Ensembl |
| `mgi_id` | MGI:88190 | Mouse Genome Informatics |
| `rgd_id` | RGD:2261 | Rat Genome Database |
| `entrez` | 673 | NCBI Entrez Gene |

## References

- Ietswaart R, Gyori BM, Bachman JA, Sorger PK, Bhatt DL, Churchman LS. GeneWalk identifies relevant gene functions for a biological context using network representation learning. *Genome Biology* 22, 55 (2021). https://doi.org/10.1186/s13059-021-02264-8
