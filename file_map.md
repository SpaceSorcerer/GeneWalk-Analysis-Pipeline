# File Map -- GeneWalk App

## Root Files

- **E:\Claude\GeneWalk-Analysis-Pipeline\app.py** -- Web app entry point (Streamlit). Provides sidebar toggle between "Run Analysis" and "View Results" modes, with GeneWalk execution and interactive dashboard.

- **E:\Claude\GeneWalk-Analysis-Pipeline\desktop.py** -- Desktop app entry point (Streamlit). Full GeneWalk runner with parameter configuration sidebar and results dashboard, bundled by PyInstaller.

- **E:\Claude\GeneWalk-Analysis-Pipeline\desktop_launcher.py** -- PyInstaller entry point that starts a local Streamlit server and opens the browser. Handles frozen bundle paths, port selection, and duplicate-tab prevention.

- **E:\Claude\GeneWalk-Analysis-Pipeline\build_desktop.py** -- Build script that invokes PyInstaller to create the desktop executable. Configures data files, collected packages, and hidden imports for the bundle.

- **E:\Claude\GeneWalk-Analysis-Pipeline\generate_sample_data.py** -- Script to regenerate the sample MAPK/ERK pathway demo data. Produces a realistic genewalk_results.csv with 20 genes and ~400 gene-GO associations.

- **E:\Claude\GeneWalk-Analysis-Pipeline\Dockerfile** -- Container deployment config for the web app. Installs web dependencies and runs Streamlit on port 8501.

- **E:\Claude\GeneWalk-Analysis-Pipeline\requirements.txt** -- Web app Python dependencies (Streamlit, Plotly, pandas, networkx, numpy). Does not include GeneWalk since the web app can visualize pre-computed results.

- **E:\Claude\GeneWalk-Analysis-Pipeline\requirements-desktop.txt** -- Desktop app Python dependencies. Includes everything in requirements.txt plus genewalk>=1.6 for running analyses locally.

- **E:\Claude\GeneWalk-Analysis-Pipeline\pyproject.toml** -- Python packaging metadata (PEP 621). Defines project name, version, dependencies, optional extras (desktop, dev), and ruff config.

- **E:\Claude\GeneWalk-Analysis-Pipeline\README.md** -- Project README with installation instructions (pip, executable, Docker), quick start guide, configuration reference, and project structure.

- **E:\Claude\GeneWalk-Analysis-Pipeline\CLAUDE.md** -- Project-level instructions for Claude Code. Describes architecture, key files, and build instructions.

- **E:\Claude\GeneWalk-Analysis-Pipeline\LICENSE** -- MIT License. Copyright 2024-2026 Brian Amburn / University of Texas Medical Branch.

- **E:\Claude\GeneWalk-Analysis-Pipeline\CITATION.cff** -- Machine-readable citation metadata (CFF format). Used by GitHub's "Cite this repository" feature.

- **E:\Claude\GeneWalk-Analysis-Pipeline\.gitignore** -- Git ignore rules for Python bytecode, virtual environments, IDE files, PyInstaller output, and GeneWalk run artifacts.

## genewalk_app/ Package

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\__init__.py** -- Package init file (empty). Marks genewalk_app as a Python package.

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\runner.py** -- Backend module for executing GeneWalk as a subprocess and parsing results. Provides gene list saving, progress streaming, result CSV loading, filtering, and per-gene summary functions.

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\_gw_wrapper.py** -- GeneWalk CLI wrapper that installs a proper HTTP User-Agent header and patches the GeneMapper class to handle MGI_EntrezGene.rpt format changes. Used as subprocess entry point.

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\dashboard.py** -- Shared results dashboard used by both app.py and desktop.py. Renders filters, summary metrics, interpretation guide, and five visualization tabs (Overview, Per-Gene Explorer, Network, Heatmap, Data Table).

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\visualizations.py** -- Seven Plotly visualization functions for GeneWalk results: volcano plot, gene bar chart, GO domain pie, heatmap, summary bar, p-value distribution, and gene-GO network graph.

- **E:\Claude\GeneWalk-Analysis-Pipeline\genewalk_app\styles.py** -- Custom CSS stylesheet for the Streamlit app. Provides theme-aware (light/dark) styling for hero banners, cards, tabs, metrics, filters, tables, and guides.

## Configuration

- **E:\Claude\GeneWalk-Analysis-Pipeline\.streamlit\config.toml** -- Streamlit theme and server configuration. Sets primary color, background colors, font, and max upload size.

## Sample Data

- **E:\Claude\GeneWalk-Analysis-Pipeline\sample_data\sample_genes.txt** -- Demo gene list with 20 MAPK/ERK pathway genes (BRAF, MAP2K1, MAPK1, etc.). Used for instant demo mode.

- **E:\Claude\GeneWalk-Analysis-Pipeline\sample_data\sample_genewalk_results.csv** -- Pre-computed GeneWalk results for the sample gene list. Contains ~400 gene-GO associations with similarity scores and p-values.

- **E:\Claude\GeneWalk-Analysis-Pipeline\sample_data\sample_deg_table.csv** -- Sample DEG table (legacy, not used by GeneWalk-only app). May be removed in a future cleanup.

## CI/CD

- **E:\Claude\GeneWalk-Analysis-Pipeline\.github\workflows\build-desktop.yml** -- GitHub Actions workflow that builds desktop executables for Windows, macOS, and Linux on version tags. Creates a GitHub Release with zipped artifacts.
