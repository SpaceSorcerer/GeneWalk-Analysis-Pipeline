# GeneWalk App

## Description
Focused desktop/web app for running GeneWalk and exploring results interactively. Provides a Streamlit-based UI for configuring and executing GeneWalk analyses, uploading previous results, and visualizing gene-GO term associations with interactive Plotly charts.

## Architecture
- **UI Framework:** Streamlit (web and desktop modes)
- **Visualization:** Plotly (volcano plots, networks, heatmaps, bar charts)
- **Analysis Engine:** GeneWalk subprocess runner with progress streaming
- **Styling:** Custom CSS theme with light/dark mode support

### Key Components
- `app.py` -- Web app entry point (run analysis + view results with sidebar toggle)
- `desktop.py` -- Desktop app entry point (bundled by PyInstaller)
- `desktop_launcher.py` -- PyInstaller launcher (starts Streamlit server + opens browser)
- `genewalk_app/runner.py` -- GeneWalk subprocess execution and result parsing
- `genewalk_app/_gw_wrapper.py` -- GeneWalk CLI wrapper (User-Agent fix + GeneMapper patch)
- `genewalk_app/dashboard.py` -- Shared results dashboard (used by both app.py and desktop.py)
- `genewalk_app/visualizations.py` -- 7 Plotly chart types for GeneWalk results
- `genewalk_app/styles.py` -- Custom CSS theme with dark mode support

## Key Files
- `requirements.txt` -- Web app dependencies (no GeneWalk needed)
- `requirements-desktop.txt` -- Desktop dependencies (includes GeneWalk)
- `build_desktop.py` -- PyInstaller build script
- `.github/workflows/build-desktop.yml` -- CI/CD for building executables
- `sample_data/` -- Demo gene list and pre-computed results

## Build Instructions

### Run locally (web)
```bash
pip install -r requirements.txt
streamlit run app.py
```

### Run locally (desktop with analysis)
```bash
pip install -r requirements-desktop.txt
streamlit run desktop.py
```

### Build executable (PyInstaller)
```bash
pip install -r requirements-desktop.txt
pip install pyinstaller
python build_desktop.py
```

### Docker
```bash
docker build -t genewalk-app .
docker run -p 8501:8501 genewalk-app
```
