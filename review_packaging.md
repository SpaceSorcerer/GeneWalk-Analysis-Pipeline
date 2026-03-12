# GeneWalk App -- Packaging & Release Readiness Review

Reviewed: 2026-03-11
Branch: `genewalk-solo`

---

## 1. README.md

**Status: PASS (clean)**

- Focused exclusively on GeneWalk. No mentions of GSEA, splicing, comparison, rMATS, DESeq, or any removed modules.
- Installation instructions are clear with three options (executable, pip, Docker).
- Quick-start guide covers both "Run Analysis" and "View Results" workflows.
- All referenced paths and commands are correct:
  - `sample_data/` exists with `sample_genes.txt` and `sample_genewalk_results.csv`
  - `requirements.txt` and `requirements-desktop.txt` exist
  - `generate_sample_data.py` exists
  - `CITATION.cff` and `LICENSE` exist
  - GitHub repo URL `SpaceSorcerer/GeneWalk-Analysis-Pipeline` is consistent throughout.
- Dashboard feature list (8 items in README) maps accurately to the 7 chart functions in `visualizations.py` plus the data table in `dashboard.py`.

**No issues found.**

---

## 2. pyproject.toml

### Issue 2.1: Invalid build-backend (SEVERITY: HIGH)

```toml
build-backend = "setuptools.backends._legacy:_Backend"
```

This import path does not exist in setuptools. It fails with `ModuleNotFoundError` even on setuptools 82.0.0. This will prevent `pip install .` and `pip install genewalk-app` from working.

**Fix:** Change to the standard setuptools build backend:
```toml
build-backend = "setuptools.build_meta"
```

### Issue 2.2: No entry points / console_scripts defined (SEVERITY: LOW)

There is no `[project.scripts]` or `[project.gui-scripts]` section. Users who `pip install genewalk-app` have no CLI command to run the app -- they must know to run `streamlit run app.py` manually. The README covers this, so it's not blocking, but a `genewalk-app` console script would improve usability.

**Suggestion:** Add:
```toml
[project.scripts]
genewalk-app = "genewalk_app._cli:main"
```
(with a small `_cli.py` that calls `streamlit run app.py`). Optional for v1.0.

### Other checks -- PASS:

- **Dependencies:** `streamlit>=1.30`, `plotly>=5.18`, `pandas>=2.0`, `networkx>=3.0`, `numpy` -- all exist on PyPI and are actively maintained. No removed-module dependencies (no gseapy, no scikit-learn, no scipy beyond what numpy covers).
- **Optional desktop dependency:** `genewalk>=1.6` -- exists on PyPI.
- **Dev dependencies:** `pytest>=7.0`, `ruff` -- both valid.
- **Python version:** `>=3.9` is reasonable. The code uses `list[str]` and `dict[str, str]` type hints (PEP 585), which work in 3.9+. Uses `str | None` union syntax (PEP 604), which requires 3.10+ at runtime.

### Issue 2.3: Python version mismatch with type syntax (SEVERITY: MEDIUM)

The codebase uses `X | Y` union type hints (e.g., `Path | None`, `str | None`) in function signatures in `runner.py`, `dashboard.py`, and `visualizations.py`. This syntax requires Python 3.10+ at runtime. The `pyproject.toml` claims `>=3.9`, and the CI matrix tests on 3.9. This will cause `TypeError` on Python 3.9.

**Fix:** Either:
- Change `requires-python = ">=3.10"` and drop 3.9 from the CI matrix, OR
- Add `from __future__ import annotations` to the affected files (runner.py, dashboard.py, visualizations.py) to make the annotations strings-only.

---

## 3. LICENSE

**Status: PASS**

- MIT License, correct full text.
- Attribution: "Copyright (c) 2024-2026 Brian Amburn / University of Texas Medical Branch"
- Properly formatted.

**No issues found.**

---

## 4. CITATION.cff

**Status: PASS**

- Valid CFF 1.2.0 format.
- Fields: `cff-version`, `message`, `type`, `title`, `authors`, `version`, `date-released`, `url`, `license`, `keywords` -- all present and correctly formatted.
- Version matches pyproject.toml (1.0.0).

**Minor note:** The `date-released` is `2026-03-11` (today). This is correct for a release being prepared now.

**No issues found.**

---

## 5. CI/CD Workflows

### `.github/workflows/test.yml` -- PASS

- Triggers on push and pull_request.
- Tests Python 3.9, 3.10, 3.11, 3.12.
- Installs via `pip install -e ".[dev]"` (uses pyproject.toml dev extras).
- Runs `pytest tests/ -v`.
- **Note:** Will fail on Python 3.9 due to Issue 2.3 above (union type syntax).
- **Note:** Will fail on all versions due to Issue 2.1 (invalid build-backend) since `pip install -e ".[dev]"` invokes the build system.

### `.github/workflows/build-desktop.yml` -- PASS (structure)

- Triggers on version tags (`v*`).
- Builds on Windows, macOS, Linux.
- Installs `requirements-desktop.txt` + `pyinstaller`.
- Runs `python build_desktop.py`.
- Uploads artifacts, creates GitHub release with zipped executables.
- Uses current action versions (checkout@v4, setup-python@v5, upload-artifact@v4, download-artifact@v4, softprops/action-gh-release@v2).

**No issues found in workflow logic.**

---

## 6. Build Script (`build_desktop.py`)

**Status: PASS (clean)**

- References only GeneWalk-related modules -- no removed modules.
- Data files: `desktop.py`, `genewalk_app/`, `sample_data/`, `.streamlit/` -- all exist.
- Hidden imports are all relevant:
  - `streamlit`, `plotly`, `pandas`, `networkx`, `numpy` -- core deps
  - `genewalk`, `genewalk.cli`, `genewalk.gene_lists`, `genewalk.resources` -- GeneWalk internals
  - `genewalk_app._gw_wrapper`, `genewalk_app.styles`, `genewalk_app.dashboard` -- app internals
  - `matplotlib` (used by GeneWalk internally for PDF plots)
  - `PIL`, `pyarrow`, `pkg_resources` -- common Streamlit hidden deps
- `--collect-all` includes: `streamlit`, `plotly`, `altair`, `pydeck`, `genewalk`, `matplotlib` -- all appropriate.
- `.streamlit/config.toml` exists and is bundled.

**No issues found.**

---

## 7. .gitignore

**Status: PASS**

Covers:
- Python bytecode (`__pycache__/`, `*.py[cod]`, `*.egg-info/`, etc.)
- Virtual environments (`venv/`, `.venv/`, `env/`)
- IDE files (`.idea/`, `.vscode/`, `*.swp`)
- OS files (`.DS_Store`, `Thumbs.db`)
- GeneWalk output (`genewalk_runs/`, `*.pkl`)
- PyInstaller (`dist/`, `*.spec`)
- Test cache (`.pytest_cache/`)
- Streamlit secrets (`.streamlit/secrets.toml`)

### Issue 7.1: Missing `build/` directory (SEVERITY: LOW)

The `build/` directory (created by `pip install -e .` or `python -m build`) is not in `.gitignore`. The `dist/` entry appears twice (lines 5 and 29) -- harmless but redundant.

**Fix:** Add `build/` to `.gitignore` (currently not listed). Remove duplicate `dist/` line.

---

## 8. requirements.txt and requirements-desktop.txt

**Status: PASS**

### requirements.txt (web viewer):
```
streamlit>=1.30
plotly>=5.18
pandas>=2.0
networkx>=3.0
numpy
```
Matches `pyproject.toml [project] dependencies` exactly.

### requirements-desktop.txt:
```
streamlit>=1.30
plotly>=5.18
pandas>=2.0
networkx>=3.0
numpy
genewalk>=1.6
```
Matches `pyproject.toml` base deps + `[desktop]` optional dependency.

- No removed-module dependencies (no gseapy, scipy, scikit-learn, etc.).
- All packages exist on PyPI.

**No issues found.**

---

## 9. Dockerfile

**Status: PASS**

- Based on `python:3.11-slim` -- appropriate.
- Installs `requirements.txt` only (web viewer, no GeneWalk) -- correct for a visualization-only container.
- Exposes port 8501 (Streamlit default).
- Healthcheck hits `/_stcore/health` -- correct Streamlit endpoint.
- Entrypoint: `streamlit run --server.port=8501 --server.address=0.0.0.0 app.py` -- correct.
- No references to removed modules.

**No issues found.**

---

## 10. Sample Data

**Status: MOSTLY PASS**

### `sample_data/sample_genes.txt` -- PASS
- 20 lines, one HGNC symbol per line (BRAF, MAP2K1, etc.).
- Used by both `app.py` and `desktop.py` for the "Use sample genes" input option.

### `sample_data/sample_genewalk_results.csv` -- PASS
- 459 data rows + 1 header = 460 lines.
- All expected columns present: `hgnc_symbol`, `hgnc_id`, `go_name`, `go_id`, `go_domain`, `sim`, `sem_sim`, `cilow`, `ciupp`, `global_padj`, `gene_padj`.
- Realistic MAPK/ERK pathway data with proper GO term IDs and domains.
- Can be regenerated with `python generate_sample_data.py`.

### Issue 10.1: Orphaned `sample_deg_table.csv` (SEVERITY: LOW)

`sample_data/sample_deg_table.csv` (51 lines) contains DEG data (baseMean, log2FoldChange, padj) that is not referenced anywhere in the GeneWalk app code. The `file_map.md` already flags it as "Legacy sample DEG table not used by the GeneWalk-only app. Candidate for removal."

**Fix:** Remove `sample_data/sample_deg_table.csv` before release.

---

## Summary of Issues

| # | Issue | Severity | File | Fix |
|---|-------|----------|------|-----|
| 2.1 | Invalid `build-backend` | **HIGH** | `pyproject.toml` | Change to `"setuptools.build_meta"` |
| 2.3 | `X \| Y` union syntax requires 3.10+, but claims 3.9 | **MEDIUM** | `pyproject.toml`, `runner.py`, `dashboard.py`, `visualizations.py` | Either bump to `>=3.10` or add `from __future__ import annotations` |
| 2.2 | No console_scripts entry point | LOW | `pyproject.toml` | Optional: add `[project.scripts]` |
| 7.1 | Missing `build/` in .gitignore; duplicate `dist/` | LOW | `.gitignore` | Add `build/`, remove duplicate `dist/` |
| 10.1 | Orphaned `sample_deg_table.csv` | LOW | `sample_data/` | Delete the file |

### Blocking for release:
- **Issue 2.1** (invalid build-backend) -- pip install will fail. Must fix.
- **Issue 2.3** (Python 3.9 incompatibility) -- CI will fail on 3.9, and any user on 3.9 will get import errors. Must fix one way or the other.

### Non-blocking:
- Issues 2.2, 7.1, 10.1 are nice-to-have cleanups.

### Clean areas (no issues):
- README.md -- focused, accurate, complete
- LICENSE -- valid MIT
- CITATION.cff -- valid format
- CI/CD workflows -- correct structure and references
- build_desktop.py -- no removed-module references, correct hidden imports
- requirements.txt / requirements-desktop.txt -- match pyproject.toml, no stale deps
- Dockerfile -- correct for web viewer
- All source code (app.py, desktop.py, desktop_launcher.py, genewalk_app/*) -- no references to GSEA, splicing, comparison, or any removed modules (one benign mention of "DESeq2" in desktop.py user-facing help text about where gene lists come from, which is appropriate context)
