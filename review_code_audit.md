# Code Audit Report -- GeneWalk App

**Date:** 2026-03-11
**Branch:** `genewalk-solo`
**Auditor:** Claude Opus 4.6
**Files audited:** 16 Python files (11 source, 5 test)

---

## Summary

The codebase is **clean and well-structured**. No critical or high-severity issues were found. There are no references to removed modules (comparison, gsea_runner, splicing, deg_visualizations, gene_investigator), no bare excepts, no deprecated Streamlit APIs, and no hardcoded secrets.

| Severity | Count |
|----------|-------|
| CRITICAL | 0 |
| HIGH | 0 |
| MEDIUM | 4 |
| LOW | 7 |

---

## Issues Found

### MEDIUM Severity

#### M1. Broad `except Exception` on CSV upload (app.py, desktop.py)

- **File:** `app.py` line 198; `desktop.py` line 198
- **Description:** When reading an uploaded CSV, `except Exception as exc` catches everything including `KeyboardInterrupt` (via BaseException subclasses that some libraries re-raise as Exception). While not a bare `except:`, catching all `Exception` on user file upload could silently swallow unexpected errors like `MemoryError`.
- **Suggested fix:** Narrow to `except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:`.

#### M2. Pandas `value_counts().reset_index()` column rename fragility (visualizations.py)

- **File:** `genewalk_app/visualizations.py` line 116-117
- **Description:**
  ```python
  counts = sig["go_domain"].value_counts().reset_index()
  counts.columns = ["go_domain", "count"]
  ```
  In pandas >= 2.0, `value_counts().reset_index()` already returns columns named `["go_domain", "count"]`, making the manual rename redundant. In pandas < 2.0, the columns were `["index", "go_domain"]`, so this rename was necessary. The code currently works on both versions but is fragile -- if pandas changes column naming again, or if the Series has no name, this could break.
- **Suggested fix:** Use the explicit pandas 2.0+ compatible pattern:
  ```python
  counts = sig["go_domain"].value_counts().reset_index(name="count")
  ```
  Since `pyproject.toml` already requires `pandas>=2.0`, this is safe.

#### M3. `run_genewalk` has no timeout for the subprocess (runner.py)

- **File:** `genewalk_app/runner.py` lines 143-198
- **Description:** The `Popen` call streams stderr line-by-line with no overall timeout. If GeneWalk hangs (e.g., network download stalls indefinitely), the subprocess runs forever. The `proc.wait(timeout=10)` on line 186 only applies after the stderr pipe is exhausted or an exception occurs in the for-loop -- it does not limit the total runtime.
- **Suggested fix:** Add an optional `timeout` parameter to `run_genewalk()`. Use a watchdog thread or `signal.alarm` to kill the process after a configurable maximum duration (e.g., 24 hours default).

#### M4. `desktop.py` writes a flag file at module import time (line 17)

- **File:** `desktop.py` line 17
- **Description:**
  ```python
  Path(tempfile.gettempdir(), ".genewalk_client_connected").write_text("")
  ```
  This runs unconditionally at import time (before any Streamlit setup), meaning importing `desktop.py` for testing or inspection always creates this file. This is a side effect on import that could interfere with the launcher's duplicate-tab detection logic in edge cases (e.g., running tests while the desktop app is launching).
- **Suggested fix:** Guard with `if __name__ == "__main__"` or move inside the Streamlit execution flow. Alternatively, check `st.runtime.exists()` before writing.

---

### LOW Severity

#### L1. Unused import: `sys` in generate_sample_data.py

- **File:** `generate_sample_data.py` line 5
- **Description:** `import sys` is present but `sys` is never used anywhere in the file.
- **Suggested fix:** Remove `import sys`.

#### L2. `except Exception` in runner.py stderr cleanup (lines 177, 183)

- **File:** `genewalk_app/runner.py` lines 177, 183
- **Description:** Two `except Exception: pass` blocks silently swallow errors during callback execution and stderr pipe cleanup. While intentional (to avoid killing the reader), they make debugging harder if something goes wrong.
- **Suggested fix:** Add `logging.debug()` calls inside these except blocks to at least log the swallowed exception.

#### L3. `except Exception` in desktop_launcher.py health check (line 80)

- **File:** `desktop_launcher.py` line 80
- **Description:** `except Exception:` catches all exceptions when polling the Streamlit health endpoint. This is acceptable for a retry loop but could mask unexpected errors.
- **Suggested fix:** Narrow to `except (urllib.error.URLError, OSError, TimeoutError):`.

#### L4. `_open_browser_once` is not thread-safe (desktop_launcher.py line 51-56)

- **File:** `desktop_launcher.py` lines 47-56
- **Description:** The `_browser_opened` flag is checked and set without a lock. In theory, two threads could both read `False` and both call `_original_wb_open`. In practice this is extremely unlikely since only one thread calls this, but it's technically a race condition.
- **Suggested fix:** Use `threading.Lock` or `threading.Event` to guard the flag.

#### L5. `st.plotly_chart` uses `width="stretch"` -- newer API (dashboard.py)

- **File:** `genewalk_app/dashboard.py` lines 242, 249, 260, 265, 285, 325, 358
- **Description:** The `width="stretch"` parameter for `st.plotly_chart` was introduced in Streamlit 1.33+. The `pyproject.toml` requires `streamlit>=1.30`, which means installations on Streamlit 1.30-1.32 would get a `TypeError` for the unexpected keyword argument.
- **Suggested fix:** Either bump the minimum Streamlit version in `pyproject.toml` to `>=1.33`, or use the older `use_container_width=True` parameter which works across all versions >=1.30.

#### L6. `st.dataframe` uses `width="stretch"` -- newer API (dashboard.py)

- **File:** `genewalk_app/dashboard.py` lines 293, 367
- **Description:** Same issue as L5 but for `st.dataframe`. The `width="stretch"` parameter was added in the same Streamlit version window.
- **Suggested fix:** Same as L5 -- bump minimum version or use `use_container_width=True`.

#### L7. `generate_sample_data.py` uses `open()` without explicit encoding (line 137)

- **File:** `generate_sample_data.py` line 137
- **Description:** `open(out, "w", newline="")` does not specify `encoding="utf-8"`. On Windows, this defaults to the system locale (cp1252). The sample data is ASCII-safe so this works in practice, but it's inconsistent with the rest of the codebase which explicitly uses UTF-8.
- **Suggested fix:** Add `encoding="utf-8"` to the `open()` call.

---

## Checks That Passed (No Issues Found)

### 1. Import Errors -- CLEAN
No file imports any removed module (`comparison`, `gsea_runner`, `splicing`, `deg_visualizations`, `gene_investigator`). All imports resolve to existing modules.

### 2. Dead Code -- CLEAN
No functions reference removed modules. No unused functions detected. All imports in source files are used (except `sys` in `generate_sample_data.py`, noted above).

### 3. Bare Excepts -- CLEAN
Zero instances of `except:` without exception type. All except clauses specify at least `Exception`.

### 4. Streamlit Deprecations -- CLEAN
No instances of `st.cache` (all would need `st.cache_data`), `st.experimental_*`, or `st.beta_*`.

### 5. Security Issues -- CLEAN
- No hardcoded secrets, tokens, API keys, or passwords.
- No `shell=True` in subprocess calls (command injection safe).
- Project names are sanitized via `_sanitize_project_name()` which strips path traversal characters.
- No user input is passed directly to shell commands.

### 6. Type Consistency -- CLEAN
- Function signatures use proper type hints (`str | None`, `list[str]`, `Path | None`).
- None checks are present where needed (e.g., `filter_results` checks if `padj_col in df.columns`).
- `render_dashboard` accepts `run_log: str | None` and checks `if run_log:` before use.

### 7. Missing Error Handling -- ACCEPTABLE
- All visualization functions handle empty DataFrames gracefully (return empty figures with descriptive titles).
- `load_results` validates required columns and raises `ValueError` with clear messages.
- CSV uploads are wrapped in try/except.
- Subprocess failures return structured error dicts rather than raising.

---

## Test Coverage Analysis

### What IS Tested (81 tests total)

| Module | Tests | Coverage |
|--------|-------|----------|
| `runner.py` | 33 | `save_gene_list`, `load_results`, `find_results_csv`, `filter_results`, `get_gene_summary`, `_sanitize_project_name`, `_parse_progress` |
| `visualizations.py` | 32 | All 7 chart functions with normal data, empty data, single gene, single GO term, missing columns, strict thresholds |
| `dashboard.py` (logic only) | 16 | Filter logic, metric computation, CSV export (via runner functions) |

### What is NOT Tested -- Coverage Gaps

#### Gap 1: `runner.py` -- `run_genewalk()` function
- **Not tested:** The `run_genewalk()` function that spawns the subprocess, streams stderr, and returns the result dict. This is the most complex function in the backend.
- **Reason:** Requires GeneWalk to be installed and takes minutes to run.
- **Recommendation:** Add a test that mocks `subprocess.Popen` to verify: (a) correct command construction, (b) progress callback invocation, (c) error handling when `FileNotFoundError` is raised, (d) cleanup when subprocess times out.

#### Gap 2: `runner.py` -- `is_genewalk_available()`
- **Not tested:** The function that checks if GeneWalk CLI is available.
- **Recommendation:** Add tests with mocked `subprocess.run` for: (a) successful case (returncode 0), (b) `FileNotFoundError`, (c) `TimeoutExpired`.

#### Gap 3: `runner.py` -- `_genewalk_base_cmd()`
- **Not tested:** Logic that switches between wrapper mode and frozen-bundle mode based on `sys.frozen`.
- **Recommendation:** Add tests mocking `sys.frozen` to verify correct command list in both modes.

#### Gap 4: `runner.py` -- `_utf8_env()`
- **Not tested:** Environment variable setup function.
- **Recommendation:** Simple function, low priority, but a one-line test would be easy to add.

#### Gap 5: `dashboard.py` -- `render_dashboard()`
- **Not tested:** The Streamlit rendering function itself.
- **Reason:** Requires a running Streamlit session. The underlying logic (filtering, metrics) is tested via runner imports.
- **Recommendation:** Consider using `streamlit.testing.v1.AppTest` (available since Streamlit 1.28) to test the full render path.

#### Gap 6: `dashboard.py` -- `_sorted_genes()`
- **Not tested:** The gene sorting function with its 5 different sort methods.
- **Recommendation:** This is a pure function that can be tested without Streamlit. Add tests for all 5 sort methods: Alphabetical, Input order, Best p-value, Most significant GO terms, Mean similarity.

#### Gap 7: `desktop_launcher.py` -- all functions
- **Not tested:** `_find_free_port()`, `_open_browser()`, `_resolve_app_path()`, `main()`, `_run_genewalk_subprocess()`.
- **Reason:** These are PyInstaller/launch infrastructure functions.
- **Recommendation:** `_find_free_port()` and `_resolve_app_path()` are testable pure functions. Add unit tests for those. The rest require integration testing.

#### Gap 8: `build_desktop.py` -- `find_package_dir()`, `build()`
- **Not tested:** Build infrastructure.
- **Recommendation:** Low priority. `find_package_dir()` could have a simple test.

#### Gap 9: `_gw_wrapper.py` -- `_install_opener()`, `_patch_gene_mapper()`
- **Not tested:** The GeneWalk CLI patches.
- **Reason:** Requires GeneWalk to be installed (imports `genewalk.gene_lists`).
- **Recommendation:** `_install_opener()` can be tested by verifying the User-Agent header is set. `_patch_gene_mapper()` could be tested with a mock if GeneWalk is available as a dev dependency.

#### Gap 10: `app.py` and `desktop.py` -- Streamlit UI logic
- **Not tested:** The main app entry points.
- **Reason:** Requires Streamlit session.
- **Recommendation:** Consider `AppTest` for smoke testing the UI flow.

---

## Overall Assessment

The GeneWalk App codebase is **production-ready** for its current scope. The code is clean, well-documented, and follows good practices:

- Clear separation of concerns (runner, dashboard, visualizations, styles)
- Defensive programming in visualization functions (empty data handling)
- Input sanitization for project names
- No security vulnerabilities
- No deprecated API usage
- Comprehensive test coverage for core logic (81 tests)

The main areas for improvement are:
1. **Bump Streamlit minimum version** to `>=1.33` in `pyproject.toml` (to match actual API usage)
2. **Add mock-based tests** for `run_genewalk()`, `is_genewalk_available()`, and `_sorted_genes()`
3. **Minor cleanup**: remove unused `sys` import in `generate_sample_data.py`, add UTF-8 encoding to its `open()` call
