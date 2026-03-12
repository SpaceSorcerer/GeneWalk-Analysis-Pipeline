from __future__ import annotations

"""Backend module for running GeneWalk and parsing results."""

import logging
import os
import re
import subprocess
import sys
import tempfile
import threading
from collections.abc import Callable
from pathlib import Path

logger = logging.getLogger(__name__)

import pandas as pd


def _utf8_env() -> dict[str, str]:
    """Return a copy of the current environment with PYTHONUTF8=1.

    On Windows, Python defaults to cp1252 encoding for file I/O.  GeneWalk's
    HTML report generation uses Unicode characters that cp1252 cannot encode,
    causing a crash.  Setting PYTHONUTF8=1 forces Python to use UTF-8 for all
    file I/O in the subprocess.
    """
    env = os.environ.copy()
    env["PYTHONUTF8"] = "1"
    return env


DEFAULT_BASE_DIR = Path(tempfile.gettempdir()) / "genewalk_runs"


_WRAPPER = str(Path(__file__).with_name("_gw_wrapper.py"))


# Mapping of GeneWalk log messages to user-friendly progress descriptions.
_PROGRESS_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"Downloading", re.I), "Downloading resources..."),
    (re.compile(r"Reading genes", re.I), "Reading gene list..."),
    (re.compile(r"Loading|Assembling.*network", re.I), "Building gene network..."),
    (re.compile(r"Get GO annotations", re.I), "Fetching GO annotations..."),
    (re.compile(r"DeepWalk.*graph", re.I), "Running DeepWalk on gene network..."),
    (re.compile(r"DeepWalk.*null", re.I), "Running DeepWalk on null model..."),
    (re.compile(r"Similarities", re.I), "Computing similarity scores..."),
    (re.compile(r"Perform.*statistic|p-value|Significance", re.I), "Computing statistics..."),
    (re.compile(r"generate_plots|Barplot generated", re.I), "Generating plots..."),
    (re.compile(r"make_html|HTML", re.I), "Writing HTML report..."),
    (re.compile(r"genewalk_results", re.I), "Writing results CSV..."),
]


def _parse_progress(line: str) -> str | None:
    """Extract a user-friendly progress message from a GeneWalk log line."""
    for pattern, message in _PROGRESS_PATTERNS:
        if pattern.search(line):
            return message
    return None


def _genewalk_base_cmd() -> list[str]:
    """Return the base command list for invoking GeneWalk.

    Always runs through our wrapper script so that a proper HTTP User-Agent
    header is set before GeneWalk tries to download any resources (the default
    ``urllib`` User-Agent is blocked by several servers with HTTP 403).

    Using ``sys.executable`` also avoids the Windows PATH issue where pip's
    ``Scripts/`` directory isn't on PATH (common with Microsoft Store Python).

    In a PyInstaller frozen bundle, ``sys.executable`` points to the ``.exe``
    rather than a Python interpreter, so we pass a ``--run-genewalk`` flag that
    tells the launcher to dispatch directly to the GeneWalk CLI instead of
    starting another Streamlit server.
    """
    if getattr(sys, "frozen", False):
        return [sys.executable, "--run-genewalk"]
    return [sys.executable, _WRAPPER]


def is_genewalk_available() -> bool:
    """Return True if the GeneWalk CLI can be invoked."""
    try:
        cmd = _genewalk_base_cmd() + ["--help"]
        result = subprocess.run(cmd, capture_output=True, timeout=15, env=_utf8_env())
        return result.returncode == 0
    except (FileNotFoundError, PermissionError, subprocess.TimeoutExpired):
        return False


def _sanitize_project_name(name: str) -> str:
    """Sanitize a project name to prevent path traversal and bad chars."""
    # Strip leading/trailing whitespace and path separators
    name = name.strip().strip("/").strip("\\")
    # Replace path separators and other problematic characters
    name = re.sub(r'[/\\:*?"<>|.\x00]', "_", name)
    # Collapse repeated underscores
    name = re.sub(r"_+", "_", name).strip("_")
    return name or "genewalk_analysis"


def save_gene_list(genes: list[str], dest: Path) -> Path:
    """Write a list of gene identifiers to a text file."""
    cleaned = [g.strip() for g in genes if g.strip()]
    if not cleaned:
        raise ValueError("Gene list is empty after removing blank entries.")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(cleaned) + "\n", encoding="utf-8")
    return dest


def run_genewalk(
    gene_file: Path,
    project: str,
    id_type: str = "hgnc_symbol",
    nproc: int = 4,
    nreps_graph: int = 3,
    nreps_null: int = 3,
    alpha_fdr: float = 1.0,
    base_folder: Path | None = None,
    on_progress: Callable[[str], None] | None = None,
    timeout: float = 7200,
) -> dict:
    """Run GeneWalk CLI and return paths to output files.

    Parameters
    ----------
    on_progress : callable, optional
        Called with a short status string whenever GeneWalk enters a new
        major step (e.g. "Building gene network...", "Computing statistics...").
    timeout : float
        Maximum wall-clock seconds before the subprocess is killed.
        Defaults to 7200 (2 hours).  Set to ``0`` or ``None`` to disable.

    Returns a dict with keys: 'return_code', 'stdout', 'stderr', 'output_dir'.
    """
    base = base_folder or DEFAULT_BASE_DIR
    base.mkdir(parents=True, exist_ok=True)

    cmd = _genewalk_base_cmd() + [
        "--project", project,
        "--genes", str(gene_file),
        "--id_type", id_type,
        "--nproc", str(nproc),
        "--nreps_graph", str(nreps_graph),
        "--nreps_null", str(nreps_null),
        "--alpha_fdr", str(alpha_fdr),
        "--base_folder", str(base),
    ]

    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            env=_utf8_env(),
        )
    except FileNotFoundError:
        return {
            "return_code": 1,
            "stdout": "",
            "stderr": (
                "ERROR: Could not find the 'genewalk' command.\n\n"
                "Install it with:  pip install genewalk\n\n"
                "If you already installed it, make sure it's in the same "
                "Python environment running this app.\n"
                f"Current Python: {sys.executable}"
            ),
            "output_dir": base / project,
        }

    # ---- Watchdog timer: kill the subprocess if it exceeds the timeout ----
    timed_out = False
    watchdog: threading.Timer | None = None
    if timeout:
        def _kill_on_timeout() -> None:
            nonlocal timed_out
            timed_out = True
            try:
                proc.kill()
            except OSError:
                pass

        watchdog = threading.Timer(timeout, _kill_on_timeout)
        watchdog.daemon = True
        watchdog.start()

    # Stream stderr line-by-line and extract progress updates.
    stderr_lines: list[str] = []
    last_status: str | None = None
    try:
        for line in proc.stderr:
            stderr_lines.append(line)
            if on_progress:
                try:
                    status = _parse_progress(line)
                    if status and status != last_status:
                        last_status = status
                        on_progress(status)
                except Exception:
                    logger.debug("Progress callback error", exc_info=True)
    finally:
        if watchdog is not None:
            watchdog.cancel()
        # Ensure the subprocess is cleaned up even on exceptions
        try:
            proc.stderr.close()
        except Exception:
            logger.debug("Error closing stderr pipe", exc_info=True)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()

    output_dir = base / project

    if timed_out:
        stderr_lines.append(
            f"\nERROR: GeneWalk subprocess killed after {timeout}s timeout.\n"
        )

    return {
        "return_code": proc.returncode,
        "stdout": "",
        "stderr": "".join(stderr_lines),
        "output_dir": output_dir,
    }


def find_results_csv(output_dir: Path) -> Path | None:
    """Locate the genewalk_results.csv file in the output directory."""
    candidate = output_dir / "genewalk_results.csv"
    if candidate.exists():
        return candidate
    # Search recursively as a fallback
    for p in output_dir.rglob("genewalk_results.csv"):
        return p
    return None


_REQUIRED_COLUMNS = {"hgnc_symbol", "go_name", "sim"}


def load_results(csv_path: Path) -> pd.DataFrame:
    """Load and clean GeneWalk results CSV.

    Raises
    ------
    ValueError
        If the CSV is empty or missing required columns.
    """
    try:
        df = pd.read_csv(csv_path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError) as exc:
        raise ValueError(f"Could not parse results CSV: {exc}") from exc
    if df.empty:
        raise ValueError("Results CSV is empty (0 rows).")
    missing = _REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Results CSV is missing required columns: {', '.join(sorted(missing))}. "
            f"Found columns: {', '.join(df.columns.tolist())}"
        )
    return df


def filter_results(
    df: pd.DataFrame,
    padj_col: str = "gene_padj",
    padj_threshold: float = 0.05,
    go_domain: str | None = None,
) -> pd.DataFrame:
    """Filter results by adjusted p-value and optional GO domain."""
    filtered = df.copy()
    if padj_col in filtered.columns:
        filtered = filtered[filtered[padj_col] <= padj_threshold]
    if go_domain and "go_domain" in filtered.columns:
        filtered = filtered[filtered["go_domain"] == go_domain]
    return filtered.sort_values(padj_col, ascending=True) if padj_col in filtered.columns else filtered


def get_gene_summary(df: pd.DataFrame, padj_col: str = "gene_padj", padj_threshold: float = 0.05) -> pd.DataFrame:
    """Summarize results per gene: count of significant GO terms."""
    if padj_col not in df.columns:
        return pd.DataFrame()
    required = {"hgnc_symbol", "go_name", "sim"}
    if not required.issubset(df.columns):
        return pd.DataFrame()
    sig = df[df[padj_col] <= padj_threshold].sort_values(padj_col, ascending=True)
    if sig.empty:
        return pd.DataFrame()
    summary = (
        sig.groupby("hgnc_symbol")
        .agg(
            significant_go_terms=("go_name", "count"),
            top_go_term=("go_name", "first"),
            best_padj=(padj_col, "min"),
            mean_similarity=("sim", "mean"),
        )
        .sort_values("significant_go_terms", ascending=False)
        .reset_index()
    )
    return summary
