"""Splicing analysis parsers for vast-tools and rMATS outputs.

Converts outputs from both tools into a common internal DataFrame format
with standardized columns so that downstream visualizations and the
Gene Investigator work identically regardless of input source.

Common internal columns
-----------------------
gene        : str   – Gene symbol
event_id    : str   – Unique event identifier
event_type  : str   – SE / RI / A3SS / A5SS / MXE / MIC / other
coord       : str   – Genomic coordinates (chr:start-end)
psi_1       : float – Mean PSI in condition 1 (control / reference)
psi_2       : float – Mean PSI in condition 2 (treatment / test)
dpsi        : float – Delta PSI (psi_2 - psi_1)
pvalue      : float – Raw p-value
fdr         : float – Adjusted p-value / FDR
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


# ───────────────────────────────────────────────────────────────────────
# vast-tools
# ───────────────────────────────────────────────────────────────────────

# EVENT id prefix → event type mapping (species-agnostic)
_VT_EVENT_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"MIC", re.I), "MIC"),
    (re.compile(r"ALTA", re.I), "A3SS"),
    (re.compile(r"ALTD", re.I), "A5SS"),
    (re.compile(r"INT", re.I), "RI"),
    (re.compile(r"EX", re.I), "SE"),
]


def _classify_vt_event(event_id: str) -> str:
    """Map a vast-tools EVENT id to a canonical event type."""
    # Strip species prefix (e.g. Hsa, Mmu, Dre, ...)
    stripped = re.sub(r"^[A-Z][a-z]{2}", "", event_id)
    for pat, etype in _VT_EVENT_PATTERNS:
        if pat.search(stripped):
            return etype
    return "other"


def _detect_vasttools_columns(df: pd.DataFrame) -> dict[str, str | None]:
    """Auto-detect column names from a vast-tools table.

    Returns a dict mapping our internal names to the actual column names,
    or None if not found.
    """
    cols = {c.lower().strip(): c for c in df.columns}

    def _find(candidates: list[str]) -> str | None:
        for c in candidates:
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    return {
        "gene": _find(["gene", "genename", "gene_name", "symbol"]),
        "event_id": _find(["event", "event_id", "eventid"]),
        "coord": _find(["coord", "coordinates", "genomic_coord"]),
        "dpsi": _find([
            "dpsi", "delta_psi", "deltapsi", "dpsi_val",
            "e[dpsi]_per_rep",  # vast-tools compare output
        ]),
        "pvalue": _find(["p", "pvalue", "p-value", "p_value", "pval"]),
        "psi_1": _find([
            "meanpsi_grp1", "mean_psi_group1", "meanpsi_control",
            "psi_1", "psi1", "meanpsi_1",
        ]),
        "psi_2": _find([
            "meanpsi_grp2", "mean_psi_group2", "meanpsi_treatment",
            "psi_2", "psi2", "meanpsi_2",
        ]),
    }


def parse_vasttools(
    filepath_or_df: Path | str | pd.DataFrame,
    sep: str = "\t",
) -> pd.DataFrame:
    """Parse a vast-tools diff/compare output into the common format.

    Accepts either:
    * A file path to a tab-separated vast-tools table
    * A pre-loaded DataFrame
    """
    if isinstance(filepath_or_df, pd.DataFrame):
        raw = filepath_or_df.copy()
    else:
        raw = pd.read_csv(filepath_or_df, sep=sep)

    detected = _detect_vasttools_columns(raw)

    if not detected["gene"]:
        raise ValueError(
            "Cannot find a GENE column in this file. "
            f"Columns found: {', '.join(raw.columns.tolist())}"
        )

    out = pd.DataFrame()
    out["gene"] = raw[detected["gene"]].astype(str).str.strip()

    if detected["event_id"]:
        out["event_id"] = raw[detected["event_id"]].astype(str)
        out["event_type"] = out["event_id"].apply(_classify_vt_event)
    else:
        out["event_id"] = [f"vt_{i}" for i in range(len(raw))]
        out["event_type"] = "other"

    out["coord"] = (
        raw[detected["coord"]].astype(str) if detected["coord"]
        else ""
    )

    for internal, key in [("psi_1", "psi_1"), ("psi_2", "psi_2"),
                          ("dpsi", "dpsi"), ("pvalue", "pvalue")]:
        if detected[key]:
            out[internal] = pd.to_numeric(raw[detected[key]], errors="coerce")
        else:
            out[internal] = float("nan")

    # If dpsi missing but both PSIs present, compute it
    if out["dpsi"].isna().all() and not out["psi_1"].isna().all():
        out["dpsi"] = out["psi_2"] - out["psi_1"]

    # vast-tools doesn't always have FDR; fill with pvalue if missing
    out["fdr"] = out["pvalue"]
    out["source"] = "vast-tools"

    return out.dropna(subset=["gene"])


# ───────────────────────────────────────────────────────────────────────
# rMATS
# ───────────────────────────────────────────────────────────────────────

_RMATS_EVENT_TYPES = ["SE", "MXE", "A3SS", "A5SS", "RI"]


def _detect_rmats_columns(df: pd.DataFrame) -> dict[str, str | None]:
    """Auto-detect column names from an rMATS results file."""
    cols = {c.lower().strip(): c for c in df.columns}

    def _find(candidates: list[str]) -> str | None:
        for c in candidates:
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    return {
        "gene": _find([
            "genesymbol", "gene_symbol", "genename", "gene", "symbol",
        ]),
        "gene_id": _find(["geneid", "gene_id"]),
        "event_id": _find(["id"]),
        "chr": _find(["chr"]),
        "strand": _find(["strand"]),
        "pvalue": _find(["pvalue", "p-value"]),
        "fdr": _find(["fdr"]),
        "inc_level_1": _find(["inclevel1"]),
        "inc_level_2": _find(["inclevel2"]),
        "dpsi": _find(["incleveldifference"]),
    }


def _mean_inc_level(series: pd.Series) -> pd.Series:
    """Compute mean PSI from rMATS comma-separated IncLevel strings."""
    def _avg(val):
        if pd.isna(val):
            return float("nan")
        parts = str(val).split(",")
        nums = []
        for p in parts:
            p = p.strip()
            if p and p.lower() not in ("na", "nan", ""):
                try:
                    nums.append(float(p))
                except ValueError:
                    pass
        return sum(nums) / len(nums) if nums else float("nan")
    return series.apply(_avg)


def parse_rmats_file(
    filepath_or_df: Path | str | pd.DataFrame,
    event_type: str = "SE",
    sep: str = "\t",
) -> pd.DataFrame:
    """Parse a single rMATS event-type output file.

    Parameters
    ----------
    filepath_or_df : path or DataFrame
        e.g. SE.MATS.JC.txt
    event_type : str
        One of SE, MXE, A3SS, A5SS, RI.
    """
    if isinstance(filepath_or_df, pd.DataFrame):
        raw = filepath_or_df.copy()
    else:
        raw = pd.read_csv(filepath_or_df, sep=sep)

    detected = _detect_rmats_columns(raw)

    out = pd.DataFrame()

    # Gene name — prefer symbol, fall back to gene id
    if detected["gene"]:
        out["gene"] = raw[detected["gene"]].astype(str).str.strip()
    elif detected["gene_id"]:
        out["gene"] = raw[detected["gene_id"]].astype(str).str.strip()
    else:
        raise ValueError(
            "Cannot find gene column in rMATS file. "
            f"Columns: {', '.join(raw.columns.tolist())}"
        )

    # Event identity
    if detected["event_id"]:
        out["event_id"] = (
            event_type + "_" + raw[detected["event_id"]].astype(str)
        )
    else:
        out["event_id"] = [f"rmats_{event_type}_{i}" for i in range(len(raw))]

    out["event_type"] = event_type

    # Coordinates
    if detected["chr"]:
        chr_col = raw[detected["chr"]].astype(str)
        out["coord"] = chr_col
    else:
        out["coord"] = ""

    # PSI values
    if detected["inc_level_1"]:
        out["psi_1"] = _mean_inc_level(raw[detected["inc_level_1"]])
    else:
        out["psi_1"] = float("nan")

    if detected["inc_level_2"]:
        out["psi_2"] = _mean_inc_level(raw[detected["inc_level_2"]])
    else:
        out["psi_2"] = float("nan")

    # Delta PSI
    if detected["dpsi"]:
        out["dpsi"] = pd.to_numeric(raw[detected["dpsi"]], errors="coerce")
    elif not out["psi_1"].isna().all():
        out["dpsi"] = out["psi_2"] - out["psi_1"]
    else:
        out["dpsi"] = float("nan")

    # Statistics
    for internal, key in [("pvalue", "pvalue"), ("fdr", "fdr")]:
        if detected[key]:
            out[internal] = pd.to_numeric(raw[detected[key]], errors="coerce")
        else:
            out[internal] = float("nan")

    out["source"] = "rMATS"
    return out.dropna(subset=["gene"])


def parse_rmats_directory(
    directory: Path | str,
    use_jcec: bool = False,
) -> pd.DataFrame:
    """Parse all rMATS event-type files from an output directory.

    Parameters
    ----------
    directory : path
        rMATS output directory containing SE.MATS.JC.txt etc.
    use_jcec : bool
        If True, prefer JCEC (junction + exon body counts) over JC.
    """
    directory = Path(directory)
    suffix = "JCEC" if use_jcec else "JC"
    frames = []

    for event_type in _RMATS_EVENT_TYPES:
        filename = f"{event_type}.MATS.{suffix}.txt"
        filepath = directory / filename
        if filepath.exists():
            try:
                df = parse_rmats_file(filepath, event_type=event_type)
                frames.append(df)
            except (ValueError, pd.errors.EmptyDataError):
                pass

    if not frames:
        raise ValueError(
            f"No rMATS result files found in {directory}. "
            f"Expected files like SE.MATS.JC.txt"
        )

    return pd.concat(frames, ignore_index=True)


# ───────────────────────────────────────────────────────────────────────
# Unified helpers
# ───────────────────────────────────────────────────────────────────────

def filter_splicing_events(
    df: pd.DataFrame,
    dpsi_threshold: float = 0.1,
    fdr_threshold: float = 0.05,
    min_psi: float = 0.0,
    max_psi: float = 1.0,
) -> pd.DataFrame:
    """Filter splicing events by delta-PSI, FDR, and PSI range.

    If the FDR/p-value column is entirely NaN (common with vast-tools
    output that lacks statistical testing), the FDR filter is skipped
    and events are filtered by |delta-PSI| only.
    """
    out = df.copy()

    if "dpsi" in out.columns:
        out = out[out["dpsi"].abs() >= dpsi_threshold]

    # Only apply FDR filter when FDR values actually exist.
    # vast-tools diff output often lacks p-values entirely.
    if "fdr" in out.columns and out["fdr"].notna().any():
        # Keep events that either have FDR <= threshold OR have no FDR
        # (so that events without p-values aren't silently dropped).
        out = out[out["fdr"].isna() | (out["fdr"] <= fdr_threshold)]

    # PSI range filter (at least one condition must be in range)
    for col in ("psi_1", "psi_2"):
        if col in out.columns and min_psi > 0:
            mask = out[col].isna() | (out[col] >= min_psi)
            out = out[mask]
        if col in out.columns and max_psi < 1.0:
            mask = out[col].isna() | (out[col] <= max_psi)
            out = out[mask]

    return out


def splicing_summary(df: pd.DataFrame) -> dict:
    """Compute summary statistics for a splicing events DataFrame."""
    total = len(df)
    if total == 0:
        return {
            "total_events": 0,
            "unique_genes": 0,
            "event_types": {},
            "mean_abs_dpsi": 0,
        }

    event_counts = df["event_type"].value_counts().to_dict() if "event_type" in df.columns else {}

    return {
        "total_events": total,
        "unique_genes": df["gene"].nunique(),
        "event_types": event_counts,
        "mean_abs_dpsi": float(df["dpsi"].abs().mean()) if "dpsi" in df.columns and df["dpsi"].notna().any() else 0,
    }


def get_gene_splicing_events(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    """Get all splicing events for a specific gene."""
    if df.empty or "gene" not in df.columns:
        return pd.DataFrame()
    mask = df["gene"].str.upper() == gene.upper()
    return df[mask].copy()


def auto_detect_format(
    filepath: Path | str,
) -> str:
    """Detect whether a file is vast-tools or rMATS output.

    Returns "vasttools", "rmats", or "unknown".
    """
    filepath = Path(filepath)
    name = filepath.name.lower()

    # rMATS filenames are distinctive
    if any(name.startswith(et.lower()) for et in _RMATS_EVENT_TYPES):
        if ".mats." in name:
            return "rmats"

    # Try reading header
    try:
        header = pd.read_csv(filepath, sep="\t", nrows=0).columns.tolist()
        header_lower = [c.lower() for c in header]
    except Exception:
        return "unknown"

    # rMATS has IncLevelDifference; vast-tools has EVENT
    if "incleveldifference" in header_lower:
        return "rmats"
    if "event" in header_lower:
        return "vasttools"

    return "unknown"
