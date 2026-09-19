"""Step 4 - flag the events worth reviewing and highlight them in Excel.

An event is *of interest* when:

* its length ``end - start + 1`` is at least ``min_length`` (3 bp), and
* its c. coordinate is in the coding sequence, i.e. contains none of ``*``
  (3' UTR), ``-`` (5' UTR or intron) and ``+`` (intron).

By default only the start coordinate is tested (``mode="start"``, the validated
rule); ``mode="any"`` accepts events whose start *or* end is coding, and
``mode="both"`` requires both. Without c. coordinates (step 3 skipped), use
``coding_only=False`` to flag on length alone.

Added columns: ``length`` and ``of_interest``; flagged rows are filled in yellow.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from eyeoncnv.common import (
    is_missing,
    log,
    pick_column,
    read_workbook,
    require_column,
    tagged_path,
    write_excel,
)

DEFAULT_MIN_LENGTH = 3
DEFAULT_MODE = "start"
DEFAULT_FILL = "FFFF00"
MODES = ("start", "any", "both")

_NON_CODING = re.compile(r"[*+\-]")
# c_hgvs_eq / c_coord are the column names written by the original notebooks.
START_C_CANDIDATES = ("c_start", "c_coord", "c_hgvs_eq")
END_C_CANDIDATES = ("c_end",)


def is_coding(c_coordinate) -> bool:
    """True for a c. coordinate inside the CDS (``467``), False for ``-45``, ``*23``,
    ``467+5``, ``100-12`` or a missing value. Full HGVS strings are accepted."""
    if is_missing(c_coordinate):
        return False
    text = str(c_coordinate).strip()
    if ":c." in text:
        text = text.split(":c.", 1)[1]
    return bool(text) and not _NON_CODING.search(text)


def event_length(df: pd.DataFrame) -> pd.Series:
    start = pd.to_numeric(df[require_column(df, ["start"], "start")], errors="coerce")
    end = pd.to_numeric(df[require_column(df, ["end"], "end")], errors="coerce")
    return end - start + 1


def interest_mask(
    df: pd.DataFrame,
    min_length: int = DEFAULT_MIN_LENGTH,
    mode: str = DEFAULT_MODE,
    coding_only: bool = True,
) -> pd.Series:
    """Boolean Series: which rows of ``df`` are events of interest."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    long_enough = (event_length(df) >= min_length).fillna(False).astype(bool)
    if not coding_only:
        return long_enough
    start_col = pick_column(df.columns, START_C_CANDIDATES)
    if start_col is None:
        raise KeyError(
            "No c. coordinate column (c_start): run step 3 (eyeoncnv hgvs) first, "
            "or disable the coding filter (--no-coding-filter)"
        )
    coding = df[start_col].map(is_coding).astype(bool)
    if mode != "start":
        end_col = pick_column(df.columns, END_C_CANDIDATES)
        end_coding = (
            df[end_col].map(is_coding).astype(bool)
            if end_col is not None
            else pd.Series(False, index=df.index)
        )
        coding = (coding | end_coding) if mode == "any" else (coding & end_coding)
    return (long_enough & coding).astype(bool)


def highlight_table(
    df: pd.DataFrame,
    min_length: int = DEFAULT_MIN_LENGTH,
    mode: str = DEFAULT_MODE,
    coding_only: bool = True,
) -> tuple[pd.DataFrame, pd.Series]:
    """Copy of ``df`` with ``length`` (after ``end``) and ``of_interest``, plus the mask."""
    out = df.copy()
    mask = interest_mask(out, min_length=min_length, mode=mode, coding_only=coding_only)
    length = event_length(out)
    if "length" in out.columns:
        out["length"] = length
    else:
        end_col = require_column(out, ["end"], "end")
        out.insert(out.columns.get_loc(end_col) + 1, "length", length)
    out["of_interest"] = mask.to_numpy()
    return out, mask


def run_highlight(
    input_path,
    output_path=None,
    min_length: int = DEFAULT_MIN_LENGTH,
    mode: str = DEFAULT_MODE,
    coding_only: bool = True,
    fill_color: str = DEFAULT_FILL,
) -> Path:
    """Step 4 on every sheet of ``input_path``; writes ``<input>_highlight.xlsx``."""
    input_path = Path(input_path)
    sheets = read_workbook(input_path)
    tables, masks = {}, {}
    for name, df in sheets.items():
        tables[name], mask = highlight_table(df, min_length=min_length, mode=mode, coding_only=coding_only)
        masks[name] = mask.tolist()
        log.info("Sheet %s: %d/%d event(s) of interest", name, int(mask.sum()), len(df))
    output = Path(output_path) if output_path else tagged_path(input_path, "highlight")
    write_excel(tables, output, highlights=masks, fill_color=fill_color)
    rule = f"length >= {min_length}"
    if coding_only:
        rule += f" and coding c. coordinate ({mode})"
    log.info("Written: %s (rule: %s)", output, rule)
    return output
