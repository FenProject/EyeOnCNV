"""Step 1 - from depth-jump positions to events, with recurrence across samples.

Input: the TSV files of step 0 (``chromosome, position, depth``; no header).

Positions of a sample are grouped into intervals (a new interval starts when two
consecutive positions of a chromosome are ``gap`` bp apart or more). Step 0 prints
positions in pairs around each breakpoint - the last base before the depth change
and the first base after it - so for an interval of 3 positions or more the event
spans from the 2nd to the second-to-last position. The depth of the 2nd position
compared with the 1st tells whether depth goes up (duplication) or down (deletion).

For every event, ``occurrences`` is the number of *other* samples having at least
one event on the same chromosome that covers ``min_overlap`` (5 %) of its length.
Samples from previous runs (per-sample Excel files already present in the history
directories) are included in that count, so recurrent artefacts stand out.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from eyeoncnv.common import chrom_key, collect_files, log, write_excel

DEFAULT_GAP = 150
DEFAULT_MIN_OVERLAP = 0.05
DEFAULT_MAX_OCCURRENCES = 1

DUPLICATION = "duplication"
DELETION = "deletion"

PER_SAMPLE_SUFFIX = "_intervals_with_variation.xlsx"
INTERVAL_COLUMNS = ["chr", "start", "end", "variation"]
SUMMARY_COLUMNS = ["sample", "chr", "start", "end", "variation", "occurrences"]

# Upper bound of the (query x other) overlap matrix computed at once.
_MAX_MATRIX_CELLS = 4_000_000


def read_depth_tsv(path) -> pd.DataFrame:
    """Read a step-0 TSV into columns ``chr``, ``position`` (int) and ``depth`` (float)."""
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            dtype={0: "string", 1: "float64", 2: "float64"},
        )
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(
            {
                0: pd.Series(dtype="string"),
                1: pd.Series(dtype="float64"),
                2: pd.Series(dtype="float64"),
            }
        )
    df.columns = ["chr", "position", "depth"]
    df = df.dropna(subset=["chr", "position"]).copy()
    df["position"] = df["position"].round().astype("int64")
    df["depth"] = df["depth"].astype("float64")
    return df


def _interval_record(chrom, positions: np.ndarray, depths: np.ndarray) -> dict:
    n = len(positions)
    if n >= 3:
        start, end = positions[1], positions[-2]
    elif n == 2:
        start, end = positions[0], positions[1]
    else:
        start = end = positions[0]
    if n >= 2:
        variation = DUPLICATION if depths[1] > depths[0] else DELETION
    else:
        variation = ""
    return {"chr": chrom, "start": int(start), "end": int(end), "variation": variation}


def build_intervals(df: pd.DataFrame, gap: int = DEFAULT_GAP) -> pd.DataFrame:
    """Group depth-jump positions into events (see module docstring for the rules).

    ``df`` needs the columns ``chr``, ``position`` and ``depth``. Two consecutive
    positions belong to the same event when their distance is strictly below ``gap``.
    """
    records = []
    for chrom, sub in df.groupby("chr", sort=False):
        sub = sub.sort_values("position", kind="mergesort")
        positions = sub["position"].to_numpy(dtype="int64")
        depths = sub["depth"].to_numpy(dtype="float64")
        if positions.size == 0:
            continue
        breaks = np.flatnonzero(np.diff(positions) >= gap) + 1
        for block in np.split(np.arange(positions.size), breaks):
            records.append(_interval_record(chrom, positions[block], depths[block]))
    return pd.DataFrame(records, columns=INTERVAL_COLUMNS)


def _by_chromosome(df: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    keys = df["chr"].map(chrom_key)
    groups: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for key, sub in df.assign(_key=keys).dropna(subset=["_key"]).groupby("_key", sort=False):
        groups[key] = (
            sub["start"].to_numpy(dtype="int64"),
            sub["end"].to_numpy(dtype="int64"),
        )
    return groups


def count_recurrence(
    query: pd.DataFrame,
    others: Mapping[str, pd.DataFrame],
    min_overlap: float = DEFAULT_MIN_OVERLAP,
) -> np.ndarray:
    """Number of samples in ``others`` with an event covering each ``query`` event.

    Intervals are 1-based inclusive. A sample counts once per query event when at
    least one of its events on the same chromosome overlaps the query event by
    ``min_overlap`` x (query length) bases or more (and by at least 1 bp).
    """
    grouped = [_by_chromosome(df) for df in others.values()]
    return _count_recurrence_grouped(query, grouped, min_overlap)


def _count_recurrence_grouped(
    query: pd.DataFrame,
    others: Iterable[dict[str, tuple[np.ndarray, np.ndarray]]],
    min_overlap: float,
) -> np.ndarray:
    n = len(query)
    counts = np.zeros(n, dtype="int64")
    if n == 0:
        return counts
    q_start = query["start"].to_numpy(dtype="int64")
    q_end = query["end"].to_numpy(dtype="int64")
    q_len = np.maximum(q_end - q_start + 1, 1)
    q_rows: dict[str, list[int]] = {}
    for row, key in enumerate(query["chr"].map(chrom_key)):
        if key is not None:
            q_rows.setdefault(key, []).append(row)

    for other_groups in others:
        hit = np.zeros(n, dtype=bool)
        for key, rows in q_rows.items():
            if key not in other_groups:
                continue
            o_start, o_end = other_groups[key]
            rows_arr = np.asarray(rows)
            chunk = max(1, _MAX_MATRIX_CELLS // max(o_start.size, 1))
            for lo in range(0, rows_arr.size, chunk):
                sel = rows_arr[lo : lo + chunk]
                left = np.maximum(q_start[sel, None], o_start[None, :])
                right = np.minimum(q_end[sel, None], o_end[None, :])
                overlap = np.maximum(right - left + 1, 0)
                covered = (overlap > 0) & (overlap / q_len[sel, None] >= min_overlap)
                hit[sel] = covered.any(axis=1)
        counts += hit
    return counts


def load_history(directories: Iterable, exclude_samples: Iterable[str] = ()) -> dict[str, pd.DataFrame]:
    """Events of previous runs, read from ``*_intervals_with_variation.xlsx`` files.

    Files of samples listed in ``exclude_samples`` are skipped: they are about to be
    rewritten by the current run and must not be counted against themselves.
    """
    excluded = set(exclude_samples)
    pool: dict[str, pd.DataFrame] = {}
    seen_dirs: set[Path] = set()
    for directory in directories:
        directory = Path(directory)
        if not directory.is_dir() or directory.resolve() in seen_dirs:
            continue
        seen_dirs.add(directory.resolve())
        for xlsx in sorted(directory.glob(f"*{PER_SAMPLE_SUFFIX}")):
            if xlsx.name.startswith("~$"):
                continue
            sample = xlsx.name[: -len(PER_SAMPLE_SUFFIX)]
            if sample in excluded:
                log.debug("History: %s skipped (sample processed in this run)", xlsx.name)
                continue
            try:
                df = pd.read_excel(xlsx, dtype={"chr": "string"})
            except Exception as exc:  # unreadable file: warn and go on
                log.warning("History: cannot read %s (%s)", xlsx, exc)
                continue
            if not {"chr", "start", "end"}.issubset(df.columns):
                log.warning("History: %s has no chr/start/end columns, ignored", xlsx.name)
                continue
            events = pd.DataFrame(
                {
                    "chr": df["chr"].astype("string"),
                    "start": pd.to_numeric(df["start"], errors="coerce"),
                    "end": pd.to_numeric(df["end"], errors="coerce"),
                }
            ).dropna()
            pool[f"history:{xlsx.resolve()}"] = events.astype({"start": "int64", "end": "int64"})
    return pool


@dataclass
class IntervalsResult:
    """Paths written by :func:`run_intervals`."""

    summary: Path
    per_sample: dict[str, Path] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)
    history_samples: int = 0


def summary_filename(max_occurrences: int) -> str:
    return f"summary_occurrences_0-{max_occurrences}.xlsx"


def run_intervals(
    inputs: Sequence,
    out_dir,
    gap: int = DEFAULT_GAP,
    min_overlap: float = DEFAULT_MIN_OVERLAP,
    max_occurrences: int = DEFAULT_MAX_OCCURRENCES,
    history_dirs: Sequence | None = None,
    use_history: bool = True,
) -> IntervalsResult:
    """Run step 1 on TSV files (or directories of TSV files).

    Writes ``<sample>_intervals_with_variation.xlsx`` for every sample and
    ``summary_occurrences_0-<max_occurrences>.xlsx`` with one sheet per occurrence
    value (``occ_0``, ``occ_1``...). ``history_dirs`` defaults to ``out_dir``.
    """
    if max_occurrences < 0:
        raise ValueError("max_occurrences must be >= 0")
    out_dir = Path(out_dir)
    tsv_files = collect_files(inputs, (".tsv",))
    if not tsv_files:
        raise FileNotFoundError("No .tsv file found in: " + ", ".join(map(str, inputs)))

    intervals: dict[str, pd.DataFrame] = {}
    errors: dict[str, str] = {}
    for tsv in tsv_files:
        sample = tsv.stem
        if sample in intervals:
            raise ValueError(f"Two input files share the sample name '{sample}'")
        try:
            log.info("Reading %s", tsv)
            intervals[sample] = build_intervals(read_depth_tsv(tsv), gap=gap)
            if intervals[sample].empty:
                log.info("%s: no event", sample)
        except Exception as exc:  # keep going with the other samples
            errors[sample] = str(exc)
            log.error("%s: %s", tsv, exc)
    if not intervals:
        raise RuntimeError("No TSV file could be processed")

    history: dict[str, pd.DataFrame] = {}
    if use_history:
        directories = list(history_dirs) if history_dirs else [out_dir]
        history = load_history(directories, exclude_samples=intervals.keys())
        if history:
            log.info("%d sample(s) from previous runs added to the recurrence count", len(history))

    out_dir.mkdir(parents=True, exist_ok=True)
    result = IntervalsResult(
        summary=out_dir / summary_filename(max_occurrences),
        errors=errors,
        history_samples=len(history),
    )
    grouped = {name: _by_chromosome(df) for name, df in intervals.items()}
    grouped_history = [_by_chromosome(df) for df in history.values()]
    kept = []
    for sample, events in intervals.items():
        others = [g for name, g in grouped.items() if name != sample] + grouped_history
        table = events.copy()
        table["occurrences"] = _count_recurrence_grouped(events, others, min_overlap)
        path = write_excel({"intervals": table}, out_dir / f"{sample}{PER_SAMPLE_SUFFIX}")
        result.per_sample[sample] = path
        log.info("%s: %d event(s) -> %s", sample, len(table), path.name)
        selected = table[table["occurrences"] <= max_occurrences]
        if not selected.empty:
            kept.append(selected.assign(sample=sample)[SUMMARY_COLUMNS])

    summary = pd.concat(kept, ignore_index=True) if kept else pd.DataFrame(columns=SUMMARY_COLUMNS)
    sheets = {
        f"occ_{k}": summary[summary["occurrences"] == k].reset_index(drop=True)
        for k in range(max_occurrences + 1)
    }
    write_excel(sheets, result.summary)
    log.info(
        "Summary: %d event(s) with occurrences <= %d -> %s",
        len(summary),
        max_occurrences,
        result.summary,
    )
    if errors:
        log.warning("%d file(s) could not be processed: %s", len(errors), ", ".join(errors))
    return result
