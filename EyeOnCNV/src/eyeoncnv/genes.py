"""Step 2 - annotate events with the regions of the capture BED.

Every event receives the names (4th BED column) of the regions it overlaps,
sorted and joined with ``;``. With a design BED whose names follow the
``GENE-NMxxxxxx`` convention, this gives the gene and transcript used by step 3.

Coordinates: events are 1-based inclusive ``[start, end]``; BED regions are
0-based half-open ``[start0, end0)``. The event is converted to
``[start - 1, end)`` before the overlap test.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from eyeoncnv.common import (
    chrom_key,
    is_missing,
    log,
    read_workbook,
    require_column,
    tagged_path,
    write_excel,
)

BED_SKIP_PREFIXES = ("#", "track", "browser")
CHR_CANDIDATES = ("chr", "chrom", "chromosome")
START_CANDIDATES = ("start",)
END_CANDIDATES = ("end",)
GENE_COLUMN = "gene"


def read_bed(path) -> pd.DataFrame:
    """Read a BED file (plain or gzipped) into ``chrom, start0, end0, name, chrom_key``.

    ``track``/``browser`` lines, comments and blank lines are skipped. Regions
    without a name get ``Unknown``; invalid regions (start < 0 or end <= start)
    are dropped.
    """
    path = Path(path)
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    records = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(BED_SKIP_PREFIXES):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                fields = stripped.split()
            if len(fields) < 3:
                raise ValueError(f"{path.name}, line {lineno}: fewer than 3 columns")
            try:
                start0, end0 = int(fields[1]), int(fields[2])
            except ValueError as exc:
                raise ValueError(f"{path.name}, line {lineno}: start/end are not integers") from exc
            name = fields[3].strip() if len(fields) > 3 else ""
            records.append((fields[0].strip(), start0, end0, name or "Unknown"))
    bed = pd.DataFrame(records, columns=["chrom", "start0", "end0", "name"])
    bed = bed[(bed["start0"] >= 0) & (bed["end0"] > bed["start0"])].reset_index(drop=True)
    bed["chrom_key"] = bed["chrom"].map(chrom_key)
    return bed


def _index_bed(bed: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    return {
        key: (
            sub["start0"].to_numpy(dtype="int64"),
            sub["end0"].to_numpy(dtype="int64"),
            sub["name"].to_numpy(dtype=object),
        )
        for key, sub in bed.groupby("chrom_key", sort=False)
    }


def annotate_genes(
    df: pd.DataFrame,
    bed: pd.DataFrame,
    out_column: str = GENE_COLUMN,
    separator: str = ";",
) -> pd.DataFrame:
    """Return a copy of ``df`` with the overlapping BED names in ``out_column``."""
    chr_col = require_column(df, CHR_CANDIDATES, "chr")
    start_col = require_column(df, START_CANDIDATES, "start")
    end_col = require_column(df, END_CANDIDATES, "end")
    index = _index_bed(bed)
    names = []
    for chrom, start, end in zip(df[chr_col], df[start_col], df[end_col]):
        key = chrom_key(chrom)
        start = pd.to_numeric(start, errors="coerce")
        end = pd.to_numeric(end, errors="coerce")
        if key not in index or is_missing(start) or is_missing(end):
            names.append("")
            continue
        start0, end0 = int(start) - 1, int(end)
        if end0 <= start0:
            names.append("")
            continue
        b_start, b_end, b_name = index[key]
        hits = (b_start < end0) & (b_end > start0)
        names.append(separator.join(sorted(set(b_name[hits]))))
    out = df.copy()
    out[out_column] = names
    return out


def run_genes(input_path, bed_path, output_path=None) -> Path:
    """Step 2 on every sheet of ``input_path``; writes ``<input>_genes.xlsx`` by default."""
    input_path = Path(input_path)
    bed = read_bed(bed_path)
    log.info("BED %s: %d region(s)", Path(bed_path).name, len(bed))
    sheets = read_workbook(input_path)
    annotated = {}
    for name, df in sheets.items():
        annotated[name] = annotate_genes(df, bed)
        n_hit = int((annotated[name][GENE_COLUMN] != "").sum())
        log.info("Sheet %s: %d/%d event(s) inside the design", name, n_hit, len(df))
    output = Path(output_path) if output_path else tagged_path(input_path, "genes")
    write_excel(annotated, output)
    log.info("Written: %s", output)
    return output
