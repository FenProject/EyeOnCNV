"""Steps 1 to 4 chained, as run by ``eyeoncnv run``."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from eyeoncnv import genes, hgvs_coords, highlight, intervals
from eyeoncnv.common import log


@dataclass
class PipelineResult:
    summary: Path
    genes: Path
    hgvs: Path | None
    final: Path


def run_pipeline(
    inputs: Sequence,
    out_dir,
    bed,
    transcript_map=None,
    *,
    gap: int = intervals.DEFAULT_GAP,
    min_overlap: float = intervals.DEFAULT_MIN_OVERLAP,
    max_occurrences: int = intervals.DEFAULT_MAX_OCCURRENCES,
    history_dirs: Sequence | None = None,
    use_history: bool = True,
    assembly: str = hgvs_coords.DEFAULT_ASSEMBLY,
    db_url: str | None = None,
    min_length: int = highlight.DEFAULT_MIN_LENGTH,
    mode: str = highlight.DEFAULT_MODE,
    fill_color: str = highlight.DEFAULT_FILL,
    projector: hgvs_coords.Projector | None = None,
) -> PipelineResult:
    """Run steps 1 -> 4. Without ``transcript_map``, step 3 is skipped and step 4
    flags events on length only."""
    log.info("Step 1/4 - events and recurrence")
    step1 = intervals.run_intervals(
        inputs,
        out_dir,
        gap=gap,
        min_overlap=min_overlap,
        max_occurrences=max_occurrences,
        history_dirs=history_dirs,
        use_history=use_history,
    )
    log.info("Step 2/4 - genes from the capture BED")
    step2 = genes.run_genes(step1.summary, bed)
    step3 = None
    if transcript_map is not None:
        log.info("Step 3/4 - c. coordinates (hgvs)")
        step3 = hgvs_coords.run_hgvs(
            step2, transcript_map, assembly=assembly, db_url=db_url, projector=projector
        )
    else:
        log.warning("Step 3/4 skipped (no transcript table): events are flagged on length only")
    log.info("Step 4/4 - events of interest")
    final = highlight.run_highlight(
        step3 or step2,
        min_length=min_length,
        mode=mode,
        coding_only=step3 is not None,
        fill_color=fill_color,
    )
    log.info("Final table: %s", final)
    return PipelineResult(summary=step1.summary, genes=step2, hgvs=step3, final=final)
