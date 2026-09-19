#!/usr/bin/env python3
"""Create small synthetic BAM files to try step 0 (bam_to_depth_jumps.sh).

The reads are fake (poly-A sequences on a toy capture design); no real sample
is involved. Each sample has a uniform depth of 200x inside the BED, and
heterozygous events are simulated on one haplotype:

* deletions: reads carry a CIGAR ``D`` over the deleted bases (depth 200 -> 100);
* duplications: extra soft-clipped reads cover the duplicated bases (200 -> 300).

Requires ``samtools`` on the PATH. Usage::

    python examples/make_example_bams.py [OUTPUT_DIR]

OUTPUT_DIR defaults to ``examples/small_events/bam``.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

READ_LENGTH = 100
HERE = Path(__file__).resolve().parent

# Contig lengths only matter for the BAM header.
CONTIGS = {"chr13": 115_169_878, "chr17": 81_195_210}

# Toy capture design (BED, 0-based half-open) - names follow the GENE-NM format.
BED_REGIONS = [
    ("chr13", 1_000_000, 1_001_000, "DEMOA-NM_900001"),
    ("chr13", 1_005_000, 1_006_000, "DEMOA-NM_900001"),
    ("chr17", 2_000_000, 2_001_000, "DEMOB-NM_900002"),
    ("chr17", 2_003_000, 2_004_000, "DEMOB-NM_900002"),
]

# Simulated events: (kind, chromosome, first base, last base), 1-based inclusive.
ARTEFACT = ("del", "chr13", 1_005_701, 1_005_730)  # present in every sample
SHARED_DEL = ("del", "chr17", 2_003_301, 2_003_340)  # present in 2 samples
SAMPLES = {
    "SAMPLE_A": [("del", "chr13", 1_000_401, 1_000_460), SHARED_DEL, ARTEFACT],
    "SAMPLE_B": [("dup", "chr17", 2_000_501, 2_000_580), SHARED_DEL, ARTEFACT],
    "SAMPLE_C": [("del", "chr17", 2_000_201, 2_000_202), ARTEFACT],
}


def _sam_line(name: str, chrom: str, pos: int, cigar: str, query_length: int) -> str:
    return "\t".join([name, "0", chrom, str(pos), "60", cigar, "*", "0", "0", "A" * query_length, "*"])


def _haplotype2_read(chrom: str, p: int, events: list) -> tuple[int, str] | None:
    """Read starting at p on the haplotype carrying the deletions (None = dropped)."""
    end = p + READ_LENGTH - 1
    for kind, ev_chrom, start, stop in events:
        if kind != "del" or ev_chrom != chrom:
            continue
        if start <= p <= stop:
            return None  # would start inside the deleted sequence
        if p < start <= end:
            left = start - p
            return p, f"{left}M{stop - start + 1}D{READ_LENGTH - left}M"
    return p, f"{READ_LENGTH}M"


def sample_sam_lines(sample: str, events: list) -> list[str]:
    lines = []
    n = 0
    for chrom, start0, end0, _name in BED_REGIONS:
        first, last = start0 + 1, end0  # 1-based inclusive
        for p in range(first - READ_LENGTH + 1, last + 1):
            n += 1
            lines.append(_sam_line(f"{sample}_r{n}", chrom, p, f"{READ_LENGTH}M", READ_LENGTH))
            read = _haplotype2_read(chrom, p, events)
            if read is not None:
                n += 1
                lines.append(_sam_line(f"{sample}_r{n}", chrom, read[0], read[1], READ_LENGTH))
    for kind, chrom, start, stop in events:
        if kind != "dup":
            continue
        for p in range(start - READ_LENGTH + 1, stop + 1):
            left = max(0, start - p)
            right = max(0, p + READ_LENGTH - 1 - stop)
            matched = READ_LENGTH - left - right
            cigar = (f"{left}S" if left else "") + f"{matched}M" + (f"{right}S" if right else "")
            n += 1
            lines.append(_sam_line(f"{sample}_r{n}", chrom, max(p, start), cigar, READ_LENGTH))
    return lines


def main(argv: list[str]) -> int:
    out_dir = Path(argv[1]) if len(argv) > 1 else HERE / "small_events" / "bam"
    samtools = shutil.which("samtools")
    if samtools is None:
        print("[ERROR] samtools is required to build the example BAM files", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)
    header = ["@HD\tVN:1.6\tSO:unsorted"]
    header += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in CONTIGS.items()]
    for sample, events in SAMPLES.items():
        sam = out_dir / f"{sample}.sam"
        bam = out_dir / f"{sample}.bam"
        sam.write_text("\n".join(header + sample_sam_lines(sample, events)) + "\n")
        subprocess.run([samtools, "sort", "-o", str(bam), str(sam)], check=True)
        subprocess.run([samtools, "index", str(bam)], check=True)
        sam.unlink()
        print(f"[INFO] {bam}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
