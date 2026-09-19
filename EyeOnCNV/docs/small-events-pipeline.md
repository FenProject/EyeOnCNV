# Small-event pipeline - method

This page describes what each step computes, with the default parameters. Coordinates
are 1-based unless stated otherwise.

- [Principle](#principle)
- [Step 0 - depth jumps](#step-0---depth-jumps-bam_to_depth_jumpssh)
- [Step 1 - events and recurrence](#step-1---events-and-recurrence-eyeoncnv-intervals)
- [Step 2 - genes](#step-2---genes-eyeoncnv-genes)
- [Step 3 - c. coordinates](#step-3---c-coordinates-eyeoncnv-hgvs)
- [Step 4 - events of interest](#step-4---events-of-interest-eyeoncnv-highlight)
- [Scope and limitations](#scope-and-limitations)
- [Parameter history](#parameter-history)

## Principle

In a capture panel, read depth is fairly smooth from one base to the next. A
heterozygous deletion halves it over the deleted bases (−50 %), a heterozygous
duplication multiplies it by 1.5 (+50 %). The edges of such an event are therefore two
*depth jumps* between consecutive bases, a few tens of bases apart - the steps a
biologist follows when reviewing a BAM by eye. The pipeline finds these jumps, pairs them
into events and keeps the events that are rare across samples.

The target is the window where bin-based callers usually miss an event, roughly 30 to
100 bp; the parameters below were tuned on events of that size.

Worked example (sample `SAMPLE_A` of the synthetic data, a 60-bp heterozygous deletion):

```
samtools depth                step 0 output              step 1 event
chr13  1000399  200
chr13  1000400  200    ->     chr13  1000400  200   \
chr13  1000401  100    ->     chr13  1000401  100    |    chr13  1000401  1000460  deletion
...    (100x)                                        |    (start = 2nd position,
chr13  1000460  100    ->     chr13  1000460  100    |     end = second-to-last)
chr13  1000461  200    ->     chr13  1000461  200   /
chr13  1000462  200
```

## Step 0 - depth jumps (`bam_to_depth_jumps.sh`)

For each BAM:

1. **Depth.** `samtools depth -b BED -q 0 -Q 0 -d 0 -g SECONDARY` gives the depth of every
   covered position of the design, with no base or mapping quality filter and counting
   secondary alignments (duplicates, QC-failed and unmapped reads stay excluded, which is
   the samtools default). `-d 0` removes the depth cap of old samtools versions.
2. **Jumps.** Positions with depth < `--min-depth` (30) are ignored. Each remaining position
   is compared with the previous remaining position of the same chromosome:

   relative change = |d₂ − d₁| / d₁ (infinite when d₁ = 0 and d₂ > 0)

   The threshold is `--th-low` (0.15) when both depths are ≥ `--depth-cut` (150), and
   `--th-high` (0.30) otherwise, to tolerate the larger random fluctuations of low depth.
   When the change reaches the threshold, **both** positions are printed (a position is
   never printed twice in a row).
3. **Isolated positions.** A printed position is kept only if the previous or next printed
   position of the same chromosome lies within `--win` (50) bp.

Output: `<sample>.tsv` (chromosome, position, depth; no header). The file is written under
a temporary name and renamed at the end, so an interrupted run leaves no truncated TSV.

| Option | Default | Environment variable |
|---|---|---|
| `--bed` | required | `BED` |
| `--in-dir` and/or BAM arguments | required | `IN_DIR` |
| `--out-dir` | required | `OUT_DIR` |
| `--min-depth` | 30 | `MIN_DEPTH` |
| `--th-low` | 0.15 | `TH_LOW` |
| `--th-high` | 0.30 | `TH_HIGH` |
| `--depth-cut` | 150 | `DEPTH_CUT` |
| `--win` | 50 | `WIN` |
| `--mode` | `ge` (keep changes ≥ threshold; `le` keeps the others) | `MODE` |
| `--samtools` | `samtools` | `SAMTOOLS` |

The awk code is POSIX: results are identical with gawk (Linux, Git Bash), mawk and the BSD
awk of macOS (checked by the tests).

## Step 1 - events and recurrence (`eyeoncnv intervals`)

**Events.** Per sample and chromosome, positions sorted by coordinate are split wherever
two consecutive positions are ≥ `--gap` (150) bp apart. For each group of *n* positions:

| n | start | end | variation |
|---|---|---|---|
| ≥ 3 | 2nd position | second-to-last position | `duplication` if depth(2nd) > depth(1st), else `deletion` |
| 2 | 1st position | 2nd position | same rule |
| 1 | the position | the position | empty |

Because step 0 prints the base before and the base after each jump, the 2nd and the
second-to-last positions are the first and last bases *inside* the event.

**Recurrence.** For each event, `occurrences` is the number of *other* samples having at
least one event on the same chromosome that overlaps it by ≥ `--min-overlap` (5 %) of
its length (and by at least 1 bp). Each other sample counts at most once. Chromosome names
are compared without the `chr` prefix.

The other samples are the TSV files of the run **and** the per-sample Excel files
already present in the history folder (by default the output folder). Files of samples
processed in the current run are not used as history, so re-running a batch does not
inflate the counts. A recurrent event is most likely a technical artefact (capture edge,
homopolymer, GC-rich target...) and its detection improves as samples accumulate.

**Outputs.**

- `<sample>_intervals_with_variation.xlsx`: every event of the sample, with `occurrences`;
- `summary_occurrences_0-<N>.xlsx`: events with `occurrences` ≤ `--max-occurrences` (1),
  one sheet per value (`occ_0`, `occ_1`), with a `sample` column.

## Step 2 - genes (`eyeoncnv genes`)

Each event `[start, end]` is converted to 0-based half-open `[start − 1, end)` and compared
with the BED regions `[start0, end0)`; the names (4th column) of all overlapping regions are
sorted, de-duplicated and joined with `;` in the `gene` column (empty outside the design).
BED header lines (`track`, `browser`, `#`) and gzipped files are accepted. Every sheet of
the input workbook is processed.

## Step 3 - c. coordinates (`eyeoncnv hgvs`)

1. **Transcript.** The first name of `gene` is parsed as `GENE-NMxxxxxx` (`NM` with or without
   underscore, version ignored); a plain gene symbol is accepted too. The versioned
   transcript comes from the `--nm-map` table (columns `genes` and `NM`): by NM number
   first (`nm_source = by_nm_base`), otherwise by gene symbol (`by_gene`), otherwise `missing`.
2. **Projection.** Start and end are written as identity variants on the chromosome RefSeq
   accession (`NC_000013.10:g.32890572=` for GRCh37) and mapped with the hgvs
   `AssemblyMapper` (splign alignments, UTA database). The results are stored as the
   coordinate (`c_start`: `467`, `467+5`, `-45`, `*23`...) and the full description
   (`hgvs_start`: `NM_000059.4:c.467=`), and the same for `end`.
3. Failures (position outside the transcript, unknown transcript, network error) are
   written in `hgvs_error`; the other rows are still processed.

`--assembly` accepts GRCh37/hg19 (default) and GRCh38/hg38.

## Step 4 - events of interest (`eyeoncnv highlight`)

`length = end − start + 1`. An event is flagged (`of_interest = True`, row filled in yellow)
when:

- `length ≥ --min-length` (3 bp), and
- the c. coordinate of its **start** contains none of `-` (5′ UTR or intron), `+` (intron)
  and `*` (3′ UTR), i.e. lies in the coding sequence.

`--mode any` accepts events whose start *or* end is coding, `--mode both` requires both.
`--no-coding-filter` flags on length only (used by `eyeoncnv run` when step 3 is skipped).
`--color` changes the fill colour.

## Scope and limitations

- **Event size.** The two breakpoints of an event must be less than `--gap` (150 bp) apart
  to be paired. A longer event appears as two 2-bp records at its edges, with opposite
  directions (e.g. a 300-bp deletion gives a `deletion` record at its start and a
  `duplication` record at its end); these are too short to be highlighted.
- **Homozygous deletions are not detected.** Positions under `--min-depth` are skipped, so
  depth is compared across the deleted segment and no jump is seen. The same applies to
  any event whose inner depth falls under 30x.
- **Targets.** Only positions inside the BED are seen. An event overlapping the edge of a
  capture region may lose one of its breakpoints.
- **Artefacts.** Capture edges, repeats or GC-rich regions also create depth jumps. The
  recurrence count removes them only if enough samples are in the history.
- **Strand.** `c_start` is the c. coordinate of the *genomic* start; on minus-strand genes it
  is the 3′ end of the event in transcript orientation. The default rule tests that position.
- **Several regions.** When an event overlaps several BED regions, step 3 uses the first name.
- **Validation.** Parameters were tuned on one capture design. Check sensitivity and
  specificity on known positive samples before using the pipeline on another panel.

## Parameter history

The lab notes (`Read_me_cnv_petit_evenement.docx`) tried several filters before the
final folder script. They can be reproduced with the options of step 0:

| Variant in the notes | Equivalent options |
|---|---|
| depth ≥ 30 only ("like Alamut") | `samtools depth ... \| awk '$3>=30'` (no jump filter) |
| jump ≥ 15 % | `--th-low 0.15 --th-high 0.15 --win 1000000000` |
| jump ≥ 15 % + 50-bp neighbour filter | `--th-low 0.15 --th-high 0.15 --win 50` |
| 15 %, or 50 % under 100x, + neighbour filter | `--th-high 0.50 --depth-cut 100` |
| **final script (default)**: 15 %, or 30 % under 150x, + neighbour filter | *(defaults)* |
