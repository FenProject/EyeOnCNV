# Changes from the original notebooks

This repository replaces two working files, which are kept unchanged outside the
repository: `Read_me_cnv_petit_evenement.docx` (command lines) and `CNV_eye_of_medah.ipynb`.
The notebook `CNV_kit_complement.ipynb` (dbVar annotation and plots of CNVkit results) is a
separate tool and is not part of this repository.

## Where each part went

| Original | Now |
|---|---|
| docx, section "Pour tout le dossier" (folder script) | `scripts/bam_to_depth_jumps.sh` |
| docx, earlier single-BAM variants | options of the same script ([parameter history](small-events-pipeline.md#parameter-history)) |
| `CNV_eye_of_medah.ipynb`, cell 1 - TSV to intervals, occurrences, summary | `eyeoncnv intervals` (`intervals.py`) |
| cell 2 - genes from the BED | `eyeoncnv genes` (`genes.py`) |
| cell 3 - g. to c. with hgvs | `eyeoncnv hgvs` (`hgvs_coords.py`) |
| cell 4 - length and highlighting | `eyeoncnv highlight` (`highlight.py`) |
| *(new)* | `eyeoncnv run` - steps 1 to 4 (`pipeline.py`) |

## General changes

- Hard-coded Windows paths are replaced by command-line arguments, and input files are
  given explicitly. Files are no longer picked as "the most recent file starting with a
  prefix", which had made step 3 read its own output (`..._with_c_coords__with_c_coords.xlsx`).
- A single code base runs on macOS, Linux and Windows; the steps still exchange Excel
  files, so they can run on different computers (e.g. step 3 on a Mac).
- Messages, column names and values are in English; Excel outputs are formatted.
- Default parameters are those of the notebook and of the lab notes.

## Results checked against the originals

The original cells were run next to the new code on synthetic data:

| Part | Result |
|---|---|
| Step 0 | byte-identical TSV files on toy BAM files and on 40 random depth streams, with gawk, mawk and BSD awk |
| Step 1 | identical events and occurrence counts (25 samples, 656 events) on a first run |
| Step 2 | identical `gene` values on 3,000 random events |
| Step 3 | identical transcripts and c. coordinates (mocked hgvs), except the fixes 5 and 6 below |
| Step 4 | identical flags on 2,000 rows |

## Bugs fixed

1. **Step 1, re-runs.** The per-sample files of a previous run of the *same* samples were
   counted as other samples: every re-run in the same folder increased the occurrences.
   Those files are now left out of the history.
2. **Step 2 read only the first sheet** (`occ_0`) of the summary, so `occ_1` never reached
   steps 2-4. Every sheet is processed.
3. **Step 2 BED reading** failed on `track`/`browser` header lines (as in Agilent
   `*_Covered.bed` files); these lines and `.bed.gz` files are now accepted. The helper
   columns `_chrom_norm`, `_start0`, `_end0` are no longer written in the output.
4. **Step 3, end position.** Only the start was projected; `c_end`/`hgvs_end` are added.
   The highlight rule still tests the start by default (`--mode` for other choices).
5. **Step 3, NM parsing.** `GENE-NM_000059` became `NM__000059` (found only through the
   gene symbol), and `GENE-NM_000059.4` or a lowercase `nm` were not parsed at all.
6. **Step 3, transcript table.** Empty `NM` cells were read as the text `nan` with pandas 2
   and could be used as a transcript; they are ignored.
7. **Step 0.** A failed run left a truncated TSV that step 1 would read; the file is now
   written under a temporary name. An index named `sample.bai` or `.csi` is recognised.
8. The comment of cell 4 said "orange" while the code filled rows in yellow: yellow is
   kept, `--color` changes it.

## Minor behaviour changes

- `--min-overlap 0` requires at least 1 bp of overlap (before: any event of the same
  chromosome counted). The default 5 % is unaffected.
- Chromosome names are compared without the `chr` prefix in the recurrence count.
- An empty TSV is a sample without events (before: an error).
- A BED name made of a gene symbol only is resolved through the transcript table.

## Renamed files and columns

| Before | Now |
|---|---|
| `resume_occurrences_0_1.xlsx` | `summary_occurrences_0-1.xlsx` |
| `..._with_gene.xlsx`, `...__with_c_coords.xlsx`, `...__with_diff_and_highlight.xlsx` | `..._genes.xlsx`, `..._hgvs.xlsx`, `..._highlight.xlsx` |
| `délétion` | `deletion` |
| `Gene` | `gene` |
| `nm_base_input`, `nm_version_used` | `nm_input`, `nm_version` (`nm_source` unchanged) |
| `c_coord`, `c_hgvs_eq`, `error` | `c_start`, `hgvs_start`, `hgvs_error` (+ `c_end`, `hgvs_end`) |
| `diff_end_start` | `length` (+ `of_interest`) |

Unchanged on purpose: `<sample>_intervals_with_variation.xlsx` (per-sample files of
earlier runs keep counting in the recurrence), the sheets `occ_0`/`occ_1`, and the
`c_hgvs_eq` column of old files, which `eyeoncnv highlight` still understands.
