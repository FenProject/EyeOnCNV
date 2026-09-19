# File formats

All Excel files written by the tool have a bold, frozen header row and filters enabled.
Positions are 1-based unless stated otherwise.

## Inputs

### Capture BED (`--bed`)

Tab-separated, 0-based half-open, plain or `.gz`. Columns: chromosome, start, end,
name. `track`, `browser` and `#` lines are skipped. For step 3, the name should follow the
`GENE-NMxxxxxx` convention (e.g. `BRCA2-NM000059` or `BRCA2-NM_000059`); a plain gene symbol
also works if it is in the transcript table.

```
chr13	32890558	32890664	BRCA2-NM000059
```

### Transcript table (`--nm-map`)

Excel, CSV or TSV with a gene column (`genes`, `gene` or `symbol`) and a transcript column
(`NM`, `transcript` or `refseq`). The first row of a gene or NM number wins.

| genes | NM |
|---|---|
| BRCA2 | NM_000059.4 |
| BRCA1 | NM_007294.4 |

## Step 0 - `<sample>.tsv`

No header, tab-separated: `chromosome`, `position`, `depth`. The sample name is the BAM
file name without `.bam`.

## Step 1

`<sample>_intervals_with_variation.xlsx` (sheet `intervals`) - every event of a sample:

| Column | Content |
|---|---|
| `chr` | chromosome as in the BAM |
| `start`, `end` | first and last base of the event |
| `variation` | `deletion`, `duplication` or empty (single position) |
| `occurrences` | number of other samples (run + history) with an overlapping event |

`summary_occurrences_0-1.xlsx` - events with `occurrences` ≤ 1, sheets `occ_0` and `occ_1`,
columns `sample, chr, start, end, variation, occurrences`.

## Steps 2 to 4

Each step keeps the sheets and columns of its input and adds:

| Step | File suffix | Added columns |
|---|---|---|
| 2 | `_genes.xlsx` | `gene` - names of the overlapping BED regions, `;`-separated |
| 3 | `_hgvs.xlsx` | `gene_symbol`, `nm_input` (NM parsed from `gene`), `nm_version` (versioned transcript used), `nm_source` (`by_nm_base`, `by_gene`, `missing`), `c_start`, `hgvs_start`, `c_end`, `hgvs_end`, `hgvs_error` |
| 4 | `_highlight.xlsx` | `length` (after `end`), `of_interest` (True/False); flagged rows are filled |

With `eyeoncnv run` and the default names, the final file is
`summary_occurrences_0-1_genes_hgvs_highlight.xlsx` (or `..._genes_highlight.xlsx` without
`--nm-map`).
