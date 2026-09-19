# Examples

Everything here is **synthetic**: toy capture design, generated reads, fictitious gene
names (`DEMOA`, `DEMOB`) and transcripts (`NM_9000xx`). No patient data.

## `small_events/` - steps 0 to 4

| File | Content |
|---|---|
| `make_example_bams.py` (one level up) | writes `SAMPLE_A/B/C.bam` into `small_events/bam/` (needs samtools) |
| `panel_demo.bed` | 4 target regions of 1 kb |
| `tsv/` | output of step 0 on the three BAM files |
| `transcripts_demo.tsv` | transcript table for step 3 (fictitious transcripts: step 3 cannot map them) |

Simulated events (depth 200x, heterozygous):

| Event | SAMPLE_A | SAMPLE_B | SAMPLE_C | Expected result |
|---|---|---|---|---|
| deletion chr13:1,000,401-1,000,460 | ✓ | | | `occ_0`, highlighted |
| duplication chr17:2,000,501-2,000,580 | | ✓ | | `occ_0`, highlighted |
| deletion chr17:2,000,201-2,000,202 (2 bp) | | | ✓ | `occ_0`, too short |
| deletion chr17:2,003,301-2,003,340 | ✓ | ✓ | | `occ_1`, highlighted |
| deletion chr13:1,005,701-1,005,730 | ✓ | ✓ | ✓ | occurs in 2 other samples: removed |

```bash
python examples/make_example_bams.py
bash scripts/bam_to_depth_jumps.sh --bed examples/small_events/panel_demo.bed \
    --in-dir examples/small_events/bam --out-dir results/tsv
eyeoncnv run results/tsv --bed examples/small_events/panel_demo.bed --out-dir results
```
