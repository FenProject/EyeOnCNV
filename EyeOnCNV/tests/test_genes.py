from __future__ import annotations

import gzip

import pandas as pd

from conftest import read_sheets
from eyeoncnv import genes


def test_read_bed_skips_headers_and_handles_gzip(tmp_path):
    content = (
        'browser position chr13:1-100\ntrack name="design"\n# comment\n\n'
        "chr13\t100\t200\tBRCA2-NM000059\n"
        "13\t300\t400\n"
        "chr13\t500\t500\tempty\n"
        "chr17 10 20 SPACES-NM1\n"
    )
    plain = tmp_path / "design.bed"
    plain.write_text(content)
    gz = tmp_path / "design.bed.gz"
    with gzip.open(gz, "wt") as handle:
        handle.write(content)
    for path in (plain, gz):
        bed = genes.read_bed(path)
        assert bed["name"].tolist() == ["BRCA2-NM000059", "Unknown", "SPACES-NM1"]
        assert bed["chrom_key"].tolist() == ["13", "13", "17"]


def test_coordinate_conventions():
    bed = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start0": [99, 100, 150],
            "end0": [100, 200, 160],
            "name": ["B", "A", "A"],
            "chrom_key": ["1", "1", "1"],
        }
    )
    df = pd.DataFrame(
        {"chr": ["1", "chr1", "chr1", "chr2"], "start": [100, 101, 90, 100], "end": [100, 155, 99, 200]}
    )
    out = genes.annotate_genes(df, bed)
    # [100,100] (1-based) is base 99 (0-based): overlaps [99,100) only
    assert out["gene"].tolist() == ["B", "A", "", ""]


def test_run_genes_keeps_every_sheet(tmp_path, example_bed):
    source = tmp_path / "summary.xlsx"
    with pd.ExcelWriter(source) as writer:
        pd.DataFrame({"sample": ["S"], "chr": ["chr13"], "start": [1000401], "end": [1000460]}).to_excel(
            writer, sheet_name="occ_0", index=False
        )
        pd.DataFrame(columns=["sample", "chr", "start", "end"]).to_excel(
            writer, sheet_name="occ_1", index=False
        )
    output = genes.run_genes(source, example_bed)
    assert output.name == "summary_genes.xlsx"
    sheets = read_sheets(output)
    assert list(sheets) == ["occ_0", "occ_1"]
    assert sheets["occ_0"]["gene"].tolist() == ["DEMOA-NM_900001"]
    assert sheets["occ_1"].empty and "gene" in sheets["occ_1"].columns
