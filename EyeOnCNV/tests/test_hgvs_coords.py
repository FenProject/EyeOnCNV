from __future__ import annotations

import importlib.util

import pandas as pd
import pytest

from conftest import FakeProjector, read_sheets
from eyeoncnv import hgvs_coords as hc


@pytest.mark.parametrize(
    ("field", "expected"),
    [
        ("BRCA2-NM000059", ("BRCA2", "NM_000059")),
        ("BRCA2-NM_000059", ("BRCA2", "NM_000059")),
        ("BRCA2-NM_000059.4", ("BRCA2", "NM_000059")),
        (" BRCA1 - NM007294 ", ("BRCA1", "NM_007294")),
        ("HLA-A-NM_002116", ("HLA-A", "NM_002116")),
        ("RAD51C-NM058216;XRCC2-NM005431", ("RAD51C", "NM_058216")),
        ("CHEK2", ("CHEK2", None)),
        ("", (None, None)),
        (None, (None, None)),
        (float("nan"), (None, None)),
        ("not a gene!", (None, None)),
    ],
)
def test_parse_gene_transcript(field, expected):
    assert hc.parse_gene_transcript(field) == expected


def test_chrom_to_nc():
    assert hc.chrom_to_nc("chr13") == "NC_000013.10"
    assert hc.chrom_to_nc("13", "hg38") == "NC_000013.11"
    assert hc.chrom_to_nc("chrM", "GRCh37") == "NC_012920.1"
    assert hc.chrom_to_nc("X", "GRCh38") == "NC_000023.11"
    with pytest.raises(ValueError):
        hc.chrom_to_nc("chr99")
    with pytest.raises(ValueError):
        hc.normalize_assembly("mm10")


def test_load_transcript_map(tmp_path):
    table = pd.DataFrame(
        {
            "Genes": ["BRCA2", "brca1", "BRCA2", "TP53", "PALB2"],
            "NM": ["NM_000059.4", "NM007294.4", "NM_000059.3", None, "garbage"],
        }
    )
    path = tmp_path / "nm.xlsx"
    table.to_excel(path, index=False)
    by_nm, by_gene = hc.load_transcript_map(path)
    assert by_nm == {"NM_000059": "NM_000059.4", "NM_007294": "NM_007294.4"}
    assert by_gene == {"BRCA2": "NM_000059.4", "BRCA1": "NM_007294.4"}
    tsv = tmp_path / "nm.tsv"
    tsv.write_text("gene\ttranscript\nBRCA2\tNM_000059.4\n")
    assert hc.load_transcript_map(tsv)[1] == {"BRCA2": "NM_000059.4"}
    bad = tmp_path / "bad.csv"
    bad.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError):
        hc.load_transcript_map(bad)


def test_add_c_coordinates_with_errors():
    df = pd.DataFrame(
        {
            "chr": ["13", "13", "13", "17", "13"],
            "start": [1467, 2005, 3000, 4000, "x"],
            "end": [1470, 2005, 3001, 4010, 10],
            "gene": ["BRCA2-NM000059", "BRCA2-NM_000059", "UNKNOWN-NM999", "CHEK2", "BRCA2-NM000059"],
        }
    )
    projector = FakeProjector(table={1467: "467", 1470: "470+2"}, fail_on={3001})
    by_nm = {"NM_000059": "NM_000059.4"}
    by_gene = {"CHEK2": "NM_007194.4"}
    out = hc.add_c_coordinates(df, projector, by_nm, by_gene)
    versions = [None if pd.isna(v) else v for v in out["nm_version"]]
    assert versions[:4] == ["NM_000059.4", "NM_000059.4", None, "NM_007194.4"]
    assert out["nm_source"].tolist() == ["by_nm_base", "by_nm_base", "missing", "by_gene", "by_nm_base"]
    assert out.loc[0, "c_start"] == "467" and out.loc[0, "c_end"] == "470+2"
    assert out.loc[0, "hgvs_start"] == "NM_000059.4:c.467="
    assert out.loc[0, "hgvs_error"] == ""
    assert out.loc[2, "hgvs_error"] == "no versioned NM found (neither by NM nor by gene)"
    assert out.loc[3, "c_start"] == "0" and out.loc[3, "hgvs_error"] == ""
    assert out.loc[4, "hgvs_error"].startswith("start: position is not an integer")
    assert pd.isna(out.loc[4, "c_start"]) and out.loc[4, "c_end"] == "10"


def test_projection_failure_is_reported():
    df = pd.DataFrame({"chr": ["13"], "start": [3000], "end": [3001], "gene": ["BRCA2-NM000059"]})
    out = hc.add_c_coordinates(df, FakeProjector(fail_on={3001}), {"NM_000059": "NM_000059.4"}, {})
    assert out.loc[0, "c_start"] == "0"
    assert out.loc[0, "hgvs_error"] == "end: position 3001 outside NM_000059.4"


def test_run_hgvs_with_fake_projector(tmp_path):
    source = tmp_path / "events_genes.xlsx"
    pd.DataFrame({"chr": ["13"], "start": [1467], "end": [1470], "gene": ["BRCA2-NM000059"]}).to_excel(
        source, index=False, sheet_name="occ_0"
    )
    nm = tmp_path / "nm.tsv"
    nm.write_text("genes\tNM\nBRCA2\tNM_000059.4\n")
    output = hc.run_hgvs(source, nm, projector=FakeProjector(table={1467: "467"}))
    assert output.name == "events_genes_hgvs.xlsx"
    assert read_sheets(output)["occ_0"].loc[0, "hgvs_start"] == "NM_000059.4:c.467="


@pytest.mark.skipif(importlib.util.find_spec("hgvs") is not None, reason="hgvs is installed")
def test_missing_hgvs_gives_install_hint():
    with pytest.raises(ImportError, match="eyeoncnv\\[hgvs\\]"):
        hc.HgvsProjector()
