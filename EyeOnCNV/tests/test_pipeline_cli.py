from __future__ import annotations

import os
import subprocess
import sys

import pandas as pd

from conftest import ROOT, SMALL_EVENTS, FakeProjector, read_sheets
from eyeoncnv import __version__, cli
from eyeoncnv.pipeline import run_pipeline


def test_cli_run_without_transcripts(example_tsv_dir, example_bed, tmp_path):
    out = tmp_path / "results"
    code = cli.main(["-q", "run", str(example_tsv_dir), "-o", str(out), "--bed", str(example_bed)])
    assert code == 0
    final = out / "summary_occurrences_0-1_genes_highlight.xlsx"
    sheets = read_sheets(final)
    occ0 = sheets["occ_0"]
    assert occ0["gene"].tolist() == ["DEMOA-NM_900001", "DEMOB-NM_900002", "DEMOB-NM_900002"]
    assert occ0["of_interest"].tolist() == [True, True, False]  # the 2 bp event is too short


def test_pipeline_with_projector(example_tsv_dir, example_bed, tmp_path):
    projector = FakeProjector(table={1000401: "120", 2000501: "300+4", 2003301: "*40"})
    result = run_pipeline(
        [example_tsv_dir],
        tmp_path / "r",
        example_bed,
        SMALL_EVENTS / "transcripts_demo.tsv",
        projector=projector,
    )
    assert result.hgvs is not None and result.final.name.endswith("_genes_hgvs_highlight.xlsx")
    occ0 = pd.read_excel(result.final, sheet_name="occ_0")
    assert occ0["nm_version"].tolist() == ["NM_900001.1", "NM_900002.1", "NM_900002.1"]
    assert occ0["of_interest"].tolist() == [True, False, False]


def test_verbosity_flags_after_command(example_tsv_dir, tmp_path):
    assert cli.main(["intervals", str(example_tsv_dir), "-o", str(tmp_path / "o"), "-q"]) == 0


def test_errors_return_1(tmp_path):
    assert cli.main(["-q", "genes", str(tmp_path / "missing.xlsx"), "--bed", str(tmp_path / "x.bed")]) == 1


def test_module_entry_point():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT / "src"), env.get("PYTHONPATH")]))
    result = subprocess.run(
        [sys.executable, "-m", "eyeoncnv", "--version"], capture_output=True, text=True, env=env
    )
    assert result.returncode == 0 and __version__ in result.stdout
