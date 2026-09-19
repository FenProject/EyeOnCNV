from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest

from eyeoncnv.common import setup_logging

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
SMALL_EVENTS = EXAMPLES / "small_events"
SCRIPT = ROOT / "scripts" / "bam_to_depth_jumps.sh"


@pytest.fixture(autouse=True, scope="session")
def _logging():
    setup_logging(-1)


@pytest.fixture
def example_tsv_dir(tmp_path) -> Path:
    target = tmp_path / "tsv"
    shutil.copytree(SMALL_EVENTS / "tsv", target)
    return target


@pytest.fixture
def example_bed() -> Path:
    return SMALL_EVENTS / "panel_demo.bed"


def write_tsv(path: Path, rows) -> Path:
    path.write_text("".join(f"{c}\t{p}\t{d}\n" for c, p, d in rows))
    return path


class FakeProjector:
    """Stand-in for HgvsProjector: c. = (position % 1000) with a few non-coding forms."""

    def __init__(self, table=None, fail_on=()):
        self.table = table or {}
        self.fail_on = set(fail_on)
        self.calls = []

    def g_to_c(self, transcript, chrom, position):
        self.calls.append((transcript, chrom, position))
        if position in self.fail_on:
            raise ValueError(f"position {position} outside {transcript}")
        coord = self.table.get(position, str(position % 1000))
        return coord, f"{transcript}:c.{coord}="


def read_sheets(path) -> dict[str, pd.DataFrame]:
    return pd.read_excel(path, sheet_name=None)
