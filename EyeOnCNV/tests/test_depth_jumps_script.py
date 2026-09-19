"""Step 0 shell script, checked against a Python re-implementation of its filters."""

from __future__ import annotations

import os
import random
import shutil
import subprocess
import sys

import pytest

from conftest import ROOT, SCRIPT, SMALL_EVENTS

BASH = os.environ.get("EYEONCNV_BASH", "bash")
AWKS = [a for a in ("awk", "gawk", "mawk", "original-awk", "nawk") if shutil.which(a)]

pytestmark = pytest.mark.skipif(shutil.which(BASH) is None, reason="bash not available")


def reference_filter(lines, min_depth=30, th_low=0.15, th_high=0.30, depth_cut=150, win=50, mode="ge"):
    """What the two awk programs of the script do, written in Python."""
    step1, last = [], None
    prev = None  # (chrom, depth, line)
    for line in lines:
        chrom, _pos, depth = line.split("\t")
        depth = float(depth)
        if depth < min_depth:
            continue
        if prev is not None and chrom == prev[0]:
            a, b = prev[1], depth
            rc = (0 if b == 0 else 1e99) if a == 0 else abs(b - a) / a
            th = th_high if (a < depth_cut or b < depth_cut) else th_low
            if (rc >= th) if mode == "ge" else (rc <= th):
                if prev[2] != last:
                    step1.append(prev[2])
                    last = prev[2]
                if line != last:
                    step1.append(line)
                    last = line
        else:
            last = ""
        prev = (chrom, depth, line)

    out = []
    for i, line in enumerate(step1):
        chrom, pos = line.split("\t")[0], int(line.split("\t")[1])
        near = False
        if i > 0:
            c, p = step1[i - 1].split("\t")[0], int(step1[i - 1].split("\t")[1])
            near |= c == chrom and pos - p <= win
        if i + 1 < len(step1):
            c, p = step1[i + 1].split("\t")[0], int(step1[i + 1].split("\t")[1])
            near |= c == chrom and p - pos <= win
        if near:
            out.append(line)
    return out


def random_depth_stream(seed: int) -> list[str]:
    rng = random.Random(seed)
    lines = []
    for chrom in ("chr1", "chr2", "chrX")[: rng.randint(1, 3)]:
        pos = rng.randint(1, 1000)
        depth = rng.choice([0, 20, 30, 80, 150, 400])
        for _ in range(rng.randint(0, 400)):
            pos += rng.choice([1, 1, 1, 2, 30, 50, 51, 400])
            r = rng.random()
            if r < 0.08:
                depth = rng.choice([0, 10, 29, 30, 60, 100, 149, 150, 200, 500])
            elif r < 0.5:
                depth = max(0, int(depth * rng.uniform(0.8, 1.25)))
            lines.append(f"{chrom}\t{pos}\t{depth}")
    return lines


@pytest.fixture
def fake_samtools(tmp_path):
    """A `samtools` that prints <FAKE_DEPTH_DIR>/<bam name>.depth for `depth`."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "samtools"
    exe.write_text(
        "#!/usr/bin/env bash\n"
        'case "$1" in\n'
        "  index) exit 0 ;;\n"
        '  depth) bam="${@: -1}"; cat "$FAKE_DEPTH_DIR/$(basename "$bam" .bam).depth" ;;\n'
        "  *) exit 1 ;;\n"
        "esac\n"
    )
    exe.chmod(0o755)
    return exe


def run_script(args, env_extra=None, cwd=None):
    env = dict(os.environ)
    env.update(env_extra or {})
    return subprocess.run(
        [BASH, str(SCRIPT), *map(str, args)], capture_output=True, text=True, env=env, cwd=cwd
    )


def make_inputs(tmp_path, streams):
    bams, depth_dir = tmp_path / "bams", tmp_path / "depth"
    bams.mkdir()
    depth_dir.mkdir()
    for name, lines in streams.items():
        (bams / f"{name}.bam").write_bytes(b"")
        (bams / f"{name}.bam.bai").write_bytes(b"")
        (depth_dir / f"{name}.depth").write_text("".join(f"{line}\n" for line in lines))
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t0\t10\n")
    return bams, depth_dir, bed


@pytest.mark.parametrize("awk", AWKS)
def test_matches_reference_on_random_streams(tmp_path, fake_samtools, awk):
    streams = {f"S{i}": random_depth_stream(i) for i in range(12)}
    bams, depth_dir, bed = make_inputs(tmp_path, streams)
    out = tmp_path / "out"
    result = run_script(
        ["--bed", bed, "--in-dir", bams, "--out-dir", out, "--samtools", fake_samtools],
        {"FAKE_DEPTH_DIR": str(depth_dir), "AWK": awk},
    )
    assert result.returncode == 0, result.stderr
    for name, lines in streams.items():
        assert (out / f"{name}.tsv").read_text().splitlines() == reference_filter(lines)


def test_custom_parameters(tmp_path, fake_samtools):
    streams = {f"S{i}": random_depth_stream(100 + i) for i in range(6)}
    bams, depth_dir, bed = make_inputs(tmp_path, streams)
    out = tmp_path / "out"
    params = dict(min_depth=20, th_low=0.1, th_high=0.5, depth_cut=100, win=10, mode="ge")
    result = run_script(
        [
            "-b",
            bed,
            "-i",
            bams,
            "-o",
            out,
            "--samtools",
            fake_samtools,
            "--min-depth",
            "20",
            "--th-low",
            "0.1",
            "--th-high",
            "0.5",
            "--depth-cut",
            "100",
            "--win",
            "10",
        ],
        {"FAKE_DEPTH_DIR": str(depth_dir)},
    )
    assert result.returncode == 0, result.stderr
    for name, lines in streams.items():
        assert (out / f"{name}.tsv").read_text().splitlines() == reference_filter(lines, **params)


def test_environment_variables_and_positional_bams(tmp_path, fake_samtools):
    streams = {"A": random_depth_stream(1), "B": random_depth_stream(2)}
    bams, depth_dir, bed = make_inputs(tmp_path, streams)
    out = tmp_path / "out"
    env = {
        "FAKE_DEPTH_DIR": str(depth_dir),
        "BED": str(bed),
        "OUT_DIR": str(out),
        "SAMTOOLS": str(fake_samtools),
    }
    result = run_script([bams / "B.bam"], env)
    assert result.returncode == 0, result.stderr
    assert sorted(p.name for p in out.iterdir()) == ["B.tsv"]


def test_help_and_errors(tmp_path):
    assert run_script(["--help"]).returncode == 0
    missing_bed = run_script(["--out-dir", tmp_path])
    assert missing_bed.returncode == 1 and "--bed" in missing_bed.stderr
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t0\t10\n")
    bad_mode = run_script(["-b", bed, "-o", tmp_path, "--mode", "xx"])
    assert bad_mode.returncode == 1 and "--mode" in bad_mode.stderr
    unknown = run_script(["--nope"])
    assert unknown.returncode == 1 and "Unknown option" in unknown.stderr


def test_no_bam_is_not_an_error(tmp_path, fake_samtools):
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t0\t10\n")
    (tmp_path / "empty").mkdir()
    result = run_script(
        ["-b", bed, "-i", tmp_path / "empty", "-o", tmp_path / "out", "--samtools", fake_samtools]
    )
    assert result.returncode == 0
    assert "No .bam file" in result.stderr


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not installed")
def test_example_bams_reproduce_committed_tsv(tmp_path):
    bam_dir = tmp_path / "bam"
    subprocess.run(
        [sys.executable, str(ROOT / "examples" / "make_example_bams.py"), str(bam_dir)],
        check=True,
        capture_output=True,
    )
    out = tmp_path / "tsv"
    result = run_script(["--bed", SMALL_EVENTS / "panel_demo.bed", "--in-dir", bam_dir, "--out-dir", out])
    assert result.returncode == 0, result.stderr
    for expected in sorted((SMALL_EVENTS / "tsv").glob("*.tsv")):
        assert (out / expected.name).read_text() == expected.read_text()
