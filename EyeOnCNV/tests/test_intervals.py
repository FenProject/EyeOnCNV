from __future__ import annotations

import pandas as pd
import pytest

from conftest import read_sheets, write_tsv
from eyeoncnv import intervals


def events(rows):
    df = pd.DataFrame(rows, columns=["chr", "position", "depth"])
    return intervals.build_intervals(df)


def test_interval_rules_by_number_of_positions():
    out = events(
        [
            ("chr1", 100, 200),
            ("chr1", 101, 100),
            ("chr1", 160, 100),
            ("chr1", 161, 200),  # 4
            ("chr1", 1000, 200),
            ("chr1", 1001, 300),
            ("chr1", 1010, 300),  # 3
            ("chr1", 5000, 300),
            ("chr1", 5001, 200),  # 2
            ("chr1", 9000, 100),  # 1
        ]
    )
    assert out.to_dict("records") == [
        {"chr": "chr1", "start": 101, "end": 160, "variation": "deletion"},
        {"chr": "chr1", "start": 1001, "end": 1001, "variation": "duplication"},
        {"chr": "chr1", "start": 5000, "end": 5001, "variation": "deletion"},
        {"chr": "chr1", "start": 9000, "end": 9000, "variation": ""},
    ]


def test_gap_is_strict_and_input_is_sorted():
    out = events([("chr2", 300, 50), ("chr2", 1, 100), ("chr2", 150, 90), ("chr2", 449, 80)])
    # 1 -> 150 and 300 -> 449 are 149 bp apart (same event); 150 -> 300 is 150 bp (new event)
    assert out[["start", "end"]].values.tolist() == [[1, 150], [300, 449]]


def test_chromosomes_are_kept_apart():
    out = events([("chr1", 10, 100), ("chr2", 11, 50), ("chr1", 12, 60)])
    assert out["chr"].tolist() == ["chr1", "chr2"]


def test_read_empty_tsv(tmp_path):
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    df = intervals.read_depth_tsv(empty)
    assert df.empty and list(df.columns) == ["chr", "position", "depth"]
    assert intervals.build_intervals(df).empty


def test_count_recurrence_threshold_and_naming():
    query = pd.DataFrame({"chr": ["chr1", "chr1", "chrX"], "start": [1, 1001, 50], "end": [100, 1100, 60]})
    others = {
        "exact_5_percent": pd.DataFrame({"chr": ["1"], "start": [96], "end": [200]}),  # 5 bp of 100
        "below_threshold": pd.DataFrame({"chr": ["chr1"], "start": [1097], "end": [2000]}),  # 4 bp
        "two_hits_count_once": pd.DataFrame({"chr": ["X", "X"], "start": [50, 55], "end": [52, 60]}),
    }
    counts = intervals.count_recurrence(query, others, min_overlap=0.05)
    assert counts.tolist() == [1, 0, 1]
    assert intervals.count_recurrence(query, others, min_overlap=0.0).tolist() == [1, 1, 1]


def test_run_on_examples(example_tsv_dir, tmp_path):
    out = tmp_path / "results"
    result = intervals.run_intervals([example_tsv_dir], out)
    assert set(result.per_sample) == {"SAMPLE_A", "SAMPLE_B", "SAMPLE_C"}
    sample_a = pd.read_excel(result.per_sample["SAMPLE_A"])
    assert sample_a[["start", "end", "variation", "occurrences"]].values.tolist() == [
        [1000401, 1000460, "deletion", 0],
        [1005701, 1005730, "deletion", 2],
        [2003301, 2003340, "deletion", 1],
    ]
    sheets = read_sheets(result.summary)
    assert result.summary.name == "summary_occurrences_0-1.xlsx"
    assert list(sheets) == ["occ_0", "occ_1"]
    assert sheets["occ_0"][["sample", "start", "variation"]].values.tolist() == [
        ["SAMPLE_A", 1000401, "deletion"],
        ["SAMPLE_B", 2000501, "duplication"],
        ["SAMPLE_C", 2000201, "deletion"],
    ]
    assert sheets["occ_1"]["sample"].tolist() == ["SAMPLE_A", "SAMPLE_B"]


def test_rerun_does_not_count_a_sample_against_itself(example_tsv_dir, tmp_path):
    out = tmp_path / "results"
    first = intervals.run_intervals([example_tsv_dir], out)
    before = {s: pd.read_excel(p)["occurrences"].tolist() for s, p in first.per_sample.items()}
    second = intervals.run_intervals([example_tsv_dir], out)
    after = {s: pd.read_excel(p)["occurrences"].tolist() for s, p in second.per_sample.items()}
    assert before == after
    assert second.history_samples == 0


def test_history_from_previous_runs(example_tsv_dir, tmp_path):
    previous = tmp_path / "batch1"
    intervals.run_intervals([example_tsv_dir / "SAMPLE_A.tsv", example_tsv_dir / "SAMPLE_B.tsv"], previous)
    new_batch = tmp_path / "batch2"
    result = intervals.run_intervals([example_tsv_dir / "SAMPLE_C.tsv"], new_batch, history_dirs=[previous])
    assert result.history_samples == 2
    sample_c = pd.read_excel(result.per_sample["SAMPLE_C"])
    # artefact shared with A and B (history) -> 2; SAMPLE_C private event -> 0
    assert sorted(sample_c["occurrences"].tolist()) == [0, 2]
    isolated = intervals.run_intervals([example_tsv_dir / "SAMPLE_C.tsv"], new_batch, use_history=False)
    assert pd.read_excel(isolated.per_sample["SAMPLE_C"])["occurrences"].tolist() == [0, 0]


def test_max_occurrences_controls_summary_sheets(example_tsv_dir, tmp_path):
    result = intervals.run_intervals([example_tsv_dir], tmp_path / "r", max_occurrences=2)
    sheets = read_sheets(result.summary)
    assert list(sheets) == ["occ_0", "occ_1", "occ_2"]
    assert len(sheets["occ_2"]) == 3


def test_bad_file_is_reported_but_others_are_processed(tmp_path):
    good = write_tsv(tmp_path / "good.tsv", [("chr1", 10, 100), ("chr1", 11, 50)])
    bad = tmp_path / "bad.tsv"
    bad.write_text("only_one_column\n")
    result = intervals.run_intervals([good, bad], tmp_path / "out")
    assert "bad" in result.errors and "good" in result.per_sample


def test_no_tsv_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        intervals.run_intervals([tmp_path], tmp_path / "out")
