"""Command-line interface: ``eyeoncnv <command> ...`` (or ``python -m eyeoncnv``)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from eyeoncnv import __version__
from eyeoncnv.common import log, setup_logging

EPILOG = "Step 0 (BAM -> TSV) is the shell script scripts/bam_to_depth_jumps.sh."


def _add_intervals_options(parser: argparse.ArgumentParser) -> None:
    from eyeoncnv import intervals as iv

    group = parser.add_argument_group("step 1 - events and recurrence")
    group.add_argument(
        "--gap",
        type=int,
        default=iv.DEFAULT_GAP,
        help="positions closer than this (bp) belong to the same event [%(default)s]",
    )
    group.add_argument(
        "--min-overlap",
        type=float,
        default=iv.DEFAULT_MIN_OVERLAP,
        help="fraction of an event another sample must cover to count as a recurrence [%(default)s]",
    )
    group.add_argument(
        "--max-occurrences",
        type=int,
        default=iv.DEFAULT_MAX_OCCURRENCES,
        help="keep events seen in at most this many other samples in the summary [%(default)s]",
    )
    group.add_argument(
        "--history-dir",
        type=Path,
        action="append",
        metavar="DIR",
        help="folder of *_intervals_with_variation.xlsx from previous runs, counted in "
        "the recurrence (repeatable; default: the output folder)",
    )
    group.add_argument(
        "--no-history",
        action="store_true",
        help="count recurrence among the input samples only",
    )


def _add_hgvs_options(parser: argparse.ArgumentParser, required: bool) -> None:
    group = parser.add_argument_group("step 3 - c. coordinates")
    group.add_argument(
        "--nm-map",
        type=Path,
        required=required,
        metavar="TABLE",
        help="transcript table with columns 'genes' and 'NM' (xlsx/csv/tsv)"
        + ("" if required else "; without it step 3 is skipped"),
    )
    group.add_argument("--assembly", default="GRCh37", help="GRCh37/hg19 or GRCh38/hg38 [%(default)s]")
    group.add_argument(
        "--uta-url",
        metavar="URL",
        help="UTA database URL (default: $UTA_DB_URL or the public biocommons server)",
    )


def _add_highlight_options(parser: argparse.ArgumentParser, standalone: bool) -> None:
    from eyeoncnv import highlight as hl

    group = parser.add_argument_group("step 4 - events of interest")
    group.add_argument(
        "--min-length",
        type=int,
        default=hl.DEFAULT_MIN_LENGTH,
        help="minimum event length in bp [%(default)s]",
    )
    group.add_argument(
        "--mode",
        choices=hl.MODES,
        default=hl.DEFAULT_MODE,
        help="which c. coordinate must be coding: start (validated rule), any, both [%(default)s]",
    )
    group.add_argument(
        "--color",
        default=hl.DEFAULT_FILL,
        help="RGB fill colour of the flagged rows [%(default)s]",
    )
    if standalone:
        group.add_argument(
            "--no-coding-filter",
            action="store_true",
            help="flag on length only (for tables without c. coordinates)",
        )


def _add_verbosity(parser: argparse.ArgumentParser, suppress: bool) -> None:
    """-v/-q accepted before or after the command name."""
    default = argparse.SUPPRESS if suppress else False
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true", default=default, help="debug messages")
    group.add_argument("-q", "--quiet", action="store_true", default=default, help="warnings and errors only")


class _Subcommands:
    """Wrap ``add_subparsers`` so that every sub-command also accepts -v/-q."""

    def __init__(self, action):
        self._action = action

    def __setattr__(self, name, value):
        if name == "_action":
            object.__setattr__(self, name, value)
        else:
            setattr(self._action, name, value)

    def add_parser(self, *args, **kwargs):
        sub = self._action.add_parser(*args, **kwargs)
        _add_verbosity(sub, suppress=True)
        return sub


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="eyeoncnv",
        description="EyeOnCNV - small copy-number events from targeted sequencing.",
        epilog=EPILOG,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    _add_verbosity(parser, suppress=False)
    commands = _Subcommands(parser.add_subparsers(dest="command", metavar="COMMAND"))
    commands.required = True

    p = commands.add_parser(
        "run",
        help="steps 1 to 4 in one go",
        description="Run steps 1 to 4 on the TSV files of step 0.",
    )
    p.add_argument("inputs", nargs="+", type=Path, help="TSV files or folders of TSV files")
    p.add_argument("-o", "--out-dir", type=Path, required=True, help="output folder")
    p.add_argument("--bed", type=Path, required=True, help="capture design BED")
    _add_intervals_options(p)
    _add_hgvs_options(p, required=False)
    _add_highlight_options(p, standalone=False)
    p.set_defaults(func=_cmd_run)

    p = commands.add_parser(
        "intervals",
        help="step 1 - events and recurrence across samples",
        description="Build events from depth-jump TSV files and count their recurrence.",
    )
    p.add_argument("inputs", nargs="+", type=Path, help="TSV files or folders of TSV files")
    p.add_argument("-o", "--out-dir", type=Path, required=True, help="output folder")
    _add_intervals_options(p)
    p.set_defaults(func=_cmd_intervals)

    p = commands.add_parser(
        "genes",
        help="step 2 - gene names from the capture BED",
        description="Add the names of the overlapping BED regions (column 'gene').",
    )
    p.add_argument("input", type=Path, help="Excel table from step 1")
    p.add_argument("--bed", type=Path, required=True, help="capture design BED")
    p.add_argument("-o", "--output", type=Path, help="output file [<input>_genes.xlsx]")
    p.set_defaults(func=_cmd_genes)

    p = commands.add_parser(
        "hgvs",
        help="step 3 - c. coordinates of the event boundaries (network)",
        description="Project start and end positions on the transcript (hgvs + UTA).",
    )
    p.add_argument("input", type=Path, help="Excel table from step 2")
    _add_hgvs_options(p, required=True)
    p.add_argument("-o", "--output", type=Path, help="output file [<input>_hgvs.xlsx]")
    p.set_defaults(func=_cmd_hgvs)

    p = commands.add_parser(
        "highlight",
        help="step 4 - flag and colour the events of interest",
        description="Add 'length' and 'of_interest' and fill the flagged rows.",
    )
    p.add_argument("input", type=Path, help="Excel table from step 3")
    p.add_argument("-o", "--output", type=Path, help="output file [<input>_highlight.xlsx]")
    _add_highlight_options(p, standalone=True)
    p.set_defaults(func=_cmd_highlight)

    return parser


def _cmd_run(args) -> None:
    from eyeoncnv.pipeline import run_pipeline

    run_pipeline(
        args.inputs,
        args.out_dir,
        args.bed,
        args.nm_map,
        gap=args.gap,
        min_overlap=args.min_overlap,
        max_occurrences=args.max_occurrences,
        history_dirs=args.history_dir,
        use_history=not args.no_history,
        assembly=args.assembly,
        db_url=args.uta_url,
        min_length=args.min_length,
        mode=args.mode,
        fill_color=args.color,
    )


def _cmd_intervals(args) -> None:
    from eyeoncnv.intervals import run_intervals

    run_intervals(
        args.inputs,
        args.out_dir,
        gap=args.gap,
        min_overlap=args.min_overlap,
        max_occurrences=args.max_occurrences,
        history_dirs=args.history_dir,
        use_history=not args.no_history,
    )


def _cmd_genes(args) -> None:
    from eyeoncnv.genes import run_genes

    run_genes(args.input, args.bed, args.output)


def _cmd_hgvs(args) -> None:
    from eyeoncnv.hgvs_coords import run_hgvs

    run_hgvs(args.input, args.nm_map, args.output, assembly=args.assembly, db_url=args.uta_url)


def _cmd_highlight(args) -> None:
    from eyeoncnv.highlight import run_highlight

    run_highlight(
        args.input,
        args.output,
        min_length=args.min_length,
        mode=args.mode,
        coding_only=not args.no_coding_filter,
        fill_color=args.color,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    setup_logging(-1 if args.quiet else 1 if args.verbose else 0)
    try:
        args.func(args)
    except KeyboardInterrupt:
        log.error("Interrupted")
        return 130
    except Exception as exc:
        if args.verbose:
            log.exception("%s", exc)
        else:
            log.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
