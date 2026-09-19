"""Step 3 - c. coordinates (HGVS) of the event boundaries.

For each event, the transcript is taken from the ``gene`` column written by step 2
(``GENE-NMxxxxxx``), completed with its version from a transcript table
(columns ``genes`` and ``NM``, e.g. ``BRCA2 | NM_000059.4``). The genomic start and
end positions are then projected on that transcript with the `hgvs` package and
the UTA database (https://github.com/biocommons/hgvs), which needs network access.

Added columns: ``gene_symbol``, ``nm_input``, ``nm_version``, ``nm_source``,
``c_start``, ``hgvs_start``, ``c_end``, ``hgvs_end`` and ``hgvs_error``.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from pathlib import Path
from typing import Protocol

import pandas as pd

from eyeoncnv.common import (
    chrom_key,
    is_missing,
    log,
    pick_column,
    read_table,
    read_workbook,
    require_column,
    tagged_path,
    write_excel,
)

DEFAULT_ASSEMBLY = "GRCh37"

# RefSeq accessions of the primary assembly chromosomes (from bioutils.assemblies).
NC_ACCESSIONS: dict[str, dict[str, str]] = {
    "GRCh37": {
        "1": "NC_000001.10", "2": "NC_000002.11", "3": "NC_000003.11", "4": "NC_000004.11",
        "5": "NC_000005.9", "6": "NC_000006.11", "7": "NC_000007.13", "8": "NC_000008.10",
        "9": "NC_000009.11", "10": "NC_000010.10", "11": "NC_000011.9", "12": "NC_000012.11",
        "13": "NC_000013.10", "14": "NC_000014.8", "15": "NC_000015.9", "16": "NC_000016.9",
        "17": "NC_000017.10", "18": "NC_000018.9", "19": "NC_000019.9", "20": "NC_000020.10",
        "21": "NC_000021.8", "22": "NC_000022.10", "X": "NC_000023.10", "Y": "NC_000024.9",
        "MT": "NC_012920.1",
    },
    "GRCh38": {
        "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
        "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
        "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
        "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
        "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
        "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10",
        "MT": "NC_012920.1",
    },
}  # fmt: skip
ASSEMBLY_ALIASES = {"grch37": "GRCh37", "hg19": "GRCh37", "grch38": "GRCh38", "hg38": "GRCh38"}

GENE_TRANSCRIPT_RE = re.compile(
    r"^\s*(?P<gene>[A-Za-z0-9._\-]+)\s*-\s*(?P<nm>NM_?\d+)(?:\.\d+)?\s*$", re.IGNORECASE
)
GENE_SYMBOL_RE = re.compile(r"^[A-Za-z0-9._\-]+$")
NM_RE = re.compile(r"^NM_?(\d+)(\.\d+)?$", re.IGNORECASE)

SOURCE_BY_NM = "by_nm_base"
SOURCE_BY_GENE = "by_gene"
SOURCE_MISSING = "missing"


def normalize_assembly(name: str) -> str:
    key = str(name).strip()
    resolved = ASSEMBLY_ALIASES.get(key.lower(), key)
    if resolved not in NC_ACCESSIONS:
        raise ValueError(f"Unsupported assembly '{name}' (use GRCh37/hg19 or GRCh38/hg38)")
    return resolved


def chrom_to_nc(chrom, assembly: str = DEFAULT_ASSEMBLY) -> str:
    """``13``/``chr13`` -> ``NC_000013.10`` (GRCh37)."""
    table = NC_ACCESSIONS[normalize_assembly(assembly)]
    key = chrom_key(chrom)
    if key not in table:
        raise ValueError(f"Unknown chromosome: {chrom}")
    return table[key]


def normalize_nm(value) -> tuple[str, str | None] | None:
    """``nm000059.4`` -> (``NM_000059``, ``.4``); None if not an NM accession."""
    match = NM_RE.match(str(value).strip())
    if not match:
        return None
    return f"NM_{match.group(1)}", match.group(2)


def parse_gene_transcript(field) -> tuple[str | None, str | None]:
    """Parse the ``gene`` column of step 2.

    ``BRCA2-NM000059``, ``BRCA2-NM_000059`` and ``BRCA2-NM_000059.4`` give
    (``BRCA2``, ``NM_000059``). A plain symbol (``BRCA2``) gives (``BRCA2``, None).
    With several regions (``A-NM_1;B-NM_2``) only the first one is used.
    """
    if is_missing(field):
        return None, None
    text = str(field).split(";")[0].strip()
    if not text:
        return None, None
    match = GENE_TRANSCRIPT_RE.match(text)
    if match:
        base, _version = normalize_nm(match.group("nm"))
        return match.group("gene"), base
    if GENE_SYMBOL_RE.match(text):
        return text, None
    return None, None


def load_transcript_map(path) -> tuple[dict[str, str], dict[str, str]]:
    """Read the transcript table: returns (NM base -> versioned NM, GENE -> versioned NM).

    Accepted column names: ``genes``/``gene``/``symbol`` and ``NM``/``transcript``/
    ``refseq``. When a gene or an NM appears several times, the first row wins.
    """
    df = read_table(path)
    gene_col = pick_column(df.columns, ["genes", "gene", "symbol", "gene_symbol"])
    nm_col = pick_column(df.columns, ["nm", "transcript", "refseq"])
    if gene_col is None or nm_col is None:
        raise ValueError(
            f"{Path(path).name}: columns 'genes' and 'NM' not found "
            f"(columns present: {', '.join(map(str, df.columns))})"
        )
    by_nm: dict[str, str] = {}
    by_gene: dict[str, str] = {}
    for gene, nm in zip(df[gene_col], df[nm_col]):
        if is_missing(nm) or not str(nm).strip():
            continue
        parsed = normalize_nm(nm)
        if parsed is None:
            log.warning("Transcript table: '%s' is not an NM accession, ignored", nm)
            continue
        base, version = parsed
        versioned = base + (version or "")
        by_nm.setdefault(base, versioned)
        if not is_missing(gene) and str(gene).strip():
            by_gene.setdefault(str(gene).strip().upper(), versioned)
    return by_nm, by_gene


def resolve_transcript(
    gene: str | None, nm_base: str | None, by_nm: Mapping[str, str], by_gene: Mapping[str, str]
) -> tuple[str | None, str]:
    """Versioned transcript for an event: by NM first, then by gene symbol."""
    if nm_base and nm_base in by_nm:
        return by_nm[nm_base], SOURCE_BY_NM
    if gene and gene.upper() in by_gene:
        return by_gene[gene.upper()], SOURCE_BY_GENE
    return None, SOURCE_MISSING


class Projector(Protocol):
    def g_to_c(self, transcript: str, chrom, position: int) -> tuple[str, str]: ...


class HgvsProjector:
    """Project genomic positions on transcripts with hgvs + UTA (network needed)."""

    def __init__(
        self,
        assembly: str = DEFAULT_ASSEMBLY,
        alt_aln_method: str = "splign",
        db_url: str | None = None,
    ):
        try:
            from hgvs.assemblymapper import AssemblyMapper
            from hgvs.dataproviders.uta import connect
            from hgvs.parser import Parser
        except ImportError as exc:
            raise ImportError(
                "Step 3 needs the 'hgvs' package: pip install 'eyeoncnv[hgvs]' "
                "(see docs/installation.md; not supported on native Windows)"
            ) from exc
        self.assembly = normalize_assembly(assembly)
        log.info("Connecting to UTA%s...", f" ({db_url})" if db_url else "")
        try:
            self._hdp = connect(db_url=db_url)
        except Exception as exc:
            raise ConnectionError(
                f"Cannot connect to the UTA database: {exc}. The public server needs "
                "outgoing access to uta.biocommons.org:5432; otherwise run a local UTA "
                "and pass --uta-url (see docs/installation.md)."
            ) from exc
        self._mapper = AssemblyMapper(
            self._hdp,
            assembly_name=self.assembly,
            alt_aln_method=alt_aln_method,
            replace_reference=True,
        )
        self._parser = Parser()
        self._cache: dict[tuple[str, str | None, int], tuple[str, str]] = {}

    def g_to_c(self, transcript: str, chrom, position: int) -> tuple[str, str]:
        """(``467+5``, ``NM_000059.4:c.467+5=``) for a genomic position."""
        key = (transcript, chrom_key(chrom), int(position))
        if key not in self._cache:
            nc = chrom_to_nc(chrom, self.assembly)
            var_g = self._parser.parse_hgvs_variant(f"{nc}:g.{int(position)}=")
            var_c = self._mapper.g_to_c(var_g, transcript)
            self._cache[key] = (str(var_c.posedit.pos), str(var_c))
        return self._cache[key]


def _as_position(value) -> int | None:
    if is_missing(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not number.is_integer():
        return None
    return int(number)


def add_c_coordinates(
    df: pd.DataFrame,
    projector: Projector,
    by_nm: Mapping[str, str],
    by_gene: Mapping[str, str],
) -> pd.DataFrame:
    """Return a copy of ``df`` with transcript and c. coordinate columns."""
    gene_col = require_column(df, ["gene"], "gene (run step 2 first)")
    chr_col = require_column(df, ["chr", "chrom", "chromosome"], "chr")
    start_col = require_column(df, ["start"], "start")
    end_col = require_column(df, ["end"], "end")

    out = df.copy()
    parsed = [parse_gene_transcript(value) for value in out[gene_col]]
    resolved = [resolve_transcript(g, nm, by_nm, by_gene) for g, nm in parsed]
    out["gene_symbol"] = [gene for gene, _ in parsed]
    out["nm_input"] = [nm for _, nm in parsed]
    out["nm_version"] = [transcript for transcript, _ in resolved]
    out["nm_source"] = [source for _, source in resolved]

    transcripts = [transcript for transcript, _ in resolved]
    errors: list[list[str]] = [[] for _ in range(len(out))]
    for side, column in (("start", start_col), ("end", end_col)):
        coords, full = [], []
        for row, (transcript, chrom, raw_pos) in enumerate(
            zip(transcripts, out[chr_col].tolist(), out[column].tolist())
        ):
            coord = hgvs = None
            position = _as_position(raw_pos)
            if position is None:
                errors[row].append(f"{side}: position is not an integer ({raw_pos})")
            elif is_missing(transcript):
                errors[row].append(f"{side}: no versioned NM found (neither by NM nor by gene)")
            else:
                try:
                    coord, hgvs = projector.g_to_c(transcript, chrom, position)
                except Exception as exc:  # hgvs raises many different error types
                    errors[row].append(f"{side}: {exc}")
            coords.append(coord)
            full.append(hgvs)
        out[f"c_{side}"] = coords
        out[f"hgvs_{side}"] = full
    out["hgvs_error"] = [_merge_errors(messages) for messages in errors]
    return out


def _merge_errors(messages: list[str]) -> str:
    if not messages:
        return ""
    bodies = [m.split(": ", 1)[1] for m in messages]
    if len(messages) == 2 and bodies[0] == bodies[1]:
        return bodies[0]
    return "; ".join(messages)


def run_hgvs(
    input_path,
    transcript_map,
    output_path=None,
    assembly: str = DEFAULT_ASSEMBLY,
    db_url: str | None = None,
    projector: Projector | None = None,
) -> Path:
    """Step 3 on every sheet of ``input_path``; writes ``<input>_hgvs.xlsx`` by default."""
    input_path = Path(input_path)
    by_nm, by_gene = load_transcript_map(transcript_map)
    log.info("Transcript table: %d NM, %d gene(s)", len(by_nm), len(by_gene))
    sheets = read_workbook(input_path)
    if projector is None:
        projector = HgvsProjector(assembly=assembly, db_url=db_url)
    projected = {}
    for name, df in sheets.items():
        projected[name] = add_c_coordinates(df, projector, by_nm, by_gene)
        n_ok = int(projected[name]["c_start"].notna().sum())
        log.info("Sheet %s: %d/%d start position(s) projected", name, n_ok, len(df))
    output = Path(output_path) if output_path else tagged_path(input_path, "hgvs")
    write_excel(projected, output)
    log.info("Written: %s", output)
    return output
