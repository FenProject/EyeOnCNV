#!/usr/bin/env bash
#
# EyeOnCNV - step 0: BAM -> depth-jump TSV
#
# For each BAM, list the positions of the capture design that flank an abrupt
# change of read depth. These positions are the input of the Python pipeline
# (`eyeoncnv intervals` / `eyeoncnv run`).
#
# Output: one headerless, tab-separated file per BAM, <OUT_DIR>/<sample>.tsv
#         columns: chromosome, position (1-based), depth
#
# Run `bam_to_depth_jumps.sh --help` for usage. Method: docs/small-events-pipeline.md

set -euo pipefail

VERSION="0.1.0"

usage() {
  cat <<'EOF'
Usage:
  bam_to_depth_jumps.sh --bed PANEL.bed --out-dir TSV_DIR --in-dir BAM_DIR
  bam_to_depth_jumps.sh --bed PANEL.bed --out-dir TSV_DIR SAMPLE1.bam [SAMPLE2.bam ...]

Required:
  -b, --bed FILE        Capture BED (.bed or .bed.gz). Depth is only computed inside it.
  -o, --out-dir DIR     Output directory for the TSV files (created if missing).
  -i, --in-dir DIR      Process every *.bam of DIR. BAM files can also be given
                        as positional arguments (both can be combined).

Filter parameters (defaults are the values of the validated lab script):
      --min-depth N     Ignore positions with depth < N.                        [30]
      --th-low X        Report a jump when |d2 - d1| / d1 >= X ...              [0.15]
      --th-high X       ... or >= X when d1 or d2 is below --depth-cut.        [0.30]
      --depth-cut N     Depth under which --th-high replaces --th-low.         [150]
      --win N           Drop a reported position that has no other reported
                        position within N bp on the same chromosome.            [50]
      --mode ge|le      ge: keep changes >= threshold; le: keep changes <= it.  [ge]

Other options:
      --samtools PATH   samtools executable.                                    [samtools]
  -h, --help            Show this help and exit.
      --version         Show the version and exit.

Every option can also be set with an environment variable of the same name in
upper case (BED, IN_DIR, OUT_DIR, MIN_DEPTH, TH_LOW, TH_HIGH, DEPTH_CUT, WIN,
MODE, SAMTOOLS); command-line options win. Example:
  export BED=panel_Covered.bed IN_DIR=bams OUT_DIR=tsv
  bash scripts/bam_to_depth_jumps.sh
EOF
}

log()  { printf '[INFO] %s\n' "$*" >&2; }
warn() { printf '[WARNING] %s\n' "$*" >&2; }
die()  { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Defaults (environment variables are honoured for backward compatibility)
# ---------------------------------------------------------------------------
BED="${BED:-}"
IN_DIR="${IN_DIR:-}"
OUT_DIR="${OUT_DIR:-}"
MIN_DEPTH="${MIN_DEPTH:-30}"
TH_LOW="${TH_LOW:-0.15}"
TH_HIGH="${TH_HIGH:-0.30}"
DEPTH_CUT="${DEPTH_CUT:-150}"
WIN="${WIN:-50}"
MODE="${MODE:-ge}"
SAMTOOLS="${SAMTOOLS:-samtools}"
AWK="${AWK:-awk}"   # only meant for testing other awk implementations

BAMS=()

need_value() {
  # $1 = option name, $2 = number of remaining arguments
  [[ "$2" -ge 2 ]] || die "Option $1 needs a value (see --help)"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -b|--bed)       need_value "$1" "$#"; BED="$2"; shift 2 ;;
    -i|--in-dir)    need_value "$1" "$#"; IN_DIR="$2"; shift 2 ;;
    -o|--out-dir)   need_value "$1" "$#"; OUT_DIR="$2"; shift 2 ;;
    --min-depth)    need_value "$1" "$#"; MIN_DEPTH="$2"; shift 2 ;;
    --th-low)       need_value "$1" "$#"; TH_LOW="$2"; shift 2 ;;
    --th-high)      need_value "$1" "$#"; TH_HIGH="$2"; shift 2 ;;
    --depth-cut)    need_value "$1" "$#"; DEPTH_CUT="$2"; shift 2 ;;
    --win)          need_value "$1" "$#"; WIN="$2"; shift 2 ;;
    --mode)         need_value "$1" "$#"; MODE="$2"; shift 2 ;;
    --samtools)     need_value "$1" "$#"; SAMTOOLS="$2"; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
    --version)      echo "bam_to_depth_jumps.sh ${VERSION}"; exit 0 ;;
    --)             shift; while [[ $# -gt 0 ]]; do BAMS+=("$1"); shift; done ;;
    -*)             die "Unknown option: $1 (see --help)" ;;
    *)              BAMS+=("$1"); shift ;;
  esac
done

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
is_number() { [[ "$1" =~ ^[0-9]+([.][0-9]+)?$ ]]; }

[[ -n "$BED" ]]     || die "Missing --bed (capture BED file)"
[[ -f "$BED" ]]     || die "BED file not found: $BED"
[[ -n "$OUT_DIR" ]] || die "Missing --out-dir (output directory)"
for spec in "min-depth=$MIN_DEPTH" "th-low=$TH_LOW" "th-high=$TH_HIGH" \
            "depth-cut=$DEPTH_CUT" "win=$WIN"; do
  is_number "${spec#*=}" || die "--${spec%%=*} must be a non-negative number (got '${spec#*=}')"
done
[[ "$MODE" == "ge" || "$MODE" == "le" ]] || die "--mode must be 'ge' or 'le' (got '$MODE')"
command -v "$SAMTOOLS" >/dev/null 2>&1 \
  || die "samtools not found ('$SAMTOOLS'). Install it or pass --samtools PATH"
command -v "$AWK" >/dev/null 2>&1 || die "awk not found ('$AWK')"

if [[ -n "$IN_DIR" ]]; then
  [[ -d "$IN_DIR" ]] || die "Input directory not found: $IN_DIR"
  shopt -s nullglob
  for bam in "$IN_DIR"/*.bam; do
    BAMS+=("$bam")
  done
  shopt -u nullglob
fi

if [[ ${#BAMS[@]} -eq 0 ]]; then
  warn "No .bam file to process (use --in-dir DIR and/or give BAM files)"
  exit 0
fi

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Filter 1 - depth jumps.
# Positions with depth >= min_depth are compared with the previous such
# position of the same chromosome. When the relative change |d2 - d1| / d1
# reaches the threshold (th_high if either depth is below depth_cut, th_low
# otherwise), both positions are printed (each line at most once).
# ---------------------------------------------------------------------------
read -r -d '' AWK_DEPTH_JUMPS <<'AWK' || true
BEGIN {
  FS = OFS = "\t"
  min_depth += 0; th_low += 0; th_high += 0; depth_cut += 0
}
function abs(x) { return x < 0 ? -x : x }
function relchange(a, b) {
  if (a == 0) return (b == 0 ? 0 : 1e99)
  return abs(b - a) / a
}
$3 >= min_depth {
  chr = $1; d = $3 + 0
  if (have_prev && chr == prev_chr) {
    rc = relchange(prev_d, d)
    th = ((prev_d < depth_cut || d < depth_cut) ? th_high : th_low)
    keep = (mode == "ge" ? rc >= th : rc <= th)
    if (keep) {
      if (prev_line != last) { print prev_line; last = prev_line }
      if ($0 != last)        { print $0;        last = $0 }
    }
  } else {
    last = ""
  }
  prev_chr = chr; prev_d = d; prev_line = $0; have_prev = 1
}
AWK

# ---------------------------------------------------------------------------
# Filter 2 - isolated positions.
# A position from filter 1 is kept only if the previous or the next one lies
# within `win` bp on the same chromosome.
# ---------------------------------------------------------------------------
read -r -d '' AWK_NEIGHBOURS <<'AWK' || true
BEGIN { FS = OFS = "\t"; win += 0 }
{
  chr = $1; pos = $2 + 0
  if (!have_prev) {
    prev_chr = chr; prev_pos = pos; prev_line = $0; prev_has = 0; have_prev = 1
    next
  }
  if (chr == prev_chr) {
    if (pos - prev_pos <= win) { prev_has = 1; cur_has = 1 } else { cur_has = 0 }
    if (prev_has) print prev_line
    prev_chr = chr; prev_pos = pos; prev_line = $0; prev_has = cur_has
  } else {
    if (prev_has) print prev_line
    prev_chr = chr; prev_pos = pos; prev_line = $0; prev_has = 0
  }
}
END { if (have_prev && prev_has) print prev_line }
AWK

CURRENT_TMP=""
cleanup() {
  if [[ -n "$CURRENT_TMP" && -f "$CURRENT_TMP" ]]; then
    rm -f "$CURRENT_TMP"
  fi
}
trap cleanup EXIT

process_bam() {
  local bam="$1" base out
  [[ -f "$bam" ]] || die "BAM file not found: $bam"
  base="$(basename "$bam" .bam)"
  out="$OUT_DIR/${base}.tsv"

  if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
    log "No index found, running samtools index: $bam"
    "$SAMTOOLS" index "$bam"
  fi

  log "$base -> $out"
  CURRENT_TMP="${out}.part"
  # -q 0 -Q 0: no base/mapping quality filter; -g SECONDARY: also count
  # secondary alignments; -d 0: no depth cap on old samtools (ignored >= 1.13).
  "$SAMTOOLS" depth -b "$BED" -q 0 -Q 0 -d 0 -g SECONDARY "$bam" \
    | "$AWK" -v min_depth="$MIN_DEPTH" -v th_low="$TH_LOW" -v th_high="$TH_HIGH" \
             -v depth_cut="$DEPTH_CUT" -v mode="$MODE" "$AWK_DEPTH_JUMPS" \
    | "$AWK" -v win="$WIN" "$AWK_NEIGHBOURS" \
    > "$CURRENT_TMP"
  mv -f "$CURRENT_TMP" "$out"
  CURRENT_TMP=""
}

log "BED: $BED | output: $OUT_DIR | ${#BAMS[@]} BAM file(s)"
log "min-depth=$MIN_DEPTH th-low=$TH_LOW th-high=$TH_HIGH depth-cut=$DEPTH_CUT win=$WIN mode=$MODE"
for bam in "${BAMS[@]}"; do
  process_bam "$bam"
done
log "Done: TSV files written to $OUT_DIR"
