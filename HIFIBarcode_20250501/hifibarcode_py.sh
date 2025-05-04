#!/usr/bin/env bash
set -euo pipefail

# ─────────── Edit these two lines ───────────
PRIMER_ASSIGN="/home/liulab/lizicong/software/metabarcoding_Li/96_pairs_primer.txt"   # primer file for Assign20250501.py
PRIMER_POLISH="/home/liulab/lizicong/software/metabarcoding_Li/primer"                 # primer file for Polish20250501.py
# ────────────────────────────────────────────

usage() {
  echo "Usage: $0 -fq1 FASTQ_R1 -fq2 FASTQ_R2 -n THREADS -i INDEX -o OUTPREFIX"
  echo
  echo "  -fq1    forward reads file (e.g. test_1.fq.gz)"
  echo "  -fq2    reverse reads file (e.g. test_2.fq.gz)"
  echo "  -n       number of threads (used for -num_threads, -threads, and -t)"
  echo "  -i       index (tag length)"
  echo "  -o       output prefix (will be used to name all folders/files)"
  exit 1
}

# parse args
if [ $# -lt 10 ]; then
  usage
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -fq1)  FQ1="$2"; shift 2;;
    -fq2)  FQ2="$2"; shift 2;;
    -num_threads|-threads|-n|-t)
          THREADS="$2"; shift 2;;
    -i|-index)
          INDEX="$2"; shift 2;;
    -o|-outpre)
          OUTPRE="$2"; shift 2;;
    *) 
          echo "Unknown argument: $1"
          usage
          ;;
  esac
done

# sanity checks
: "${FQ1:?Missing -fq1}"
: "${FQ2:?Missing -fq2}"
: "${THREADS:?Missing thread count}"
: "${INDEX:?Missing index}"
: "${OUTPRE:?Missing output prefix}"

# make output directories
ASSIGN_DIR="${OUTPRE}_assign"
BUILDEND_DIR="${OUTPRE}_buildend"
CHAIN_DIR="${OUTPRE}_chain"
GAPFILL_DIR="${OUTPRE}_gapfill"
POLISH_DIR="${OUTPRE}_polish"

mkdir -p "$ASSIGN_DIR" "$BUILDEND_DIR" "$CHAIN_DIR" "$GAPFILL_DIR" "$POLISH_DIR"

echo "Starting pipeline with:"
echo "  FQ1        = $FQ1"
echo "  FQ2        = $FQ2"
echo "  THREADS    = $THREADS"
echo "  INDEX      = $INDEX"
echo "  OUTPRE     = $OUTPRE"
echo "  PRIMER1    = $PRIMER_ASSIGN"
echo "  PRIMER2    = $PRIMER_POLISH"
echo

# 1. Assign
echo ">> Running Assign20250501.py"
Assign20250501.py \
  -fq1 "$FQ1" \
  -fq2 "$FQ2" \
  -num_threads "$THREADS" \
  -index "$INDEX" \
  -primer "$PRIMER_ASSIGN" \
  -outpre "$OUTPRE" \
  |& tee "${ASSIGN_DIR}/${OUTPRE}.assign.log"

# 2. Buildend
echo ">> Running Buildend20250501.py"
Buildend20250501.py \
  -outpre "$OUTPRE" \
  -index "$INDEX" \
  -threads "$THREADS" \
  -list "${ASSIGN_DIR}/assign.list" \
  |& tee "${BUILDEND_DIR}/${OUTPRE}.buildend.log"

# 3. Chain
echo ">> Running Chain20250501.py"
Chain20250501.py \
  -outpre "$OUTPRE" \
  -num_threads "$THREADS" \
  -mf "${ASSIGN_DIR}/assign_middle.ssam" \
  |& tee "${CHAIN_DIR}/${OUTPRE}.chain.log"

# 4. Gapfill
echo ">> Running Gapfill20250501.py"
Gapfill20250501.py \
  -outpre "$OUTPRE" \
  -ends "${BUILDEND_DIR}/buildends.fasta" \
  -middle "${CHAIN_DIR}/chain.fasta" \
  -t "$THREADS" \
  |& tee "${GAPFILL_DIR}/${OUTPRE}.gapfill.log"

# 5. Polish
echo ">> Running Polish20250501.py"
Polish20250501.py \
  -d "${GAPFILL_DIR}" \
  -outpre "$OUTPRE" \
  -primer "$PRIMER_POLISH" \
  |& tee "${POLISH_DIR}/${OUTPRE}.polish.log"

echo "Pipeline completed successfully!"
