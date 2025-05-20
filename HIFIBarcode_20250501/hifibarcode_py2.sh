#!/usr/bin/env bash
set -euo pipefail

# ─────────── Edit these two lines ───────────
PRIMER_ASSIGN="/home/liulab/lizicong/software/metabarcoding_Li/96_pairs_primer.txt"
PRIMER_POLISH="/home/liulab/lizicong/software/metabarcoding_Li/primer"
# ────────────────────────────────────────────

usage() {
  echo "Usage: $0 -fq1 FASTQ_R1 -fq2 FASTQ_R2 -n THREADS -i INDEX -kdb KRAKEN2_DB -o OUTPREFIX"
  exit 1
}

# parse args
if [ $# -lt 12 ]; then
  usage
fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -fq1)        FQ1="$2"; shift 2;;
    -fq2)        FQ2="$2"; shift 2;;
    -num_threads|-n) THREADS="$2"; shift 2;;
    -i)          INDEX="$2"; shift 2;;
    -kdb)        KRAKEN2_DB="$2"; shift 2;;
    -o)          OUTPRE="$2"; shift 2;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

# sanity checks
: "${FQ1:?Missing -fq1}"
: "${FQ2:?Missing -fq2}"
: "${THREADS:?Missing -n}"
: "${INDEX:?Missing -i}"
: "${KRAKEN2_DB:?Missing -kdb}"
: "${OUTPRE:?Missing -o}"

# create output directories
ASSIGN_DIR="${OUTPRE}_assign"
BUILDEND_DIR="${OUTPRE}_buildend"
CHAIN_DIR="${OUTPRE}_chain"
GAPFILL_DIR="${OUTPRE}_gapfill"
POLISH_DIR="${OUTPRE}_polish"
mkdir -p "$ASSIGN_DIR" "$BUILDEND_DIR" "$CHAIN_DIR" "$GAPFILL_DIR" "$POLISH_DIR"

echo ">> Step 1: Assign reads to primers"
Assign20250501.py -fq1 "$FQ1" -fq2 "$FQ2" -num_threads "$THREADS" -index "$INDEX" \
  -primer "$PRIMER_ASSIGN" -outpre "$OUTPRE" |& tee "${ASSIGN_DIR}/${OUTPRE}.assign.log"

echo ">> Step 2: Build ends"
Buildend20250501.py -outpre "$OUTPRE" -index "$INDEX" -threads "$THREADS" \
  -list "${ASSIGN_DIR}/assign.list" |& tee "${BUILDEND_DIR}/${OUTPRE}.buildend.log"

echo ">> Step 2.5: Kraken2 classification & filtering"
# backup original ends
cp "${BUILDEND_DIR}/buildends.fasta" "${BUILDEND_DIR}/raw_buildends.fasta"
# run kraken2
if ! command -v kraken2 >/dev/null 2>&1; then
  echo "Error: kraken2 not found in PATH" >&2
  exit 1
fi
kraken2 --db "$KRAKEN2_DB" --output "${BUILDEND_DIR}/${OUTPRE}_kraken2_out.txt" \
        --report "${BUILDEND_DIR}/${OUTPRE}_kraken2_report.txt" "${BUILDEND_DIR}/buildends.fasta"
# extract unclassified IDs
echo ">> Extracting unclassified IDs"
awk '$1=="U" { print $2 }' "${BUILDEND_DIR}/${OUTPRE}_kraken2_out.txt" \
  > "${BUILDEND_DIR}/keep_ids.txt"
# filter FASTA using two-file awk
if [ ! -s "${BUILDEND_DIR}/keep_ids.txt" ]; then
  echo "Warning: no unclassified IDs found; buildends.fasta will be emptied." >&2
  > "${BUILDEND_DIR}/buildends.fasta"
else
  echo ">> Filtering buildends.fasta"
  awk 'NR==FNR { keep[$1]; next } \
       /^>/ { split(substr($0,2),a,/\s+/); hdr=a[1]; flag=(hdr in keep) } \
       flag { print }' \
    "${BUILDEND_DIR}/keep_ids.txt" "${BUILDEND_DIR}/raw_buildends.fasta" \
    > "${BUILDEND_DIR}/buildends.fasta"
  echo "Kept $(wc -l < "${BUILDEND_DIR}/keep_ids.txt") sequences."
fi

echo ">> Step 3: Chain sequences"
Chain20250501.py -outpre "$OUTPRE" -num_threads "$THREADS" \
  -mf "${ASSIGN_DIR}/assign_middle.ssam" |& tee "${CHAIN_DIR}/${OUTPRE}.chain.log"

echo ">> Step 4: Gap fill"
Gapfill20250501.py -outpre "$OUTPRE" -ends "${BUILDEND_DIR}/buildends.fasta" \
  -middle "${CHAIN_DIR}/chain.fasta" -t "$THREADS" |& tee "${GAPFILL_DIR}/${OUTPRE}.gapfill.log"

echo ">> Step 5: Polish sequences"
Polish20250501.py -d "$GAPFILL_DIR" -outpre "$OUTPRE" -primer "$PRIMER_POLISH" \
  |& tee "${POLISH_DIR}/${OUTPRE}.polish.log"

echo "Pipeline completed successfully!"
