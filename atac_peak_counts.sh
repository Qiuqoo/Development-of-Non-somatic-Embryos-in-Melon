#!/usr/bin/env bash

BED="/data3/users/Qiuqoo/CM_ZJU/hintatac/genome/gene_regions.bed"

BAM_DIR="/data3/users/Qiuqoo/CM_ZJU/hintatac/bam"

SAMPLES=(NSE GE HE PTE LTE PCE LCE)

REPS=(1 2 3)

mkdir -p region_counts

export BED BAM_DIR

parallel -j 6 '
  SAMPLE={1}
  REP={2}
  BAM_FILE="$BAM_DIR/${SAMPLE}_${REP}.bam"
  OUT="region_counts/${SAMPLE}_${REP}_read_counts.txt"
  if [[ -f "$BAM_FILE" ]]; then
    echo "Processing $BAM_FILE..."
    bedtools coverage -a "$BED" -b "$BAM_FILE" -counts > "$OUT"
  else
    echo "File not found: $BAM_FILE" >&2
  fi
' ::: "${SAMPLES[@]}" ::: "${REPS[@]}"
