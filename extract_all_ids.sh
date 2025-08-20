#!/usr/bin/env bash

# Usage: ./extract_all_ids.sh [-f] /path/to/phmmer_outputs
#   -f   Force overwrite of existing _hits.txt files

set -euo pipefail

force=false
while getopts ":f" opt; do
  case $opt in
    f) force=true ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

indir=${1:-}
if [[ -z "$indir" ]]; then
  echo "Usage: $0 [-f] /path/to/phmmer_outputs"
  exit 1
fi

shopt -s nullglob

for file in "$indir"/*_phmmer_output.txt; do
  base=$(basename "$file" _phmmer_output.txt)
  outfile="$indir/${base}_hits.txt"

  if [[ -f "$outfile" && $force == false ]]; then
    echo "Skipping $file (output already exists: $outfile)"
    continue
  fi

  echo "Processing $file -> $outfile"

  awk '
    BEGIN {
      inblock=0
      header_seen=0
      seqcol=0
    }

    /Scores[[:space:]]+for[[:space:]]+complete[[:space:]]+sequences/ {
      inblock=1
      header_seen=0
      seqcol=0
      next
    }

    inblock && /inclusion[[:space:]]+threshold/ {
      inblock=0
      next
    }

    inblock && /^[[:space:]]*$/ {
      inblock=0
      next
    }

    inblock {
      sub(/\r$/, "", $0)

      if ($0 ~ /---[[:space:]]+full[[:space:]]+sequence[[:space:]]+---/) next
      if ($0 ~ /^[[:space:]]*-{2,}([[:space:]]+-{2,})+/) next

      if (header_seen==0 && $1=="E-value") {
        for (i=1; i<=NF; i++) if ($i=="Sequence") seqcol=i
        header_seen=1
        next
      }

      if (header_seen==1 && seqcol>0) {
        if (NF>=seqcol) {
          val = $(seqcol)
          if (val !~ /^-+$/ && val != "Sequence" && val != "") print val
        }
      }
    }
  ' "$file" | sort -u > "$outfile"
done

echo "All ID lists generated."
