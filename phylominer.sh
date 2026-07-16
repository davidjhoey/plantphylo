#!/usr/bin/env bash

set -euo pipefail

ORIGINAL_ARGS=("$@")

usage() {
  cat <<'EOF'
Usage:
  phylominer.sh [options] query.fasta /path/to/databases

Required:
  query.fasta              Protein or CDS FASTA query
  /path/to/databases       Directory containing FASTA databases

Options:
  -f                       Overwrite existing outputs
  --threads N              Number of threads for phmmer and hmmsearch (default: 1)
  --keep-temp              Keep translated files, hit lists, and search tables
  --motif-hmm HMM_ID       Filter extracted proteins by HMM profile
  --pfam-db PATH           HMM database for --motif-hmm
  -h, --help               Show this help message and exit

Notes:
  - If the query is nucleotide, it is translated with transeq using frame 1.
  - Database files are treated as protein or CDS automatically.
  - If a database is nucleotide, a protein FASTA and a CDS FASTA are both written.
  - HMM filtering is applied after phmmer hit extraction.
  - Input databases should be .fa format (these are protected from deletion).
  - DEPENDENCIES: HMMER (phmmer, hmmsearch, hmmfetch), EMBOSS (transeq), SeqKit (seqkit). 
  - If using HMM profile filtering, download Pfam-A.hmm from https://www.ebi.ac.uk/interpro/download/Pfam/ (last accessed 02-07-2026).
EOF
}

safe_rm() {
    local file="$1"

    [[ -n "$file" ]] || return

    case "$file" in
        *.fa)
            echo "[WARN] Refusing to delete .fa file:"
            echo "       $file"
            return
            ;;
    esac

    rm -f -- "$file"
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: required command not found: $cmd" >&2
    exit 1
  fi
}

safe_count() {
  local file="$1"
  if [[ -s "$file" ]]; then
    grep -c '^>' "$file" || echo 0
  else
    echo 0
  fi
}

standardize_fasta_headers() {
  local input="$1"
  local output="$2"
  sed -E 's/^>([^[:space:]]+).*/>\1/' "$input" > "$output"
}

is_nucleotide_fasta() {
  local file="$1"
  local seq
  seq=$(grep -v '^>' "$file" | tr -d '\n\r' | head -c 50000)
  [[ -n "$seq" ]] || return 1

  awk -v seq="$seq" '
    BEGIN {
      n = length(seq)
      if (n == 0) exit 1
      gsub(/[^ACGTUNacgtun]/, "", seq)
      if ((length(seq) / n) >= 0.90) exit 0
      exit 1
    }
  '
}

translate_cds_to_protein() {

  local cds_fa="$1"
  local prot_fa="$2"

  transeq \
    -sequence "$cds_fa" \
    -outseq "$prot_fa" \
    -frame 1

  sed -i 's/^\(>[^[:space:]]*\)_1\([[:space:]]\|$\)/\1\2/' "$prot_fa"

  standardize_fasta_headers "$prot_fa" "$prot_fa.tmp"

  mv "$prot_fa.tmp" "$prot_fa"
}

extract_ids_from_phmmer() {
    local input="$1"
    local output="$2"

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

        if (!header_seen && $1=="E-value") {
            for (i=1; i<=NF; i++)
                if ($i=="Sequence")
                    seqcol=i
            header_seen=1
            next
        }

        if (header_seen && seqcol>0 && NF>=seqcol) {
            id=$(seqcol)
            if (id != "" && id !~ /^-+$/)
                print id
        }
    }
    ' "$input" | sort -u > "$output"

    [[ -s "$output" ]] || touch "$output"
}

extract_ids_from_hmmsearch() {
  local tbl="$1"
  local out="$2"

  awk '!/^#/ {print $1}' "$tbl" | sort -u > "$out"
}

extract_cds_with_failsafe_ids() {

  local prot_out="$1"
  local cds_db="$2"
  local cds_out="$3"
  local id_list="$4"

  echo "[INFO] Standard CDS extraction failed. Trying alternative ID strategies."

  # Failsafe 1: use complete protein headers
  grep "^>" "$prot_out" | sed 's/^>//' > "${id_list}.failsafe1.txt"

  seqkit grep -n -r -f "${id_list}.failsafe1.txt" "$cds_db" > "$cds_out" || true

  local cds_final
  cds_final=$(safe_count "$cds_out")

  if [[ "$cds_final" -gt 0 ]]; then
    echo "[INFO] Failsafe ID extraction 1 recovered $cds_final CDS sequences"
    return 0
  fi


  # Failsafe 2: use ID after final pipe
  echo "[INFO] Failsafe 1 failed. Trying pipe-based ID extraction."

  grep "^>" "$prot_out" \
    | sed 's/^>.*|//; s/[[:space:]].*$//' \
    > "${id_list}.failsafe2.txt"

  seqkit grep -n -f "${id_list}.failsafe2.txt" "$cds_db" > "$cds_out" || true

  cds_final=$(safe_count "$cds_out")

  if [[ "$cds_final" -gt 0 ]]; then
    echo "[INFO] Failsafe ID extraction 2 recovered $cds_final CDS sequences"
    return 0
  fi


  echo "[WARN] CDS extraction failed after alternative ID strategies"
  return 1
}


resolve_pfam_acc() {
  local acc="$1"
  local pfam_db="$2"
  local resolved=""

  resolved=$(grep -m1 -E "^ACC[[:space:]]+${acc}(\.[0-9]+)?$" "$pfam_db" | awk '{print $2}' || true)
  if [[ -n "$resolved" ]]; then
    printf '%s\n' "$resolved"
    return 0
  fi

  resolved=$(grep -m1 -E "^ACC[[:space:]]+${acc}\.[0-9]+" "$pfam_db" | awk '{print $2}' || true)
  if [[ -n "$resolved" ]]; then
    printf '%s\n' "$resolved"
    return 0
  fi

  return 1
}

extract_pfam_hmm() {
  local acc="$1"
  local pfam_db="$2"
  local out_hmm="$3"

  if [[ ! -f "${pfam_db}.h3i" || ! -f "${pfam_db}.h3m" ]]; then
    echo "ERROR: Pfam database must be pressed with hmmpress: $pfam_db" >&2
    return 1
  fi

  if hmmfetch -o "$out_hmm" "$pfam_db" "$acc" >/dev/null 2>&1; then
    return 0
  fi

  return 1
}

filter_proteins_by_hmms() {
  local fasta_in="$1"
  local fasta_out="$2"
  local workdir="$3"
  shift 3
  local hmm_ids=("$@")

  cp "$fasta_in" "$fasta_out"

  local hmm_id hmm_file tblout hits
  local current_ids

  current_ids=$(mktemp "$workdir/current_ids.XXXXXX")

  grep '^>' "$fasta_in" | sed 's/^>//' > "$current_ids"

  for hmm_id in "${hmm_ids[@]}"; do

    hmm_file="$workdir/${hmm_id}.hmm"
    tblout="$workdir/${hmm_id}.tbl"
    hits="$workdir/${hmm_id}.ids"

    if ! extract_pfam_hmm "$hmm_id" "$PFAM_DB" "$hmm_file"; then
      echo "WARN: could not fetch HMM $hmm_id"
      continue
    fi

    if ! hmmsearch --cpu "$THREADS" --tblout "$tblout" "$hmm_file" "$fasta_out" >/dev/null 2>&1; then
      echo "WARN: hmmsearch failed for $hmm_id"
      continue
    fi

    extract_ids_from_hmmsearch "$tblout" "$hits"

    if [[ ! -s "$hits" ]]; then
      echo "WARN: no hits returned for $hmm_id"
      : > "$fasta_out"
      safe_rm "$current_ids"
      return 0
    fi

    grep -Fxf "$hits" "$current_ids" > "$current_ids.tmp" || true
    mv "$current_ids.tmp" "$current_ids"

    if [[ ! -s "$current_ids" ]]; then
      : > "$fasta_out"
      safe_rm "$current_ids"
      return 0
    fi

  done

  seqkit grep -n -f "$current_ids" "$fasta_in" > "$fasta_out" || true

  safe_rm "$current_ids"
}

cleanup_temp_files() {
  if [[ "$KEEP_TEMP" == true ]]; then
    return 0
  fi

  for f in \
    "$PHMMER_TXT" \
    "$PHMMER_TBL" \
    "$ID_LIST" \
    "$PROT_DB_CLEAN" \
    "$CDS_DB_CLEAN" \
    "$CDS_ID_LIST" \
    "$PROT_OUT.filtered" \
    "$WORKDIR"/*.tbl \
    "$WORKDIR"/*.ids \
    "${CDS_ID_LIST}.failsafe1.txt" \
    "${CDS_ID_LIST}.failsafe2.txt" \
    "$WORKDIR"/*.hmm 
   do
    safe_rm "$f"
   done 
}

FORCE=false
KEEP_TEMP=false
THREADS=1
PFAM_DB=""
MOTIF_HMMS=()
VERSION="1.2.1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f) FORCE=true; shift ;;
    --threads) THREADS="${2:?Missing value for --threads}"; shift 2 ;;
    --keep-temp) KEEP_TEMP=true; shift ;;
    --motif-hmm) MOTIF_HMMS+=("${2:?Missing value for --motif-hmm}"); shift 2 ;;
    --pfam-db) PFAM_DB="${2:?Missing value for --pfam-db}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "ERROR: unknown option: $1" >&2; usage; exit 1 ;;
    *) break ;;
  esac
done

if [[ $# -ne 2 ]]; then
  usage
  exit 1
fi

QUERY="$1"
DB_DIR="$2"

if [[ ! -f "$QUERY" ]]; then
  echo "ERROR: query file not found: $QUERY" >&2
  exit 1
fi

if [[ ! -d "$DB_DIR" ]]; then
  echo "ERROR: database directory not found: $DB_DIR" >&2
  exit 1
fi

require_cmd phmmer
require_cmd hmmsearch
require_cmd seqkit
require_cmd transeq

if [[ ${#MOTIF_HMMS[@]} -gt 0 && -z "$PFAM_DB" ]]; then
  echo "ERROR: --motif-hmm requires --pfam-db" >&2
  exit 1
fi

if [[ -n "$PFAM_DB" && ! -f "$PFAM_DB" ]]; then
  echo "ERROR: Pfam database not found: $PFAM_DB" >&2
  exit 1
fi

if [[ ${#MOTIF_HMMS[@]} -gt 0 ]]; then
  RESOLVED_MOTIFS=()
  for motif in "${MOTIF_HMMS[@]}"; do
    if resolved=$(resolve_pfam_acc "$motif" "$PFAM_DB" 2>/dev/null); then
      RESOLVED_MOTIFS+=("$resolved")
    else
      RESOLVED_MOTIFS+=("$motif")
    fi
  done
  MOTIF_HMMS=("${RESOLVED_MOTIFS[@]}")
fi



BASE_DIR="$DB_DIR"
WORKDIR="$BASE_DIR/phylominer_work"
PROT_DIR="$BASE_DIR/proteins"
CDS_DIR="$BASE_DIR/CDS"
LOG_FILE="$BASE_DIR/phylominer.txt"

{
  echo "PhyloMiner v${VERSION}"
  echo "Date: $(date)"
  echo "Working directory: $(pwd)"
  echo "Command:"
  printf '  %q' "$0" "${ORIGINAL_ARGS[@]}"
  echo
  echo
  echo "Run settings:"
  echo "  Threads      : $THREADS"
  echo "  Force        : $FORCE"
  echo "  Keep temp    : $KEEP_TEMP"
  echo "  Pfam database: ${PFAM_DB:-None}"

  if [[ ${#MOTIF_HMMS[@]} -gt 0 ]]; then
    echo "  Motif HMMs   : ${MOTIF_HMMS[*]}"
  else
    echo "  Motif HMMs   : None"
  fi
} > "$LOG_FILE"

exec > >(tee -a "$LOG_FILE") 2>&1

mkdir -p "$WORKDIR" "$PROT_DIR" "$CDS_DIR"

if [[ "$KEEP_TEMP" == true ]]; then
  LOG_DIR="$BASE_DIR/logfiles"
  mkdir -p "$LOG_DIR"
else
  LOG_DIR="$WORKDIR"
fi

SUMMARY_CSV="$DB_DIR/$(basename "${QUERY%.*}")_hit_summary.csv"
printf 'database,total_protein_seqs,phmmer_hits_total,phmmer_below_threshold,extracted_proteins,hmmsearch_filtered_out,proteins_final,cds_final,cds_output_exists,prot_output_exists\n' > "$SUMMARY_CSV"


QUERY_PROT="$WORKDIR/query.prot.fasta"
if is_nucleotide_fasta "$QUERY"; then
  echo "[INFO] Query detected as nucleotide; translating with transeq frame 1"
  translate_cds_to_protein "$QUERY" "$QUERY_PROT"
else
  echo "[INFO] Query detected as protein"
  cp "$QUERY" "$QUERY_PROT"
  standardize_fasta_headers "$QUERY_PROT" "$QUERY_PROT.tmp"
  mv "$QUERY_PROT.tmp" "$QUERY_PROT"
fi

shopt -s nullglob
DB_FILES=("$DB_DIR"/*)
shopt -u nullglob

for db in "${DB_FILES[@]}"; do
  [[ -f "$db" ]] || continue

  base=$(basename "$db")

  case "$base" in
    hit_summary.csv)
      continue
      ;;
  esac

  if [[ "$base" == *"phylominer_work"* ]]; then
    continue
  fi

  case "$base" in
    *.fa|*.fasta|*.faa|*.fna) ;;
    *) continue ;;
  esac

  name="$base"
  name="${name%.fasta}"
  name="${name%.fa}"
  name="${name%.faa}"
  name="${name%.fna}"

  PHMMER_TXT="$LOG_DIR/${name}_phmmer_output.txt"
  PHMMER_TBL="$LOG_DIR/${name}_phmmer_tbl.txt"
  ID_LIST="$LOG_DIR/${name}_hit_ids.txt"
  PROT_OUT="$PROT_DIR/${name}.homologs.prot.fasta"
  CDS_OUT="$CDS_DIR/${name}.homologs.cds.fasta"
  PROT_DB_TRANSLATED="$WORKDIR/${name}.translated.prot.fasta"
  PROT_DB_CLEAN="$WORKDIR/${name}.prot.clean.fasta"
  CDS_DB_CLEAN="$WORKDIR/${name}.cds.clean.fasta"
  CDS_ID_LIST="$WORKDIR/${name}.cds.ids"

  if [[ "$FORCE" == true ]]; then
    echo "[INFO] Force mode: removing previous outputs for $name"
    safe_rm \
      "$PROT_OUT" \
      "$CDS_OUT" \
      "$PHMMER_TXT" \
      "$PHMMER_TBL" \
      "$ID_LIST"
  fi
  if [[ "$FORCE" == false && -s "$PROT_OUT" ]]; then
    echo "[INFO] Existing output found for $name, skipping"
    continue
  fi

  if is_nucleotide_fasta "$db"; then
    DB_TYPE="cds"
    echo "[INFO] $name appears to be nucleotide"
    standardize_fasta_headers "$db" "$CDS_DB_CLEAN"
    translate_cds_to_protein "$CDS_DB_CLEAN" "$PROT_DB_TRANSLATED"
    PROT_DB="$PROT_DB_TRANSLATED"
    TOTAL_SEQS=$(safe_count "$CDS_DB_CLEAN")
  else
    DB_TYPE="protein"
    echo "[INFO] $name appears to be protein"
    standardize_fasta_headers "$db" "$PROT_DB_CLEAN"
    PROT_DB="$PROT_DB_CLEAN"
    TOTAL_SEQS=$(safe_count "$PROT_DB_CLEAN")
  fi

  if [[ "$FORCE" == false ]]; then
    if [[ "$DB_TYPE" == "protein" && -s "$PROT_OUT" ]]; then
      echo "[INFO] Existing protein output found for $name, skipping"
      continue
    fi
    if [[ "$DB_TYPE" == "cds" && -s "$PROT_OUT" && -s "$CDS_OUT" ]]; then
      echo "[INFO] Existing protein and CDS outputs found for $name, skipping"
      continue
    fi
  fi

  PHMMER_BELOW_THRESHOLD=0
  HMMSEARCH_FILTERED_OUT=0
  CDS_FINAL=0
  CDS_EXISTS="no"
  PROT_EXISTS="no"

  echo "[INFO] Searching $name"
  if ! phmmer --cpu "$THREADS" --noali --tblout "$PHMMER_TBL" "$QUERY_PROT" "$PROT_DB" > "$PHMMER_TXT" 2>&1; then
    echo "WARN: phmmer failed for $name, skipping" >&2
    echo "---- phmmer error for $name ----"
    tail -20 "$PHMMER_TXT" || true
    echo "--------------------------------"
    printf '%s,%s,0,0,0,0,0,0,no,no\n' "$name" "$TOTAL_SEQS" >> "$SUMMARY_CSV"
    continue
  fi

  PHMMER_BELOW_THRESHOLD=$(awk '/inclusion threshold/{f=1;next} f&&NF{c++} END{print c+0}' "$PHMMER_TXT" 2>/dev/null || echo 0)

  if [[ ! -s "$PHMMER_TBL" ]] || [[ $(grep -vc '^#' "$PHMMER_TBL" || true) -eq 0 ]]; then
    echo "[INFO] No phmmer hits for $name"
    printf '%s,%s,0,0,0,0,0,0,no,no\n' "$name" "$TOTAL_SEQS" >> "$SUMMARY_CSV"
    continue
  fi

  extract_ids_from_phmmer "$PHMMER_TXT" "$ID_LIST"
  PHMMER_HITS=$(wc -l < "$ID_LIST" | tr -d ' ')
  echo "[INFO] PHMMER hits: $PHMMER_HITS"

  if [[ "$PHMMER_HITS" -eq 0 ]]; then
    echo "[INFO] No phmmer hits for $name"
    printf '%s,%s,0,0,0,0,0,0,no,no\n' "$name" "$TOTAL_SEQS" >> "$SUMMARY_CSV"
    continue
  fi

  if [[ -s "$ID_LIST" ]]; then
    seqkit grep -n -f "$ID_LIST" "$PROT_DB" > "$PROT_OUT" || true
  else
    : > "$PROT_OUT"
  fi    
    
  EXTRACTED_PROTEINS=$(safe_count "$PROT_OUT")

  if [[ "$EXTRACTED_PROTEINS" -eq 0 ]]; then
    echo "[INFO] Standard protein extraction failed for $name. Trying alternative ID strategies."
    grep "^>" "$PROT_DB" | sed 's/^>//' > "${ID_LIST}.failsafe1.txt"
    seqkit grep -n -r -f "${ID_LIST}.failsafe1.txt" "$PROT_DB" > "$PROT_OUT" || true
    EXTRACTED_PROTEINS=$(safe_count "$PROT_OUT")
    if [[ "$EXTRACTED_PROTEINS" -eq 0 ]]; then
      echo "[INFO] Failsafe 1 failed. Trying pipe-based IDs."
      grep "^>" "$PROT_DB" \
      | sed 's/^>.*|//; s/[[:space:]].*$//' \
      > "${ID_LIST}.failsafe2.txt"
    seqkit grep -n -f "${ID_LIST}.failsafe2.txt" "$PROT_DB" > "$PROT_OUT" || true
    EXTRACTED_PROTEINS=$(safe_count "$PROT_OUT")
    fi
  fi

  if [[ "$EXTRACTED_PROTEINS" -eq 0 ]]; then
    echo "WARN: protein extraction failed for $name" >&2
    printf '%s,%s,0,0,0,0,0,0,no,no\n' "$name" "$TOTAL_SEQS" >> "$SUMMARY_CSV"
    continue
  fi

  PROT_FINAL=$(safe_count "$PROT_OUT")
  BEFORE_HMM=$EXTRACTED_PROTEINS
  PROTEINS_FINAL=$EXTRACTED_PROTEINS
  CDS_FINAL=0
  HMMSEARCH_FILTERED_OUT=0
  PROTEINS_FINAL=$PROT_FINAL
  
  if [[ ${#MOTIF_HMMS[@]} -gt 0 ]]; then
    echo "[INFO] Applying HMM filtering to $name"
    filter_proteins_by_hmms "$PROT_OUT" "$PROT_OUT.filtered" "$WORKDIR" "${MOTIF_HMMS[@]}"
    AFTER_HMM=$(safe_count "$PROT_OUT.filtered")
    echo "[INFO] Before HMM: $BEFORE_HMM After HMM: $AFTER_HMM"
    if [[ "$AFTER_HMM" -gt 0 ]]; then
      HMMSEARCH_FILTERED_OUT=$((BEFORE_HMM - AFTER_HMM))
    else
      echo "[WARN] HMM filtering removed all candidates for $name"
      HMMSEARCH_FILTERED_OUT=0
    fi
    mv "$PROT_OUT.filtered" "$PROT_OUT"
    PROTEINS_FINAL=$AFTER_HMM
  fi

  if [[ "$PROTEINS_FINAL" -eq 0 ]]; then
    : > "$CDS_OUT"
  else

  if [[ "$DB_TYPE" == "cds" ]]; then
    # Strategy 1: exact protein IDs
    grep '^>' "$PROT_OUT" | sed 's/^>//' > "$CDS_ID_LIST"
    seqkit grep -n -f "$CDS_ID_LIST" "$CDS_DB_CLEAN" > "$CDS_OUT" || true
    CDS_FINAL=$(safe_count "$CDS_OUT")
  if [[ "$CDS_FINAL" -eq 0 ]]; then
    echo "[INFO] Exact CDS ID match failed. Trying pipe-based IDs."
    sed -E 's/^.*\|//' "$CDS_ID_LIST" > "${CDS_ID_LIST}.pipe"
    seqkit grep -n -f "${CDS_ID_LIST}.pipe" "$CDS_DB_CLEAN" > "$CDS_OUT" || true
    CDS_FINAL=$(safe_count "$CDS_OUT")
  fi
  # Strategy 3: try protein headers directly
  if [[ "$CDS_FINAL" -eq 0 ]]; then
    echo "[INFO] Pipe-based CDS ID match failed. Trying protein headers."
    grep '^>' "$PROT_OUT" | sed 's/^>//' > "${CDS_ID_LIST}.headers"
    seqkit grep -n -r -f "${CDS_ID_LIST}.headers" "$CDS_DB_CLEAN" > "$CDS_OUT" || true
    CDS_FINAL=$(safe_count "$CDS_OUT")
  fi
  if [[ "$CDS_FINAL" -eq 0 ]]; then
    echo "[WARN] CDS extraction methods failed. Please simplify headers for this database."
  fi

  else
    safe_rm "$CDS_OUT" 2>/dev/null || true
  fi
  fi

  PROTEINS_FINAL=$(safe_count "$PROT_OUT")
  CDS_FINAL=$(safe_count "$CDS_OUT")

  CDS_EXISTS="no"
  PROT_EXISTS="no"

  [[ -s "$CDS_OUT" ]] && CDS_EXISTS="yes"
  [[ -s "$PROT_OUT" ]] && PROT_EXISTS="yes"

  printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' "$name" "$TOTAL_SEQS" "$PHMMER_HITS" "$PHMMER_BELOW_THRESHOLD" "$EXTRACTED_PROTEINS" "$HMMSEARCH_FILTERED_OUT" "$PROTEINS_FINAL" "$CDS_FINAL" "$CDS_EXISTS" "$PROT_EXISTS" >> "$SUMMARY_CSV"

  cleanup_temp_files
done

safe_rm "$QUERY_PROT"
safe_rm "$PROT_DB_TRANSLATED"

echo "[INFO] Finished. Summary written to: $SUMMARY_CSV"
