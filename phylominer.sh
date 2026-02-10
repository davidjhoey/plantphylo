#!/usr/bin/env bash

##################################################################################################################
# phylominer.sh
#
# Protein- and CDS-compatible orthologue extraction pipeline
#
# Default:
#   - Extract FASTA files of protein orthologues from a database using a query protein sequence (phmmer-based)
#   - Capable of querying protein or CDS databases.
#   - Removes intermediate files.
# Options:
#   - Filter using HMM profiles (before or after phmmer)
#   - Keep any intermediate files (translated protein databases, gene id lists.
# 
# Dependencies: EMBOSS transeq, HMMER, seqkit
# IMPORTANT: All input databases should be .fa format in order to be protected from deletion.
# Show help: -h
#
##################################################################################################################

set -euo pipefail

show_help() {
cat << EOF
Usage:
  ./phylominer.sh [options] query.fasta /path/to/databases

Required arguments:
  query.fasta              Protein or CDS FASTA query
  /path/to/databases       Directory containing FASTA databases (MUST be in .fa format)

Options:
  -f                       Force overwrite of existing outputs
  -p                       Keep phmmer output (.txt)
  -i                       Keep list of hit IDs (.ids.txt)
  --keep-prot-db           Keep translated protein databases
  --keep-prot-out          Keep translated protein homolog FASTA
  --threads N              Number of threads for phmmer (default 1)
  --phmmer-only            Legacy mode: phmmer pipeline
                           (no fasta translation or extraction)
  --hmm-first             Forego phmmer search and search using a hmmer profile instead
  -m MOTIF                Require exact motif match (can be repeated)
  --motif-any             Match any motif (default: all motifs required)
  --motif-hmm HMM_ID      Use HMMER motif/domain search (e.g. PF00249)
  --pfam-db PATH          Use local Pfam-A.hmm database (recommended). Need to unpack Pfam-A.hmm using gunzip and dos2unix.
  -h, --help               Show this help message and exit
EOF
}

resolve_pfam_acc() {
  local acc="$1"
  local pfam_db="$2"

  # Exact match first (with optional version)
  local resolved
  resolved=$(grep -m1 -E "^ACC[[:space:]]+${acc}(\.[0-9]+)?$" "$pfam_db" | awk '{print $2}')

  if [[ -n "$resolved" ]]; then
    echo "$resolved"
    return 0
  fi

  # Try any version if non-versioned accession given
  resolved=$(grep -m1 -E "^ACC[[:space:]]+${acc}\.[0-9]+" "$pfam_db" | awk '{print $2}')

  if [[ -n "$resolved" ]]; then
    echo "$resolved"
    return 0
  fi

  return 1
}


############################################
# Defaults
############################################

hmm_first=""
force=false
keep_phmmer=false
keep_ids=false
keep_prot_db=false
keep_prot_out=false
keep_cds_db_clean=false
protein_only=false
threads=1

motifs=()
motif_mode="all"
hmm_ids=()
pfam_db=""

workdir=""


############################################
# Parse options
############################################


while [[ $# -gt 0 ]]; do
  case $1 in
    -f) force=true; shift ;;
    -p) keep_phmmer=true; shift ;;
    -i) keep_ids=true; shift ;;
    --keep-prot-db) keep_prot_db=true; shift ;;
    --keep-prot-out) keep_prot_out=true; shift ;;
    --keep-cds-db-clean) keep_prot_out=true; shift ;;
    --threads) threads="$2"; shift 2 ;;
    --phmmer-only) protein_only=true; shift ;;
    --hmm-first) hmm_first="$2"; shift 2 ;;
    -m) motifs+=("$2"); shift 2 ;;
    --motif-any) motif_mode="any"; shift ;;
    --motif-hmm) hmm_ids+=("$2"); shift 2 ;;
    --pfam-db) pfam_db="$2"; shift 2 ;;
    -h|--help) show_help; exit 0 ;;
    -*) echo "Unknown option: $1" >&2; show_help; exit 1 ;;
    *) break ;;
  esac
done

# REQUIRE TWO POSITIONAL ARGS
if [[ $# -lt 2 ]]; then
  show_help
  exit 1
fi

query=$1
db_dir=$2

if [[ -z "$query" || -z "$db_dir" ]]; then
  show_help
  exit 1
fi

if [[ ! -d "$db_dir" ]]; then
  echo "ERROR: database directory not found: $db_dir"
  exit 1
fi

if [[ ${#motifs[@]} -gt 0 && ${#hmm_ids[@]} -gt 0 ]]; then
  echo "ERROR: cannot use -m and --motif-hmm together."
  exit 1
fi

if [[ ${#motifs[@]} -gt 0 || ${#hmm_ids[@]} -gt 0 ]]; then
  keep_phmmer=false
  keep_ids=false
  keep_prot_db=false
  keep_prot_out=false
fi

# Enforce local HMM database when using motif-HMM
if [[ ${#hmm_ids[@]} -gt 0 && -z "$pfam_db" ]]; then
  echo "ERROR: --motif-hmm requires --pfam-db to be set."
  exit 1
fi

# SAFETY CHECK: never delete original database files
if [[ "$db_dir" == "/" || "$db_dir" == "" ]]; then
  echo "ERROR: database directory invalid. Aborting."
  exit 1
fi

workdir="$db_dir/phmmer_work"
mkdir -p "$workdir"
echo "[INFO] workdir set to: $workdir"

echo "[INFO] pfam_db=${pfam_db}"
echo "[INFO] hmm_ids=${hmm_ids[*]}"

trap 'echo "[FATAL] Error at line $LINENO (exit code $?)"' ERR

if [[ ${#hmm_ids[@]} -gt 0 ]]; then
  resolved_ids=()
  for id in "${hmm_ids[@]}"; do
    resolved_ids+=("$(resolve_pfam_acc "$id" "$pfam_db")")
  done
  hmm_ids=("${resolved_ids[@]}")
fi

if [[ -n "$hmm_first" && "$protein_only" == true ]]; then
  echo "ERROR: --hmm-first cannot be used with --phmmer-only."
  exit 1
fi


############################################
# Pfam HMM extraction helpers
############################################

extract_pfam_hmm() {
  local acc="$1"
  local pfam="$2"
  local out="$3"

  # Try both versioned and non-versioned accession
  local acc_try=("$acc" "${acc%%.*}")

  # If Pfam is pressed, use hmmfetch
  if [[ -f "${pfam}.h3i" || -f "${pfam}.h3m" ]]; then
    for a in "${acc_try[@]}"; do
      if hmmfetch -o "$out" "$pfam" "$a" 2> "$workdir/hmmfetch.err"; then
        if [[ -s "$out" && $(grep -c "^HMMER3/f" "$out") -gt 0 ]]; then
          return 0
        fi
      fi
    done
    echo "[ERROR] Failed to extract valid HMM for $acc from $pfam"
    cat "$workdir/hmmfetch.err"
    rm -f "$out"
    return 1
  else
    echo "[ERROR] Pfam database is not pressed or missing index files (.h3i/.h3m)."
    return 1
  fi
}


############################################
# CSV summary
############################################

summary_csv="$db_dir/hit_summary.csv"

echo "database,total_protein_seqs,phmmer_hits_total,phmmer_below_threshold,extracted_proteins,motif_filtered_out,hmmsearch_filtered_out,proteins_final,cds_final,cds_output_exists,prot_output_exists" > "$summary_csv"

write_summary() {
  echo "$*" >> "$summary_csv"
}

############################################
# Legacy phmmer-only mode
############################################

if [[ $protein_only == true ]]; then
  echo "Running in legacy phmmer-only mode."
  for db in "$db_dir"/*.fa; do
    [[ -e "$db" ]] || continue
    base=$(basename "$db" .fa)
    out="$db_dir/${base}_phmmer_output.txt"
    [[ -f "$out" && $force == false ]] && continue
    phmmer --cpu "$threads" "$query" "$db" > "$out"
  done
  exit 0
fi

############################################
# CDS-aware default mode
############################################

workdir="$db_dir/phmmer_work"
# ensure workdir is within db_dir
if [[ "$workdir" != "$db_dir/"* ]]; then
  echo "ERROR: workdir is outside db_dir. Aborting."
  exit 1
fi
mkdir -p "$workdir"

if grep -v "^>" "$query" | tr -d "\n" | grep -qi "[^ACGTN]"; then
  query_prot="$query"
else
  query_prot="$workdir/query.faa"
  transeq -sequence "$query" -outseq "$query_prot" -frame 1

  # remove the _1 suffix from transeq output headers
  sed -i 's/^\(>[^[:space:]]*\)_1\([[:space:]]\|$\)/\1\2/' "$query_prot"
fi

############################################
# Helper functions
############################################

extract_ids_from_phmmer() {
  awk '!/^#/ {
    id=$1
    # Remove common HMMER prefixes but keep the rest
    sub(/^lcl\|/, "", id)
    sub(/^tr\|/, "", id)
    sub(/^sp\|/, "", id)
    sub(/^ref\|/, "", id)
    sub(/^jgi\|/, "", id)
    sub(/\s+.*/, "", id)
    print id
  }' "$1" | sort -u > "$2"
}


filter_by_motif() {
  local in_fa="$1"
  local out_fa="$2"
  local workdir="$3"

  cp "$in_fa" "$out_fa"

  # if output is empty, nothing to filter
  if [[ ! -s "$out_fa" ]]; then
    echo "[WARN] No sequences to filter in $out_fa"
    return
  fi

  if [[ ${#motifs[@]} -gt 0 ]]; then
    if [[ "$motif_mode" == "all" ]]; then
      for m in "${motifs[@]}"; do
        seqkit grep -s -r -p "$m" "$out_fa" > "${out_fa}.tmp" || true
        mv "${out_fa}.tmp" "$out_fa"
      done
    else
      seqkit grep -s -r -p "$(IFS='|'; echo "${motifs[*]}")" "$out_fa" > "${out_fa}.tmp" || true
      mv "${out_fa}.tmp" "$out_fa"
    fi
    return
  fi

  if [[ ${#hmm_ids[@]} -gt 0 ]]; then
    local ids_current ids_next
    ids_current="$(mktemp)"
    grep -s "^>" "$out_fa" | sed 's/^>//' > "$ids_current" || true

    if [[ ! -s "$ids_current" ]]; then
      rm -f "$ids_current"
      return
    fi

    for id in "${hmm_ids[@]}"; do
      local hmm_file="$workdir/${id}.hmm"
      local tbl="$workdir/${id}.tbl"
      local ids_hit="$workdir/${id}.ids"

      if [[ -n "$pfam_db" ]]; then
        if [[ -d "$pfam_db" ]]; then
          hmm_file="${pfam_db}/${id}.hmm"
          if [[ ! -s "$hmm_file" ]]; then
            echo "[WARN] Missing HMM file: $hmm_file"
            continue
          fi
        else
          extract_pfam_hmm "$id" "$pfam_db" "$hmm_file"
          if [[ $? -ne 0 ]]; then
            echo "[WARN] Failed to extract $id from $pfam_db, skipping"
            continue
          fi
        fi
      else
        echo "[ERROR] --motif-hmm requires --pfam-db to be set."
        exit 1
      fi

      [[ -s "$hmm_file" ]] || {
        echo "ERROR: HMM file $hmm_file is empty"
        exit 1
      }

      grep -q "^HMMER3/f" "$hmm_file" || {
        echo "ERROR: $hmm_file is not a valid HMM file"
        exit 1
      }

      echo "[INFO] Running hmmsearch for ${id} against ${out_fa}"
      if ! hmmsearch --cpu "$threads" --tblout "$tbl" "$hmm_file" "$out_fa" > /dev/null 2>&1; then
        echo "[WARN] hmmsearch failed for ${id}; skipping"
        continue
      fi

      if [[ ! -s "$tbl" ]]; then
        echo "[WARN] hmmsearch produced no table for ${id}; skipping motif filter for this HMM."
        continue
      fi

      awk '!/^#/ {print $1}' "$tbl" > "$ids_hit"

      ids_next="$(mktemp)"
      if [[ ! -s "$ids_current" || ! -s "$ids_hit" ]]; then
        : > "$ids_next"
      else
        sort -u "$ids_current" -o "$ids_current"
        sort -u "$ids_hit" -o "$ids_hit"
        comm -12 "$ids_current" "$ids_hit" > "$ids_next" || true
      fi
      mv "$ids_next" "$ids_current"

      if [[ ! -s "$ids_current" ]]; then
        echo "[WARN] No sequences remain after filtering with ${id}. Exiting."
        rm -f "$ids_current"
        return
      fi
    done

    seqkit grep -n -f "$ids_current" "$out_fa" > "${out_fa}.tmp" || true
    mv "${out_fa}.tmp" "$out_fa"
    rm -f "$ids_current"
  fi
}

cleanup() {
  safe_rm() {
    local file="$1"
    # Protect only the original database FASTA files
    if [[ "$file" == "$db_dir"/*.fa ]]; then
      echo "[WARN] Refusing to delete .fa database file: $file"
      return 0
    fi
    rm -f "$file"
  }

  [[ $keep_phmmer == false ]] && safe_rm "$phmmer_txt" && safe_rm "$phmmer_tbl"
  [[ $keep_ids == false ]] && safe_rm "$id_list" && safe_rm "$id_list_clean" && safe_rm "$id_list_prot"
  [[ $keep_prot_db == false ]] && safe_rm "$prot_db_trans"
  [[ $keep_cds_db_clean == false ]] && safe_rm "$cds_db_clean"
  [[ $keep_prot_out == false ]] && safe_rm "$prot_out"

  # optional alt cleanup
  safe_rm "${cds_db}.alt_clean.fasta"
  safe_rm "${id_list}.alt.txt"
}


verify_motif_filter() {
  local prot_in="$1"
  local prot_out="$2"
  local base="$3"

  if [[ ${#hmm_ids[@]} -gt 0 ]]; then
    local count_in
    local count_out

    count_in=$(grep -c "^>" "$prot_in" 2>/dev/null || echo 0)
    count_out=$(grep -c "^>" "$prot_out" 2>/dev/null || echo 0)

    if [[ "$count_out" -ge "$count_in" ]]; then
      echo "[WARN] Motif filtering did not reduce sequences for $base (count_in=$count_in count_out=$count_out)"
    else
      echo "[INFO] Motif filtering reduced sequences for $base (count_in=$count_in count_out=$count_out)"
    fi
  fi
}


############################################
# Main loop
############################################

if [[ -n "$hmm_first" ]]; then
  if [[ ! -s "$hmm_first" ]]; then
    echo "ERROR: HMM file not found or empty: $hmm_first"
    exit 1
  fi

  if ! grep -q "^HMMER3/f" "$hmm_first"; then
    echo "ERROR: $hmm_first is not a valid HMMER3 profile"
    exit 1
  fi

  echo "[INFO] Primary homology detection: HMM-first (hmmsearch)"
fi


safe_count() {
  [[ -s "$1" ]] && grep -c "^>" "$1" || echo 0
}

is_nucleotide_fasta() {
    local file="$1"
    local non_nuc
    non_nuc=$(grep -v '^>' "$file" | tr -d 'ACGTNacgtnRYSMKWBVDHrysmkwbvdh' | wc -c)
    [[ $non_nuc -eq 0 ]]
}


for cds_db in "$db_dir"/*.fa; do
  [[ -e "$cds_db" ]] || continue

  base=$(basename "$cds_db")
  base="${base%.fa}"
  base="${base%.fasta}"

  prot_db_trans="$workdir/${base}.faa"
  prot_db="$prot_db_trans"

  phmmer_txt="$db_dir/${base}_phmmer_output.txt"
  phmmer_tbl="$db_dir/${base}_phmmer_tbl.txt"
  id_list="$db_dir/${base}_hit_ids.txt"
  id_list_clean="$db_dir/${base}_hit_ids.clean.txt"
  id_list_prot="$db_dir/${base}_hit_ids.prot.txt"

  prot_out="$db_dir/${base}_homologs.prot.fasta"
  cds_out="$db_dir/${base}_homologs.cds.fasta"
  cds_db_clean=""


  if is_nucleotide_fasta "$cds_db"; then
    echo "[INFO] Translating nucleotide database: $cds_db"
    transeq -sequence "$cds_db" -outseq "$prot_db_trans" -frame 1
    # remove ONLY the transeq _1 suffix
    sed -i 's/^\(>[^[:space:]]*\)_1\([[:space:]]\|$\)/\1\2/' "$prot_db_trans"

    # remove description text immediately after translation
    sed -E -i 's/^>(lcl\||tr\||sp\||ref\||jgi\|)?([^[:space:]]+).*/>\2/' "$prot_db_trans"

    cds_db_clean="$workdir/${base}.cds.clean.fasta"
    sed -E 's/^>(lcl\||tr\||sp\||ref\||jgi\|)?([^[:space:]]+).*/>\2/' "$cds_db" > "$cds_db_clean"

    prot_db="$prot_db_trans"
    db_is_nuc=true
  else
    echo "[INFO] Database appears to be protein: $cds_db"
    prot_db="$cds_db"
    sed -i 's/^>\([^[:space:]]*\).*/>\1/' "$prot_db"
    db_is_nuc=false
  fi

    if [[ $force != true ]]; then
  if [[ -s "$prot_out" ]] && [[ "$db_is_nuc" == false || -s "$cds_out" ]]; then
    echo "[INFO] Outputs already exist for $base. Use --force to overwrite."
    continue
  fi
  fi

  total_prot=$(safe_count "$prot_db")

  if [[ -n "$hmm_first" ]]; then
  if ! hmmsearch --cpu "$threads" --noali --tblout "$phmmer_tbl" \
                 "$hmm_first" "$prot_db" > "$phmmer_txt"; then
    echo "[WARN] hmmsearch failed on $base; skipping"
    continue
  fi
  else
  if ! phmmer --cpu "$threads" --noali --tblout "$phmmer_tbl" \
              "$query_prot" "$prot_db" > "$phmmer_txt"; then
    echo "[WARN] phmmer failed on $base; skipping"
    continue
  fi
  fi


  if [[ ! -s "$phmmer_tbl" ]] || [[ $(grep -v "^#" "$phmmer_tbl" | wc -l) -eq 0 ]]; then
    echo "[WARN] phmmer produced no hits for $base"
    continue
  fi

  below_threshold=$(awk '/inclusion threshold/{f=1;next} f&&NF{c++} END{print c+0}' "$phmmer_txt" 2>/dev/null || echo 0)

  extract_ids_from_phmmer "$phmmer_tbl" "$id_list"
  cp "$id_list" "$id_list_prot"

  sed 's/[[:space:]].*$//' "$id_list" > "$id_list_clean"

  phmmer_hits_total=$(wc -l < "$id_list" 2>/dev/null || echo 0)

  seqkit grep -n -f "$id_list_prot" "$prot_db" > "$prot_out" || :
  if [[ ! -s "$prot_out" ]]; then
    echo "[WARN] No protein sequences extracted for $base. Check that phmmer IDs match the database headers."
    continue
  fi

  extracted_proteins=$(safe_count "$prot_out")

  pre_filter=$extracted_proteins
  prot_tmp="${prot_out%.fasta}"
  prot_tmp="${prot_tmp%.fa}.motif.fasta"

  filter_by_motif "$prot_out" "$prot_tmp" "$workdir"
  verify_motif_filter "$prot_out" "$prot_tmp" "$base"

  mv "$prot_tmp" "$prot_out" 2>/dev/null || :
  post_filter=$(safe_count "$prot_out")

  motif_filtered_out=$(( pre_filter - post_filter ))
  hmm_filtered_out=$([[ ${#hmm_ids[@]} -gt 0 ]] && echo "$motif_filtered_out" || echo 0)

  proteins_final=$post_filter
  
  grep "^>" "$prot_out" | sed 's/^>//' > "$id_list_clean"
  #echo "[DEBUG] id_list_clean content:"
  #cat "$id_list_clean"

  if [[ "$db_is_nuc" == true ]]; then
  #echo "[INFO] IDs count: $(wc -l < "$id_list_clean")"
  #echo "[DEBUG] CDS clean file exists: $( [[ -s "$cds_db_clean" ]] && echo yes || echo no )"
  #echo "[DEBUG] First 5 headers in cds_db_clean:"
  grep -m 5 "^>" "$cds_db_clean" | sed 's/^>//'
  echo "[INFO] First 5 IDs in id_list_clean:"
  head -n 5 "$id_list_clean"
    seqkit grep -n -f "$id_list_clean" "$cds_db_clean" > "$cds_out" || :
    echo "[INFO] cds_out count: $(grep -c "^>" "$cds_out" 2>/dev/null || echo 0)"
    cds_final=$(safe_count "$cds_out")
  else
    rm -f "$cds_out" 2>/dev/null || :
    cds_final=0
  fi

  if [[ "$db_is_nuc" == true ]]; then

  seqkit grep -n -f "$id_list_clean" "$cds_db_clean" > "$cds_out" || :
  cds_final=$(safe_count "$cds_out")

  if [[ "$cds_final" -eq 0 && -s "$prot_out" ]]; then
    echo "[INFO] CDS extraction failed for $base. Retrying with alternative ID strategy."

    # -----------------------------
    # FAILSAFE 1: use raw protein headers
    # -----------------------------
    grep "^>" "$prot_out" | sed 's/^>//' > "${id_list}.failsafe.txt"

    # use original CDS (no cleaning)
    cp "$cds_db" "${cds_db}.failsafe.fasta"

    # Use regex matching to match IDs embedded in complex headers
    seqkit grep -n -r -f "${id_list}.failsafe.txt" "${cds_db}.failsafe.fasta" > "$cds_out" || :
    cds_final=$(safe_count "$cds_out")
    safe_rm "${cds_db}.failsafe.fasta"

    # If still zero, attempt FAILSAFE 2
    if [[ "$cds_final" -eq 0 ]]; then
      echo "[INFO] Failsafe 1 did not work. Trying alternate ID extraction."

      # -----------------------------
      # FAILSAFE 2: attempt to extract ID after the last pipe
      # -----------------------------
      grep "^>" "$prot_out" | sed 's/^>.*|//; s/[[:space:]].*$//' > "${id_list}.failsafe2.txt"

      seqkit grep -n -f "${id_list}.failsafe2.txt" "${cds_db}.failsafe.fasta" > "$cds_out" || :
      cds_final=$(safe_count "$cds_out")
    fi

    if [[ "$keep_ids" == false ]]; then
    safe_rm "${id_list}.failsafe.txt"
    safe_rm "${id_list}.failsafe2.txt"
    fi

    # If still zero, final message
    if [[ "$cds_final" -eq 0 ]]; then
      echo "[WARN] CDS extraction failed for $base after both alternative strategies. Please simplify gene identifiers for this database."
    fi
    
  fi
  fi

  write_summary \
    "$base,$total_prot,$phmmer_hits_total,$below_threshold,$extracted_proteins,$motif_filtered_out,$hmm_filtered_out,$proteins_final,$cds_final,$([[ -s "$cds_out" ]] && echo yes || echo no),$([[ -s "$prot_out" ]] && echo yes || echo no)"

  cleanup
  rm -f "$workdir"/*.tbl "$workdir"/*.ids
done

set -e

echo "All searches completed."
