#!/bin/bash

# Usage: ./run_phmmer.sh [-f] query.fasta /path/to/databases
#   -f   Force overwrite of existing outputs

force=false

# Parse options
while getopts ":f" opt; do
  case $opt in
    f) force=true ;;
    \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

query=$1
db_dir=$2

if [[ -z "$query" || -z "$db_dir" ]]; then
    echo "Usage: $0 [-f] query.fasta /path/to/databases"
    exit 1
fi

# Loop over .fa and .fasta files
for db in "$db_dir"/*.fa "$db_dir"/*.fasta; do
    # Skip if no matching files are found
    [[ -e "$db" ]] || continue

    # Extract base name without extension for output naming
    db_name=$(basename "$db")
    db_base="${db_name%.*}"

    # Define output file in the same directory
    out="$db_dir/${db_base}_phmmer_output.txt"

    if [[ -f "$out" && $force == false ]]; then
        echo "Skipping $db (output already exists: $out)"
        continue
    fi

    echo "Running phmmer on $db..."
    phmmer "$query" "$db" > "$out"
done

echo "All searches completed."
