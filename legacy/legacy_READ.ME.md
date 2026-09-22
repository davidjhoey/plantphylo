
## phmmer scaled up
It may be necessary to run phmmer on many databases. If so, here are some scripts which are helpful with that. These scripts are now in the 'legacy' folder, and formed the basis for `PhyloDig`.
### 1. Run phmmer on a directory of proteomes
If you need to find homologs from many genomes, run phmmer on all of these with a query file of a closely related species. 
Use the file **run_phmmer.sh** (now in the legacy folder) to run phmmer on all proteomes in a directory, using a query file.
If running on ubuntu, need to make this file executable:
```
sed -i 's/\r$//' run_phmmer.sh
chmod +x run_phmmer.sh
```
And then the usage is as below:
```
./run_phmmer.sh query.fasta /path/to/databases
```
There is also a force option, which will overwrite the output file if it exists, -f:
```
./run_phmmer.sh -f query.fasta /path/to/databases
```
### 2. Extracting gene id lists of homologs
Use the file **extract_all_ids.sh** (now in the legacy folder) in order to extract all gene ids which phmmer has identified to be above the inclusion threshold, as homologs of your query protein.
```
sed -i 's/\r$//' extract_all_ids.sh
chmod +x extract_all_ids.sh
./extract_all_ids.sh /path/to/phmmer_outputs
```
There is also a force overwrite here,
```
./extract_all_ids.sh -f /path/to/phmmer_outputs
```
### 3. Extracting fasta files using gene id lists
This next script can be used to convert all the gene lists into fasta sequence files, by searching the proteomes for the sequences in the list. 
```
#!/bin/bash
input_dir="path/to/your/directory"
output_dir="path/to/output/directory"
output_extension="_SEQ.fasta"

for input_file in "$input_dir"/*.fasta; do
    base_name=$(basename "$input_file" .fasta)
    query_file="${base_name}_phmmer_processed.txt"
    output_file="$output_dir/${base_name}${output_extension}"
    
    seqtk subseq "$input_file" "$query_file" > "$output_file"
    
    echo "Processed $input_file -> $output_file"
done
```
The script **extract_all_ids_autoextract.sh** (also now in the legacy folder) combines steps 2. and 3. into one script.
```
./extract_all_ids_autoextract.sh /path/to/phmmer_outputs
```
The above scripts are made obsolete by `PhyloDig`, but may still be useful if you need to do these tasks individually.
