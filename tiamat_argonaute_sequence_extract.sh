#!/usr/bin/env bash

# Search Padloc results for Tiamat and Argonaute carrying contigs.


for i in `ls PATH_TO_PADLOC_RESULTS/*padloc.csv`
do
    # Print the file name to keep track of which file is being processed
    echo "Processing file: $i"
    
    # Run grep search to find lines containing "Tiamat" and include filename in result
    grep -H "argonaute" $i | while IFS= read -r line; do
        # Format output: Filename, matching line
        echo "$i,$line" >> search_outputs/argonaute_grep_results.csv
    done
done


for i in `ls PATH_TO_PADLOC_RESULTS/*padloc.csv`
do
    # Print the file name to keep track of which file is being processed
    echo "Processing file: $i"
    
    # Run grep search to find lines containing "Tiamat" and include filename in result
    grep -H "Tiamat" $i | while IFS= read -r line; do
        # Format output: Filename, matching line
        echo "$i,$line" >> search_outputs/tiamat_grep_results.csv
    done
done


# Extract lists of sequences for Tiamat and argonaute in R

argonaute_tiamat_seq_extract.R

# tarball output files and scp'd to server.

# Extract identified contigs from assemblies:

module load Biopython/1.83-foss-2023a

# Loop through all samples and save outputs:

for i in `ls ./seqs_lists/*.txt`
do
	echo $i
	base=$(basename ${i} "_tiamat_seqs.txt")
	echo ${base}
	# Update python script to handle gzipped inputs
	python extract_gzipped_seqs.py ../assemblies/${base}_FASTA.fasta.gz $i >> tiamat_sequences.fna
done

for i in `ls ./seqs_lists/*.txt`
do
	echo $i
	base=$(basename ${i} "_argonaute_seqs.txt")
	echo ${base}
	# Update python script to handle gzipped inputs
	python extract_gzipped_seqs.py ../assemblies/${base}_FASTA.fasta.gz $i >> tiamat_sequences.fna
done


# Taxonomically assign these contigs using MMseqs2 with the NCBI NR database:
# mmseqs.sh script 

