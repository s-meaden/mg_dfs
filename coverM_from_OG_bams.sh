#!/bin/bash
#SBATCH --job-name=coverage_extract
#SBATCH --time=02:00:00               
#SBATCH --mem=5GB   
#SBATCH --cpus-per-task=1   

module load SAMtools/1.17-GCC-12.2.0
module load Miniconda3/23.5.2-0
source activate coverm


# Get coverage of defence systems identified by padloc

for i in mapping/bam_files/*_sorted.bam
do
	# Check if the BAM file exists
	if [[ ! -f "$i" ]]; then
		echo "Skipping $i: BAM file not found."
		continue
	fi

	base=$(basename ${i} "_sorted.bam")
	echo "Processing sample: ${base}"

	# Check if the expected output files already exist, skip if they do
	if [[ -f def_coverm_outputs/${base}_def_count.tsv ]]; then
		echo "Skipping ${base}: Output already exists."
		continue
	fi

	# Check if the required padloc file exists
	if [[ ! -f padloc_outputs/${base}_FASTA_pad_out/${base}_FASTA.fasta_padloc.gff ]]; then
		echo "Skipping ${base}: padloc GFF file not found."
		continue
	fi

	# Run coverm commands
	coverm contig --bam-files $i --methods count --output-file def_coverm_outputs/${base}_full_count_tmp.tsv
	coverm contig --bam-files $i --methods metabat --output-file def_coverm_outputs/${base}_full_cov_tmp.tsv

	# Ensure cut command can find the file
	cut -f1 padloc_outputs/${base}_FASTA_pad_out/${base}_FASTA.fasta_padloc.gff > def_coverm_outputs/${base}_tmp_contig_list.txt

	# Filter the contigs and output the final results
	grep -Ff def_coverm_outputs/${base}_tmp_contig_list.txt def_coverm_outputs/${base}_full_count_tmp.tsv | awk '$2 != 0' > def_coverm_outputs/${base}_def_count.tsv
	grep -Ff def_coverm_outputs/${base}_tmp_contig_list.txt def_coverm_outputs/${base}_full_cov_tmp.tsv | awk '$3 != 0' > def_coverm_outputs/${base}_def_cov.tsv

	echo "Completed processing for ${base}."
done



# Do the same for the viral coverage:
# Must be run after above job completes as relies on the full count CoverM file produced above.


for i in mapping/bam_files/*_sorted.bam
do
    echo "Processing $i"
    base=$(basename "${i}" "_sorted.bam")
    echo "Base name: ${base}"

    viral_file="genomad_outputs/${base}_FASTA_pad_out/${base}_FASTA_summary/${base}_FASTA_virus_summary.tsv"

    # Check if viral summary file exists
    if [[ ! -f "${viral_file}" ]]; then
        echo "Skipping ${base}: viral summary file not found (${viral_file})."
        continue
    fi

    # Skip processing if output files already exist
    if [[ -f "viral_coverm_outputs2/${base}_vir_count.tsv" && -f "viral_coverm_outputs2/${base}_vir_cov.tsv" ]]; then
        echo "Skipping ${base}: Output files already exist."
        continue
    fi

    # Process viral file and generate contig list
    awk -F '\t' '$2 > 999 && $9 > 0' "${viral_file}" | cut -f1 | sed 's/|provirus_[0-9]*_[0-9]*//g' | tail -n +2 > viral_coverm_outputs2/${base}_tmp_contig_list.txt

    # Ensure tmp_contig_list.txt is not empty before proceeding
    if [[ ! -s viral_coverm_outputs2/${base}_tmp_contig_list.txt ]]; then
        echo "No valid contigs found for ${base}. Skipping."
        continue
    fi

    # Perform filtering on count and coverage files
    grep -Ff viral_coverm_outputs2/${base}_tmp_contig_list.txt def_coverm_outputs/${base}_full_count_tmp.tsv | awk '$2 != 0' > viral_coverm_outputs2/${base}_vir_count.tsv
    grep -Ff viral_coverm_outputs2/${base}_tmp_contig_list.txt def_coverm_outputs/${base}_full_cov_tmp.tsv | awk '$3 != 0' > viral_coverm_outputs2/${base}_vir_cov.tsv

    echo "Completed processing for ${base}."
done


# And now the same for the CRISPR array carrying contigs.

for i in ../data/mapping/bam_files/*_sorted.bam
do
    echo "Processing $i"
    base=$(basename "${i}" "_sorted.bam")
    echo "Base name: ${base}"

    array_file="../data/metaCRT_outputs/${base}_FASTA_metaCRT_out.txt"

    # Check if array summary file exists
    if [[ ! -f "${array_file}" ]]; then
        echo "Skipping ${base}: array summary file not found (${array_file})."
        continue
    fi
	
    # Skip processing if output files already exist
    if [[ -f "../data/metaCRT_outputs2/${base}_array_count.tsv" && -f "../data/metaCRT_outputs2/${base}_array_cov.tsv" ]]; then
        echo "Skipping ${base}: Output files already exist."
        continue
    fi

    # Process metaCRT file and generate contig list
    grep "SEQ" "${array_file}" | sed -e 's/SEQ: //g' > ../data/metaCRT_outputs2/${base}_tmp_contig_list.txt
   
    # Ensure tmp_contig_list.txt is not empty before proceeding
    if [[ ! -s ../data/metaCRT_outputs2/${base}_tmp_contig_list.txt ]]; then
        echo "No valid contigs found for ${base}. Skipping."
        continue
    fi

    # Perform filtering on count and coverage files
    grep -Ff ../data/metaCRT_outputs2/${base}_tmp_contig_list.txt ../data/def_coverm_outputs/${base}_full_count_tmp.tsv | awk '$2 != 0' > ../data/metaCRT_outputs2/${base}_array_count.tsv
    grep -Ff ../data/metaCRT_outputs2/${base}_tmp_contig_list.txt ../data/def_coverm_outputs/${base}_full_cov_tmp.tsv | awk '$3 != 0' > ../data/metaCRT_outputs2/${base}_array_cov.tsv

    echo "Completed processing for ${base}."
done








