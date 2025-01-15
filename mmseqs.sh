#!/bin/bash
#SBATCH --job-name=mmseqs               # job name (shows up in the queue)
#SBATCH --time=24:00:00                                        # Walltime (HH:MM:SS)
#SBATCH --mem=450G 
#SBATCH --cpus-per-task=24                                                             # memory/cpu (in MB)

# Classify contigs with MMSeqs2
# Run Mseqs taxonomy on the contigs identified as carrying Tiamat defence

# Load required module
module load MMseqs2/13-45111-gompic-2020b

mmseqs createdb tiamat_sequences.fna tiamat_query_DB

mmseqs taxonomy tiamat_query_DB PATH_TO_NR_DATABASE/nr tiamat_mmseqs_taxonomy.txt tiamat_temp --tax-lineage 1 --orf-filter 1 --threads 24 --split-memory-limit 500G 

mmseqs createtsv tiamat_query_DB tiamat_mmseqs_taxonomy.txt tiamat_mmseqs_taxonomy_final.txt

# Repeat for argonautes


mmseqs createdb argonaute_sequences.fna argonaute_query_DB

mmseqs taxonomy argonaute_query_DB PATH_TO_NR_DATABASE/nr argonaute_mmseqs_taxonomy.txt argonaute_temp --tax-lineage 1 --orf-filter 1 --threads 24 --split-memory-limit 500G 

mmseqs createtsv argonaute_query_DB argonaute_mmseqs_taxonomy.txt argonaute_mmseqs_taxonomy_final.txt

