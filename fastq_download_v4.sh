#!/bin/bash
#SBATCH --job-name=download_fastqs
#SBATCH --time=48:00:00               
#SBATCH --mem=2GB   
#SBATCH --cpus-per-task=1

module load seqtk/1.3

# Download fatsq files for associated assemblies
# and subsample to 1 million reads total.

# Need to download enaBrowserTools before running.

LOG_FILE="processed_samples.log"

# Create the log file if it doesn't exist
touch "$LOG_FILE"

# Function to check if a sample has already been processed
is_processed() {
    local sample=$1
    grep -q "^$sample$" "$LOG_FILE"
}

# Function to process a sample
process_sample() {
    local line=$1
    echo "Downloading sample $line"
    base=$(basename "$line" "_FASTA.fasta.gz")

    # Use timeout to limit the download and processing time to 5 minutes (300 seconds)
    if timeout 300 python ~/PATH_TO_PROGRAM/enaBrowserTools-1.7.1/python3/enaGroupGet.py -f fastq -a "$line" -d ~/PATH_TO_OUTPUT/fastqs; then
        # Combine all R1s and R2s:
        find PATH_TO_PROJECT/data/fastqs/${line} -type f -name "*1.fastq.gz" -exec cat {} + > PATH_TO_PROJECT/data/fastqs/${line}/${line}_R1.fq.gz
        find PATH_TO_PROJECT/data/fastqs/${line} -type f -name "*2.fastq.gz" -exec cat {} + > PATH_TO_PROJECT/fastqs/${line}/${line}_R2.fq.gz

        # Check files exist:
        if [[ -f "PATH_TO_PROJECT/data/fastqs/${line}/${line}_R1.fq.gz" && -f "PATH_TO_PROJECT/data/fastqs/${line}/${line}_R2.fq.gz" ]]; then
            # Use timeout to limit the download and processing time to 5 minutes (300 seconds)
            seqtk sample -s100 PATH_TO_PROJECT/data/fastqs/${line}/${line}_R1.fq.gz 500000 > /PATH_TO_PROJECT/data/subsampled_fastqs/${line}_subsampled_R1.fq
            seqtk sample -s100 PATH_TO_PROJECT/data/fastqs/${line}/${line}_R2.fq.gz 500000 > PATH_TO_PROJECT/data/subsampled_fastqs/${line}_subsampled_R2.fq

            rm -r PATH_TO_PROJECT/data/fastqs/${line}
            echo "$line" >> "$LOG_FILE"
        else
            echo "Files for sample $line do not exist, skipping."
        fi
    else
        echo "Processing sample $line took too long and was skipped."
    fi
}


# Read the input file and process each line
while IFS= read -r line; do
    if ! is_processed "$line"; then
        process_sample "$line"
    else
        echo "Sample $line has already been processed, skipping."
    fi
done < n100_sample_accessions_v2.txt
