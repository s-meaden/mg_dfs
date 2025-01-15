#!/bin/bash
#SBATCH --job-name=kraken_array
#SBATCH --time=05:00:00               
#SBATCH --mem=100GB   
#SBATCH --cpus-per-task=1
#SBATCH --array=0-130 # adjust to number of input samples

module load Kraken2/2.1.2-gompi-2021a

# Define variables
#INPUT_FILE_LIST="n100_sample_accessions_v2.txt"  # File containing the list of input filenames
INPUT_FILE_LIST="remaining_samples.txt" 
RESULTS_DIR="./results"        # Directory for output results
LOG_FILE="./processing.log"    # Log file to track completed files
NUMERIC_VALUE=0.01             # Numeric threshold for filtering

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# Create log file if it doesn't exist
touch "$LOG_FILE"

# Get the total number of lines (input files) in INPUT_FILE_LIST
TOTAL_FILES=$(wc -l < "$INPUT_FILE_LIST")

# Ensure the array size is correct. If SLURM_ARRAY_TASK_ID exceeds the number of input files, exit.
if [ "$SLURM_ARRAY_TASK_ID" -ge "$TOTAL_FILES" ]; then
    echo "SLURM_ARRAY_TASK_ID exceeds total files count. Exiting."
    exit 1
fi

# Extract the line corresponding to this array task from INPUT_FILE_LIST
INPUT_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$INPUT_FILE_LIST")

# Check if this file has already been processed
if grep -q "$INPUT_FILE" "$LOG_FILE"; then
    echo "Skipping already processed file: $INPUT_FILE"
    exit 0
fi

# Run the Kraken2 command and log output
echo "Processing file: $INPUT_FILE"
if kraken2 --db ~/PATH_TO_KRAKEN_DATABASE/databases/kraken2 --threads 2 --report "$RESULTS_DIR/${INPUT_FILE}_kraken_report.txt" --report-minimizer-data --minimum-hit-groups 2 --use-names PATH_TO_FASTQ_DATA/${INPUT_FILE}_subsampled_R1.fq > tmp_kraken_output.txt; then
    # If Kraken2 command succeeds, process the report file
    awk -v sample_name="$INPUT_FILE" '$6 ~ /^[URDPCOFG]$/ && $1 >= '"$NUMERIC_VALUE"' {print $0, sample_name}' "$RESULTS_DIR/${INPUT_FILE}_kraken_report.txt" >> "$RESULTS_DIR/kraken_report_master_file.txt"
    rm "$RESULTS_DIR/${INPUT_FILE}_kraken_report.txt"

    # Log the processed file
    echo "$INPUT_FILE" >> "$LOG_FILE"
    echo "Processed and logged file: $INPUT_FILE"
else
    echo "Error processing file: $INPUT_FILE" >> "$LOG_FILE"
    echo "Failed to process file: $INPUT_FILE. Check the log for details."
fi

echo "Processing completed for file: $INPUT_FILE"


