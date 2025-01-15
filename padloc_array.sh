#!/bin/bash
#SBATCH --job-name=padloc_processing
#SBATCH --output=logs/padloc_processing_%A_%a.out
#SBATCH --error=logs/padloc_processing_%A_%a.err
#SBATCH --array=1-1595 # adjust for number of input files
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

module load Miniconda3/23.5.2-0
source activate padloc

# Define the input directory, output directory, and the file list
INPUT_DIR="PATH_TO_PROJECT/data/assemblies"
OUTPUT_DIR="PATH_TO_PROJECT/data/padloc_outputs"
FILE_LIST="PATH_TO_PROJECT/data/assemblies_list.txt"
FAILED_LOG="PATH_TO_PROJECT/data/failed_files.log"
MASTER_STATUS_LOG="PATH_TO_PROJECT/data/master_status.log"


# Get the specific input file based on the SLURM_ARRAY_TASK_ID
INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILE_LIST)

# Define the output directory name (matching the input file name without extension)
OUTPUT_DIR_PER_SAMPLE="${OUTPUT_DIR}/$(basename ${INPUT_FILE} .fasta.gz)_pad_out"

# Print the input and output file names for logging
echo "Processing file: ${INPUT_FILE}"
echo "Output directory: ${OUTPUT_DIR_PER_SAMPLE}"

# Function to log failures
log_failure() {
    echo "${INPUT_FILE}" >> "${FAILED_LOG}"
}

# Trap any exit signals indicating failure (non-zero exit code)
trap 'log_failure' ERR

# Gunzip the input assembly file
gunzip -c "${INPUT_DIR}/${INPUT_FILE}" > "${INPUT_DIR}/$(basename ${INPUT_FILE} .gz)"

# Create the output directory for the current sample
mkdir -p "${OUTPUT_DIR_PER_SAMPLE}"

# Run padloc with the specified input file and output directory
if ! padloc --fna "${INPUT_DIR}/$(basename ${INPUT_FILE} .gz)" --outdir "${OUTPUT_DIR_PER_SAMPLE}" --cpu 4; then
    log_failure
    exit 1
fi

# Remove the unzipped input file to save space
rm "${INPUT_DIR}/$(basename ${INPUT_FILE} .gz)"


# Tidy up outputs

# Gzip the output directory
tar -czf "${OUTPUT_DIR_PER_SAMPLE}.tar.gz" -C "${OUTPUT_DIR}" "$(basename ${OUTPUT_DIR_PER_SAMPLE})"

# Remove the uncompressed output directory
rm -rf "${OUTPUT_DIR_PER_SAMPLE}"


# Search .out file for "State: COMPLETED (exit code 0)"
# Check if the job completed successfully

JOB_OUT_FILE="logs/padloc_processing_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
if grep -q "State: COMPLETED (exit code 0)" "$JOB_OUT_FILE"; then
    STATUS="SUCCESS"
else
    STATUS="FAILURE"
    log_failure
fi

# Log the status to the master status log
echo "${INPUT_FILE}: ${STATUS}" >> "${MASTER_STATUS_LOG}"












