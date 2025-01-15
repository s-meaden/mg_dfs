#!/bin/bash
#SBATCH --job-name=genomad_processing
#SBATCH --output=logs/genomad_processing_%A_%a.out
#SBATCH --error=logs/genomad_processing_%A_%a.err
#SBATCH --array=1-2 # Update to match number of input files
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00

module load Miniconda3/23.5.2-0
source activate genomad

# Define the input directory, output directory, and the file list
INPUT_DIR="PATH_TO_PROJECT/data/assemblies"
OUTPUT_DIR="PATH_TO_PROJECT//data/genomad_outputs"
FILE_LIST="PATH_TO_PROJECT/data/assemblies_list.txt"
FAILED_LOG="PATH_TO_PROJECT//genomad_failed_files.log"


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


# Create the output directory for the current sample
mkdir -p "${OUTPUT_DIR_PER_SAMPLE}"

# Run padloc with the specified input file and output directory
if ! genomad end-to-end --cleanup --threads 4 "${INPUT_DIR}/$(basename ${INPUT_FILE})" "${OUTPUT_DIR_PER_SAMPLE}" /users/kvw508/scratch/databases/genomad_db; then
    log_failure
    exit 1
fi