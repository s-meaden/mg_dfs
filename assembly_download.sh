#!/bin/bash
#SBATCH --job-name=fq_download_sample        # job name (shows up in the queue)
#SBATCH --time=48:00:00               
#SBATCH --mem=10GB   
#SBATCH --cpus-per-task=1

# Download assemblies based on MGnifyR queries of metadata

completed_downloads_file=PATH_TO_PROJECT/completed_downloads.txt
input_file=PATH_TO_PROJECT/data/n100_URLs.txt
output_dir=PATH_TO_PROJECT/data/assemblies

# Create the completed downloads file if it doesn't exist
touch "$completed_downloads_file"

# Read each URL from the input file
while IFS= read -r line
do
  # Check if the URL is already in the completed downloads file
  if ! grep -Fxq "$line" "$completed_downloads_file"
  then
    # If not, download the file
    echo "Downloading: $line"
    base=$(basename "$line")
    curl -o "$base" "$line"

    # Check if the download was successful
    if [ $? -eq 0 ]; then
      # If successful, add the URL to the completed downloads file
      echo "$line" >> "$completed_downloads_file"
      echo "Downloaded and recorded: $line"
    else
      echo "Failed to download: $line"
    fi
  else
    echo "Already downloaded: $line"
  fi
done < "$input_file"

# Move the downloaded files to the output directory
mv ./*.fasta.gz "$output_dir"
