
library(tidyverse)

# Collect contigs carrying Tiamat or argonaute hits.
# Loop through padloc results CSV files and grep search for "Tiamat" or "Ago"

# Read in tiamat grep search. Extract protein names for extraction.

arg<-read.csv("PATH_TO_RESULTS/argonaute_grep_results.csv")

# Make a master file for each sample with list on target contigs. 

master_sample_list<-arg %>%
  select( -best.hits, -all.domains) %>%
  mutate( system.number = stringi::stri_extract_last_regex(system.number, "[0-9]+")) %>%
  mutate( sample = gsub("_out/", "", stringi::stri_extract_last_regex(filename, "_out/.+_padloc.csv"))) %>%
  mutate( sample = gsub("_FASTA.fasta_padloc.csv", "", sample)) %>%
  select( sample, seqid)

# Divide into files for writing out unique lists (per sample)
split_data <- split(master_sample_list, master_sample_list$sample)

# Write each subset to a new file

setwd("PATH_TO_OUTPUT_DIR/")

for (sample_id in names(split_data)) {
  # Define the output file name (e.g., "sample_ERZ11466566_argonaute_seqs.txt")
  output_file <- paste0(sample_id, "_argonaute_seqs", ".txt")
  
  # Extract only the 'sequence ID' column (assuming it is named 'sequence_id')
  sequence_ids <- split_data[[sample_id]][["seqid"]]
  
  # Write the sequence IDs to a new file as a single column
  write.table(sequence_ids, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Optional: Print a message to confirm each file is written
  cat("Written file:", output_file, "\n")
}


######## Repeat for Tiamat sequences
tiam<-read.csv("PATH_TO_RESULTS//tiamat_grep_results.csv")

# Again make master list
master_sample_list<-tiam %>%
  select( -best.hits, -all.domains) %>%
  mutate( system.number = stringi::stri_extract_last_regex(system.number, "[0-9]+")) %>%
  mutate( sample = gsub("_out/", "", stringi::stri_extract_last_regex(filename, "_out/.+_padloc.csv"))) %>%
  mutate( sample = gsub("_FASTA.fasta_padloc.csv", "", sample)) %>%
  select( sample, seqid)

# Split dataset and write uniqe file per sample:

split_data <- split(master_sample_list, master_sample_list$sample)

# Step 3: Write each subset to a new file

setwd("PATH_TO_OUTPUT_FOLDER/")

for (sample_id in names(split_data)) {
  # Define the output file name (e.g., "sample_ERZ11466566_tiamat_seqs.txt")
  output_file <- paste0(sample_id, "_tiamat_seqs", ".txt")
  
  # Extract only the 'sequence ID' column (assuming it is named 'sequence_id')
  sequence_ids <- split_data[[sample_id]][["seqid"]]
  
  # Write the sequence IDs to a new file as a single column
  write.table(sequence_ids, file = output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Optional: Print a message to confirm each file is written
  cat("Written file:", output_file, "\n")
}
