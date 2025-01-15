#!/usr/bin/env Rscript

# Defence system classifier (chromosomal vs MGE)

# Load necessary libraries
library(dplyr)

########### Read in appropriate files for each sample

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure there are enough arguments
if (length(args) != 4) {
  stop("You must provide 4 arguments: <sample_id> <defence_carrying_contigs> <class_file> <pro_file>")
}

# Assign arguments to variables
sample_id <- args[1]
defence_carrying_contigs_file <- args[2]
class_file <- args[3]
pro_file <- args[4]

# Read in the data
defence_carrying_contigs <- read.csv(defence_carrying_contigs_file, header = TRUE)
class <- read.table(class_file, header = TRUE, sep = "\t")

# Check if 'pro_file' exists
if (file.exists(pro_file)) {
  # If the file exists, read it
  pro <- read.table(pro_file, header = TRUE, sep = "\t")
} else {
  # If the file does not exist, create an empty data frame
  message("Warning: The pro_file does not exist. Proceeding without it.")
  pro <- data.frame()  # Empty data frame to allow the script to continue
}

###### Debugging statements (optional)
#print("Debug stop 1")
#print(args[1])
#print(class(sample_id))
#sample_id<-as.vector(sample_id)
#print(head(defence_carrying_contigs))

# Filter the defence_carrying_contigs for the specific sample_id
tmp <- defence_carrying_contigs %>%
  select(seqid, system, start, end, Folder) %>%
  mutate(sample = gsub("_FASTA_pad_out", "", Folder)) %>%
  select(-Folder) %>%
  filter(sample == sample_id)  # Use sample_id directly

print("Debug stop 2")

######### Run the required searches / merges:

# Select max score from genomad:

tmp2 <- tmp %>%
  rename(seq_name = seqid) %>%
  left_join(class, by = "seq_name") %>%
  mutate(general_class = case_when(
    chromosome_score > plasmid_score & chromosome_score > virus_score ~ "chromosome",
    virus_score > plasmid_score & virus_score > chromosome_score ~ "phage",
    plasmid_score > virus_score & plasmid_score > chromosome_score ~ "plasmid")) %>%
  select(-chromosome_score, -plasmid_score, -virus_score)

# Add in prophage info if the pro_file exists
if (nrow(pro) > 0) {
  # Only perform these steps if the 'pro' data frame has rows (i.e., pro_file exists)
  
  # First extract original contigID and parse start stop positions of prophage
  pro <- pro %>%
    mutate(prophage_start = as.numeric(gsub("_[0-9]+", "", stringi::stri_extract_last_regex(seq_name, "[0-9]+_[0-9]+")))) %>%
    mutate(prophage_stop = as.numeric(gsub("_", "", stringi::stri_extract_last_regex(seq_name, "_[0-9]+")))) %>%
    mutate(seq_name = gsub("\\|provirus_[0-9]+_[0-9]+", "", seq_name)) %>%
    select(seq_name, prophage_start, prophage_stop)
  
  # Combine with contig defence predictions
  classifications <- tmp2 %>%
    left_join(pro, by = "seq_name", relationship = "many-to-many") %>%
    arrange(-prophage_start) %>%
    mutate(prophage = ifelse(start > prophage_start & end < prophage_stop, "prophage", general_class)) %>%
    mutate(classification = ifelse(is.na(prophage), general_class, prophage)) %>%
    select(seq_name, system, sample, classification)
} else {
  # If no prophage data, just use the general classification
  classifications <- tmp2 %>%
    mutate(classification = general_class) %>%
    select(seq_name, system, sample, classification)
}

###########

# Save the results or display them
output_file <- paste0("PATH_TO_PROJECT/contig_classifications/", sample_id, ".tsv")
# Write to output file
write.table(classifications, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
