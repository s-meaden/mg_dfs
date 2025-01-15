
library(dplyr)

# Command line get padloc files into a single directory.
#mkdir padloc_tmp_outputs/
#cp padloc_outputs/*_FASTA_pad_out/*.csv padloc_tmp_outputs/


# Function to process a pair of files (defense read count file and padloc file)
process_pair <- function(def_count_file, padloc_file, sample_name) {
  
  # Check if the defense count file is empty
  if (file.info(def_count_file)$size == 0) {
    warning(paste("Skipping", def_count_file, "- Defense count file is empty"))
    return(NULL)
  }
  
  # Read the defense count file (file1)
  file1 <- read.table(def_count_file, header = TRUE, sep = "\t")
  
  colnames(file1)<-c('Contig', 'read_count')
  
  # Check if file1 is empty (i.e., no rows of data)
  if (nrow(file1) == 0) {
    warning(paste("Skipping", def_count_file, "- Defense count file is empty"))
    return(NULL)
  }
  
  # Read the corresponding padloc file (file2)
  file2 <- read.csv(padloc_file)
  
  # Check if file2 is empty (i.e., no rows of data)
  if (nrow(file2) == 0) {
    warning(paste("Skipping", padloc_file, "- Padloc file is empty"))
    return(NULL)
  }
  
  # Process file1: rename and select columns
  file1 <- file1 %>%
    rename(seqid = Contig)
  
  # Process file2: select and deduplicate relevant columns
  file2 <- file2 %>%
    select(-all.domains, -best.hits, -target.description) %>%
    select(seqid, system) %>%
    unique()
  
  # Combine the data by joining file1 and file2
  final_dataset <- file2 %>%
    left_join(file1, by = "seqid") %>%
    mutate(sample = sample_name) %>%
    rename(read_count = 3)
  
  return(final_dataset)
}

# List of defense read recruitment files and corresponding padloc files
def_count_files <- list.files(path = "PATH_TO_COVERM_OUTPUTS/def_coverm_outputs/", 
                              pattern = "_def_count.tsv", full.names = TRUE)
padloc_files <- list.files(path = "PATH_TO_PADLOC_OUTPUTS/padloc_tmp_outputs/", 
                           pattern = "_padloc.csv", full.names = TRUE)

# Initialize an empty list to store results
results_list <- list()

# Process each file pair
for (def_count_file in def_count_files) {
  
  # Extract the sample name from the file name
  sample_name <- gsub("_def_count.tsv", "", basename(def_count_file))
  
  print(sample_name)
  # Find the corresponding padloc file
  padloc_file <- padloc_files[grep(sample_name, padloc_files)]
  
  print(padloc_file)
  # Ensure that there's a matching padloc file
  if (length(padloc_file) == 1) {
    # Process the file pair and store the result
    result <- process_pair(def_count_file, padloc_file, sample_name)
    
    # Only store non-null results (i.e., files that were not skipped)
    if (!is.null(result)) {
      results_list[[sample_name]] <- result
    }
  } else {
    warning(paste("No matching padloc file found for", sample_name))
  }
}

# Combine all non-null results into a master dataset
if (length(results_list) > 0) {
  master_dataset <- do.call(rbind, results_list)
} else {
  warning("No valid datasets processed.")
}


# Display or save the master dataset
#print(master_dataset)

#print(unique(master_dataset$sample))
# Optionally, write to CSV
write.csv(master_dataset, "PATH_TO_PROJECT/master_dataset_defence_read_count2.csv", row.names = FALSE, quote = FALSE)