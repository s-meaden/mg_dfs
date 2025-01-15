
library(MGnifyR)
library(dplyr)
library(readr)

#### Script to collect metadata from MGnify query.
### Works by splitting the accessions list (35000)
### into smaller groups then recombining.

mgclnt <- MgnifyClient(usecache = T, cache_dir = 'PATH_TO_PROJECT/tmp/MGnify_cache/')

df<-read.csv("PATH_TO_DATA/mgnify_assemblies_query_14062024.csv")

# Set working dir so checkpoints folder in correct location;
setwd( "PATH_TO_PROJECT/pipeline")

# Test run
#tmp_list<- df %>%
#  select( accession) %>%
#  head( n = 100)

# In anger
tmp_list<- df %>%
    select( accession)

# Define chunk size. Use 100 for full df. Will add a checkpoint every 100 so job can be re-started
chunk_size <- 100

# Extract the vector that getMetadata() MGnifyR script will use
vector_to_process <- tmp_list$accession

# Split the vector into chunks
split_vector <- split(vector_to_process, ceiling(seq_along(vector_to_process) / chunk_size))

# Directory to store checkpoints
checkpoint_dir <- "checkpoints"
if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir)
}

# Function to save checkpoint
save_checkpoint <- function(chunk, index) {
  saveRDS(chunk, file = file.path(checkpoint_dir, paste0("chunk_", index, ".rds")))
}

# Function to load checkpoint
load_checkpoint <- function(index) {
  file_path <- file.path(checkpoint_dir, paste0("chunk_", index, ".rds"))
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    return(NULL)
  }
}

# Function to find the last completed chunk
find_last_completed_chunk <- function() {
  files <- list.files(checkpoint_dir, pattern = "chunk_.*\\.rds", full.names = TRUE)
  if (length(files) == 0) {
    return(0)
  } else {
    indices <- as.numeric(sub(".*chunk_([0-9]+)\\.rds", "\\1", files))
    return(max(indices))
  }
}

# Find the last completed chunk
last_completed_chunk <- find_last_completed_chunk()

# Process chunks
processed_chunks <- list()
for (i in seq_along(split_vector)) {
  if (i <= last_completed_chunk) {
    cat("Skipping chunk", i, "as it has already been processed.\n")
    checkpoint <- load_checkpoint(i)
    processed_chunks[[i]] <- checkpoint
    next
  }
  
  cat("Processing chunk", i, "of", length(split_vector), "\n")
  processed_chunk <- tryCatch({
    getMetadata(mgclnt, split_vector[[i]], usecache = TRUE)
  }, error = function(e) {
    cat("Error processing chunk", i, ": ", e$message, "\n")
    NULL
  })
  
  if (!is.null(processed_chunk)) {
    # Save checkpoint
    save_checkpoint(processed_chunk, i)
    processed_chunks[[i]] <- processed_chunk
  } else {
    cat("Skipping chunk", i, "due to error.\n")
  }
}

# Combine the processed chunks back into a single dataframe
final_df <- bind_rows(processed_chunks)


# Print the final dataframe (optional)
print(final_df)


##### Read in RDS objects:

tst<-readRDS("./checkpoints/chunk_2.rds")

tst %>%
  head()

names(tst)

# Keep all columns that contain biome information:

file_paths <- paste("./checkpoints/", list.files("./checkpoints/", pattern = "rds"), sep = "")

get_headers_from_rds_files <- function(file_paths) {
  # Initialize an empty list to store headers
  headers_list <- list()
  # Loop through each file path
  for (file_path in file_paths) {
    # Read the RDS file
    df <- readRDS(file_path)
    # Get the column names (headers)
    headers <- colnames(df)
    # Store the headers in the list
    headers_list[[file_path]] <- headers
  }
  # Return the list of headers
  return(headers_list)
}


headers_list <- get_headers_from_rds_files(file_paths)

# Print the list of headers
print(headers_list)
# Find the ones that vary between files:

find_differing_columns <- function(file_paths) {
  # Initialize an empty list to store headers
  headers_list <- lapply(file_paths, function(file_path) colnames(readRDS(file_path)))
  
  # Get the unique columns across all files
  all_columns <- unique(unlist(headers_list))
  
  # Find the differing columns for each file
  differing_columns <- sapply(headers_list, function(headers) setdiff(all_columns, headers), simplify = FALSE)
  
  # Name the list elements with the corresponding file paths
  names(differing_columns) <- file_paths
  
  # Return the list of differing columns
  return(differing_columns)
}


differing_columns <- find_differing_columns(file_paths)

setdiff(headers_list)

# Need biome_string, analysis_experiment-type, analysis_pipeline-version,
# analysis_accession, analysis_instrument-platform, study_acc_type, sample_biosample, sample_accession, sample_environment-biome
# assembly_accession


# Define the function
read_and_bind_rds_files <- function(file_paths, selected_columns) {
  # Initialize an empty list to store dataframes
  df_list <- list()
  
  # Loop through each file path
  for (file_path in file_paths) {
    # Read the RDS file
    df <- readRDS(file_path)
    
    # Get the column names in the current dataframe
    current_columns <- colnames(df)
    
    # Find columns that are missing in the current dataframe
    missing_columns <- setdiff(selected_columns, current_columns)
    
    # Add the missing columns with NA values
    for (col in missing_columns) {
      df[[col]] <- NA
    }
    
    # Select the desired columns
    df <- df %>% select(all_of(selected_columns))
    
    # Append the dataframe to the list
    df_list <- append(df_list, list(df))
  }
  
  # Bind all dataframes together into a single dataframe
  combined_df <- bind_rows(df_list)
  
  # Return the combined dataframe
  return(combined_df)
}
  

# Specify the columns to select

selected_columns<-c("biome_string", "analysis_experiment-type", "analysis_pipeline-version",
  "analysis_accession", "analysis_instrument-platform", "study_acc_type", "sample_biosample",
  "sample_accession", "sample_environment-biome")

# Call the function
combined_df <- read_and_bind_rds_files(file_paths, selected_columns)


# Write out results file:

write.csv(combined_df, "PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_metadata_01072024.csv",
          row.names = FALSE, quote = FALSE)



