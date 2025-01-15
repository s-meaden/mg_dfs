
library(tidyverse)

# Script to select just viral sequences with sufficent length
# and 1+ hallmark genes

df<-read.table("PATH_TO_CONCATENATED_GENOMAD_SUMMARY_FILES/combined_viral_summaries.tsv",
               sep = "\t", header = T)

df<-df %>%
  filter( !seq_name == "seq_name") # remove leftover header from concatenating

# Filter out <1kb and lacking any hallmark genes.

df2<-df %>%
  mutate( length = as.numeric(length),
          n_hallmarks = as.numeric(n_hallmarks)) %>%
  filter( n_hallmarks > 0 & length >= 1000)

viral_seq_names<-df2$seq_name

write.table(viral_seq_names, "PATH_TO_CONCATENATED_GENOMAD_SUMMARY_FILES/filtered_viral_contig_names.txt",
            row.names = F, quote = F)
