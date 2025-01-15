
library(tidyverse)
library(MGnifyR)

# Locate URL for target assemblies (identified from metadata attributes)

# Set up client
mgclnt <- MgnifyClient(usecache = T, cache_dir = 'PATH_TO_CACHE/tmp/MGnify_cache/')

df<-read.csv("PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_list_sampled_01072024.csv")

# Download fasta and fastq data for these.

analyses_accessions<-df$analysis_accession

# Test on single
#dl_urls <- searchFile(mgclnt, "MGYA00708932", type = "analyses")

# Run for all
dl_urls <- searchFile(mgclnt, analyses_accessions, type = "analyses")

# Parse out URLs for FASTA files
dl_urls<-dl_urls %>%
  mutate( keep = ifelse( grepl("FASTA.fasta.gz", download_url), "yes", "no")) %>%
  filter( keep == "yes") %>%
  select( download_url)

# Write out so can download to server:

write.table(dl_urls, "PATH_TO_PROJECT/n100_URLs.txt",
            sep = "\t", row.names = F, quote = F)



