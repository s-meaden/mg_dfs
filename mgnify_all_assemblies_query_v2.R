
library(MGnifyR)
library(dplyr)
library(readr)

# Query MGnify database to get all possible metagenomic assembly accessions
# Must be assembly (not amplicon) and processed in pipeline 5.0

# Run locally as problems with installing MGnifyR on server

#Set up the MGnify client instance
mgclnt <- MgnifyClient(usecache = T, cache_dir = './tmp/MGnify_cache')

# Run full query

doQuery(mgclnt, type = "analyses", accession = NULL,
        asDataFrame = T, max.hits = NULL, usecache = F) %>%
  filter( `pipeline-version` == "5.0" & `experiment-type` == "assembly", `is-private` == FALSE) %>%
  write_csv( "./mgnify_assemblies_query_14062024.csv")

# Took 3.5 days to run on desktop.


