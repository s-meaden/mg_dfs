
library(tidyverse)

# Script to select samples for downstream analysis. Aiming for balanced dataset 
# with similar sample sizes for each environment

combined_df<-read.csv("PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_metadata_01072024.csv")

# Get summary stats:

combined_df %>%
  tally()

# 34,000 samples.

sorted_groups<-combined_df %>%
  group_by( biome_string) %>%
  tally() %>%
  arrange( -n)

unique(sorted_groups$biome_string)

tst<-combined_df %>%
  mutate( root_category = str_extract(biome_string, "^root:[^:]+"),
  sub_category = str_remove(biome_string, "^root:[^:]+:"),
  sub_sub_category = str_remove(str_remove(biome_string, "^root:[^:]+:[^:]+:"), "^[^:]+:") )

# Count number of samples per group
table(tst$sub_sub_category)


# Pick 100 of each group from above table() results

gut_vector<-c('Fecal',
'Intestine',
'Intestine:Fecal',
'Large intestine:Fecal',
'Digestive system',
'Gastrointestinal tract',
'Large intestine') # Many different classifications for gut samples. Collapse w/ vector

groups_w_over_100<-tst %>%
  filter( analysis_instrument.platform == "ILLUMINA") %>%
  mutate( sub_sub_category = ifelse( sub_sub_category %in% gut_vector, "gut", sub_sub_category)) %>%
  group_by( sub_sub_category) %>%
  tally() %>%
  filter( n > 99)

set.seed(187)
tst2<-tst %>%
  mutate( sub_sub_category = ifelse( sub_sub_category %in% gut_vector, "gut", sub_sub_category)) %>%
  filter( sub_sub_category %in% groups_w_over_100$sub_sub_category) %>%
  group_by( sub_sub_category) %>%
  sample_n( size = 100) %>%
  filter( !sub_sub_category == "Host-associated:Human") %>%
  filter( !sub_sub_category == "Hindgut:Rectum") %>%
  filter( !sub_sub_category == "Mixed")

unique(tst2$sub_sub_category)
dim(tst2)
table(tst2$sub_sub_category)

# Write out file with samples list

write.csv(tst2, "PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_list_sampled_01072024.csv")


