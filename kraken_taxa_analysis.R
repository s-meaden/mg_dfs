
# Load required packages. If you need to install them this can be done
# as follows:
# install.packages("tidyverse", repos = "https://www.stats.bris.ac.uk/R/") # repos part may not be necessary and can be changed to a different mirror URL.
library(tidyverse)
library(vegan)
library(stringi)
library(ecodist)

# Script to take output of an aggregated Kraken report
# and run diversity analyses (alpha, beta and clustering plot)

# Read in Kraken report (not sure if this works for a Bracken report also)
# Update the path to the correct location in your filesystem

# Format causes some issues. Run this command in terminal
# sed 's/[[:space:]]/_/g' kraken_report_master_file.txt > kraken_report_master_file_copy.txt
df<-data.table::fread("PATH_TO_KRAKEN_MASTER_OUTPUT/kraken_report_master_file_copy.txt", sep = "\t", quote = "\"",
                      header = F)

# Extract relevant info into columns using various regexes
# bit messy
df2<-df %>%
  mutate(perc = stringi::stri_extract_first_regex(V1, "[0-9]+\\.[0-9]+")) %>%
  mutate( count = gsub( "[0-9]+\\.[0-9]+_", "", stri_extract_first_regex(V1, "[0-9]+\\.[0-9]+_[0-9]+"))) %>%
  mutate(sample = stringi::stri_extract_last_regex(V1, "[A-Z]+[0-9]+")) %>%
  mutate(level = stringi::stri_extract_first_regex(V1, "[A-Z]{1}")) %>%
  mutate(taxa = gsub("_[A-Z]+[0-9]+", "", stringi::stri_extract_last_regex(V1, "([0-9]+_[A-Za-z]+).*"))) %>%
  mutate( taxa2 = sub(".*?([A-Z][a-zA-Z].*)", "\\1", taxa)) %>%
  mutate( taxa3 = case_when( grepl("unclassified", taxa2) ~ "unclassified",
                     grepl("root", taxa2) ~ "root",
                    .default = taxa2 )) %>%
  select( perc, count, sample, level, taxa3) %>%
  rename( taxa = taxa3)

head(df2)

# And remove samples that had < 1 million reads
# (consistent with other analyses)

# Repeat with  just samples that had exactly 1 million reads:

r1<-read.table("PATH_TO_PROJECT/fastq_R1_wc_stats.txt", header = F, sep = " ")
r2<-read.table("PATH_TO_PROJECT/fastq_R2_wc_stats.txt", header = F, sep = " ")

colnames(r1)<-c('R1_count', 'file')
colnames(r2)<-c('R2_count', 'file')

r1<-r1[-1,] # remove first line as was empty wc -l search result.
r2<-r2[-1,]

r1<-r1 %>% mutate(file = gsub("_R1.fq", "", file))
r2<-r2 %>% mutate(file = gsub("_R2.fq", "", file))

sample_counts<-r1 %>%
  left_join( r2, by = "file") %>%
  mutate( read_count = (R1_count + R2_count) / 4) %>%
  select( file, read_count) %>%
  mutate( file = gsub("_subsampled", "", file)) %>%
  mutate( file = gsub("subsampled_fastqs_combined/" , "", file))

table( ifelse(sample_counts$read_count >= 1000000, "YES", "NO"))

sample_counts_keep<-sample_counts %>%
  mutate( keep = ifelse(read_count >= 1000000, "YES", "NO")) %>%
  filter( keep == "YES")

sample_counts_keep<-unique(sample_counts_keep$file)

samples_to_keep<-link %>%
  filter( sample_accession %in% sample_counts_keep) %>%
  select( assembly_accession) 

samples_to_keep<-samples_to_keep$assembly_accession


# Step 1.
# Assess how may samples have extreme amounts of unclassified reads.

# Link up with metadata:
meta<-read.csv("PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_list_sampled_01072024.csv")

meta<-meta %>%
  select( -X, -analysis_accession) %>%     # Remove leftover rownames column and analysis accession (multiple analyses can be done for same sample)
  distinct()

link<-read.csv("PATH_TO_PROJECT/n100_link_file.csv")

table(meta$sub_sub_category)

meta<-meta %>%
  left_join( link, by = "sample_accession") %>%
  filter( assembly_accession %in% samples_to_keep) %>%
  select( sub_sub_category, sample_accession) %>%
  distinct()

# Make more meaningful categories- collapse redundant groups

df3<-df2 %>%
  rename( sample_accession = sample) %>%
  inner_join( meta, by = "sample_accession") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category)))

# Plot unclassified reads
df3 %>%
  filter( level == "U") %>%
  mutate( perc = as.numeric(perc)) %>%
  ggplot(., aes(perc))+
  geom_histogram( color = "black", fill = "white")+
  facet_wrap( ~ env)+
  xlab("Percentage of\nUnclassified Reads")+
  ylab("Count")+
  theme_classic()+
  theme( text = element_text( size = 15))

# So rhizosphere and sludge quite a lot of unclassified reads. Makes sense if complex, understudied envs.

#### What about eukaryotic (human) levels?

# Check samples sum to 100%. Remove any that don't...

df3 %>%
  filter( level == "U" | level == "R") %>%
  mutate( perc = as.numeric(perc)) %>%
  group_by( sample_accession) %>%
  summarise( total = sum(perc)) %>%
  arrange( total) %>%
  head()

# Maximum percentage total is 100
# but minimum can be < 100. Inspect

df3 %>%
  filter( sample_accession == "SRS11685460") %>%
  arrange( -as.numeric(perc)) %>%
  head(n = 20)

# This sample has no unclassified reads. Check original data:

df2 %>%
  filter( sample == "SRS11685460")

# Remove to be safe:

df3<-df3 %>%
  filter( !sample_accession == "SRS11685460")


# Now plot eukaryotic levels:

df3 %>%
  filter( taxa == "Eukaryota") %>%
  group_by( sample_accession) %>%
  mutate( perc = as.numeric(perc)) %>%
  mutate( med_euk = median(perc)) %>%
  ggplot(., aes( reorder(env, -med_euk), perc))+
  geom_jitter( width = 0.1)+
  geom_boxplot( alpha = 0.7)+
  xlab( "Environment")+
  ylab("Percentage of\nEukaryotic Reads")+
  theme_classic()+
  theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 15))

# Identify samples with >50% eukaryotic.

euk_50<-df3 %>%
  filter( taxa == "Eukaryota") %>%
  mutate( perc = as.numeric(perc)) %>%
  filter( perc > 50) %>%
  select( sample_accession)

#write.table(euk_50, "PATH_TO_PROJECT/samples_w_50_perc_eukaryotic.txt",
#            row.names = F, quote = F)

# Run models with these samples excluded.

## Cluster based on Genus level and look for environmental signature.

gen_df<-df3 %>%
  filter( !sample_accession %in% euk_50$sample_accession) %>%
  filter( level == "G")

gen_df_env<-gen_df %>%
  select(sample_accession, env)

# Make species abundance table
spec_mat<-gen_df %>%
  select(sample_accession, perc, taxa) %>%
  mutate( perc = as.numeric(perc)) %>%
  mutate(taxa = str_trim(taxa)) %>% # remove whitespaces from these fields.
  pivot_wider( names_from = taxa, values_from = perc)

# Sort out matrix properties.
spec_mat<-as.matrix(spec_mat)   # force R to treat this as a matrix object (not a dataframe)
row.names(spec_mat)<-spec_mat[,1] # Convert the first column to the 'row name'. This is bacterial species.
spec_mat<-spec_mat[,-1] # remove the first column as this info is now stored as the row name.
class(spec_mat) <- "numeric" # Enforce treating the data as numeric 
spec_mat[is.na(spec_mat)] <- 0 # Convert NAs (where a species was not detected) to a zero


# PCoA analysis

bray_curtis_dist <- vegdist(spec_mat, method = "bray") # This vegdist() function from the 'vegan' package 

bray_curtis_pcoa <- pco(bray_curtis_dist) # Store results of pcoa analysis

bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1],  # Extract the first axis of the analysis
                                  pcoa2 = bray_curtis_pcoa$vectors[,2],  # Extract the second. You could inspect the 3rd, 4th etc
                                  sample = row.names(bray_curtis_pcoa$vectors))

# Create a plot. This follows the tidyverse syntaxt using %>% to 'pipe' or chain commands together.
# It passes the data through some filters and into ggplot. 

bray_curtis_pcoa_df %>%
  rename( sample_accession = sample) %>%
  left_join( gen_df_env, by = "sample_accession") %>%
  ggplot(., aes(x=pcoa1, y=pcoa2)) +
  geom_point( aes( color = env), size = 4) +
  geom_point(shape = 1,size = 4,colour = "black")+
  labs(x = "PC1",
       y = "PC2") +
  theme_classic()+
  theme(text = element_text(size = 15))+
  scale_color_brewer( palette = "Set3")+
  stat_ellipse(linetype = 2, linewidth = 1, type = "t", aes(color = env))+
  labs(color = "Environment")


## Extract diversity metrics.

# Need to be aware of eukaryotic level and unclassified level

df3 %>%
  filter( !sample_accession %in% euk_50$sample_accession) %>%
  filter( level == "G") %>%
  select( sample_accession, taxa) %>%
  distinct() %>%
  group_by( sample_accession) %>%
  tally() %>%
  arrange( -n) %>%
  left_join( gen_df_env, by = "sample_accession") %>%
  distinct() %>%
  ungroup() %>%
  group_by( env) %>%
  mutate( median = median(n)) %>%
  ggplot(., aes(reorder(env, -median), n))+
  geom_jitter( width = 0.1, size = 1)+
  geom_boxplot( alpha = 0.7, outliers = F)+
  theme_classic()+
  xlab("Environment")+
  ylab("Unique Genera\nCount")+
  theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust=1))+
  theme( text = element_text(size = 15))

# Which extra genera are included here?

tmp<-df3 %>%
  filter( !sample_accession %in% euk_50$sample_accession) %>%
  filter( level == "G")

unique(tmp$taxa)

# Collect Shannin diversity index scores.

diversity(spec_mat, index = "shannon") %>%
  reshape2::melt() %>%
  mutate( index = "shannon") %>%
  mutate( sample_accession = row.names(.)) %>%
  filter( !sample_accession %in% euk_50$sample_accession) %>%
  left_join( gen_df_env, by = "sample_accession") %>%
  distinct() %>%
  group_by( env) %>%
  mutate( median = median(value)) %>%
  ggplot(., aes(reorder(env, -median), value))+
  geom_jitter( width = 0.1, size = 1)+
  geom_boxplot( alpha = 0.7, outliers = F)+
  theme_classic()+
  xlab("Environment")+
  ylab("Shannon Diversity")+
  theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust=1))+
  theme( text = element_text(size = 15))

