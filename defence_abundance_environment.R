
library(tidyverse)
library(vegan)

# Plot environmental distribution of defences.

# First compare read recruitment vs coverage.

# Updated version to fix samtools view read inflation problem

counts<-read.csv("PATH_TO_PROJECT/master_dataset_defence_read_count2.csv")

cov<-read.csv("PATH_TO_PROJECT/master_dataset_defence_abundance_v2.csv")

# Combine datasets and plot correlation:

tmp_cov<-cov %>%
  mutate( idx = paste(sample, seqid, sep = "_")) %>%
  select(coverage, idx) %>%
  unique()


counts %>%
  mutate( idx = paste(sample, seqid, sep = "_")) %>%
  select(read_count, idx) %>%
  unique() %>%
  left_join( tmp_cov, by = "idx") %>%
  sample_n(10000) %>%
  ggplot(., aes(coverage, read_count))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

##### As expect coverage and read recruitment strongly correlated (non-linear.)

counts %>%
  group_by( sample, system) %>%
  summarise( tot = sum(read_count)) %>%
  arrange( -tot) %>%
  head()

# Work with the read recruitment as count data better for ecological models.

# Filter out spurious defences and candidates:

grep("DMS", unique(counts$system))
grep("dXTP", unique(counts$system))
grep("PDC", unique(counts$system))

counts<-counts %>%
  filter( !system == "DMS_other") %>%
  filter( !system == "dXTPase") %>%
  filter( !grepl("PDC",system)) %>%
  filter( !grepl("HEC",system)) %>%
  filter( !system == "VSPR")


# And remove duplicates in case multi-gene systems on same contig.

counts %>%
  group_by(seqid, system) %>%
  tally() %>%
  arrange(-n) %>%
  head()

# Note, no system occurs twice on the same contig because we used the unique command when linking (defence_count_link_v2.R)
# Intentional choice not to count duplicates of same systems on same contig
# But up to 10 systems on one contig. Analysis here on is looking at the abundance of the contig carrying the defence,
# rather than the defence itself. See manuscript for details.

# Remove NAs from count data:
min(counts$read_count)

counts<-counts %>%
  drop_na()

# Check how many samples made it through the pipeline
length(unique(counts$sample))
# n = 1,145

# Most abundant overall:
min(counts$read_count)
max(counts$read_count)

counts %>%
  group_by( system) %>%
  summarise( tot = sum(read_count)) %>%
  arrange( -tot) %>%
  head(n = 20) %>%
  ggplot(., aes( reorder(system, -tot), tot))+
  geom_bar( stat = "identity")+
  theme_classic()+
  theme( axis.text.x = element_text(angle = 90))+
  theme( text = element_text(size = 15))+
  xlab("Defence System")+
  ylab("Total Read Recruitment")

# Most widespread
  
counts %>%
  mutate( read_count = 1) %>%
  select(system, read_count, sample) %>%
  unique() %>%
  group_by(system) %>%
  tally() %>%
  arrange( -n) %>%
  head(n = 20) %>%
  ggplot(., aes( reorder(system, -n), n))+
  geom_bar( stat = "identity")+
  theme_classic()+
  theme( axis.text.x = element_text(angle = 90))+
  theme( text = element_text(size = 15))+
  xlab("Defence System")+
  ylab("Prevalence in Samples")


# Heatmap of defence systems.

def_mat<-counts %>%
  group_by(system, sample) %>%
  summarise( count = sum(read_count)) %>%
  select(sample, system, count) %>%
  pivot_wider(id_cols = sample, names_from = system, values_from = count, values_fill = 0) %>%
  as.data.frame()

row.names(def_mat)<-def_mat$sample
def_mat<-def_mat %>% select( -sample)

# Drop empty rows
def_mat<-def_mat[rowSums(def_mat[, -1] > 0) != 0, ]

# Run ordination:
# ordination by NMDS
NMDS <- metaMDS(def_mat, distance = "bray", k = 3)
NMDS

# Stress ~1.7 Not great but < 2. 
stressplot(NMDS)

# Plot ordination (just to see if clustering exists)
plot(NMDS, "sites")
orditorp(NMDS, "sites")

# Extract ordination info to link with metadata:
ord_data<-plot(NMDS, "sites")
ord_df<-as.data.frame(ord_data$sites)

# Link up with metadata:
meta<-read.csv("~PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_list_sampled_01072024.csv")

link<-read.csv("~PATH_TO_PROJECT/n100_link_file.csv")

meta<-meta %>%
  left_join( link, by = "sample_accession") %>%
  select( sub_sub_category, assembly_accession) %>%
  unique()

# Tidy up metadata

ord_df %>%
  mutate( assembly_accession = row.names(ord_df)) %>%
  left_join( meta, by = "assembly_accession") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category))) %>%
  ggplot(., aes(x = NMDS1, y = NMDS2, colour = env)) +
  geom_point() +                                              
  coord_fixed()+                                              
  theme_classic()+
  stat_ellipse(linetype = 2, linewidth = 1, type = "t")+
  #ggtitle("Country")+
  #scale_colour_brewer(type = "qual", palette = 3)+
  theme(text = element_text(size = 15))+
  labs(color = "Environment")+
  scale_color_brewer( palette = "Paired")

#### Heatmap of defences clustered by environment

df<-counts %>%
  group_by(system, sample) %>%
  summarise( count = sum(read_count)) %>%
  rename( assembly_accession = sample) %>%
  left_join( meta, by = "assembly_accession")

# Check no duplicates:
df %>%
  ungroup() %>%
  janitor::get_dupes() %>%
  head()

# Prevalence of contigs carrying multiple defences?
counts %>%
  group_by(seqid) %>%
  tally() %>%
  arrange( -n) %>% 
  head()

# Track one of these through to final plot:
counts %>%
  filter( seqid == "ERZ3455859.2-NODE-2-length-895680-cov-12.386554") %>%
  summarise( count = sum(read_count)) %>%
  head()

# Note we are counting different systems intentionally.
# How common is this scenario where multiple systems present on one contig:

counts %>%
  group_by(seqid) %>%
  tally() %>%
  arrange( -n) %>% 
  ggplot(., aes(n))+
  geom_histogram()+
  scale_y_log10()

# Overwhelmingly 1 system count per contig.
# Proportions of contigs carrying 1,2,3 defences etc

counts %>%
  group_by(seqid, sample) %>%
  tally() %>%
  arrange( -n) %>%
  ungroup() %>%
  group_by( n) %>%
  tally() %>%
  mutate( proportion = nn / sum(nn))


# Permanova for defence composition:

# Generate distance matrix for the samples:
def_dist<-as.matrix(vegdist(x = def_mat, method = "bray"))
str(def_dist)
tmp<-as.data.frame(def_dist)

# Combine with metadata.
# This ensures row orders match up for following PERMANOVA analyses
tmp<-tmp %>%
  mutate( assembly_accession = row.names(.)) %>%
  left_join( meta, by = "assembly_accession") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category)))


# Split data into species matric and metadata
adonis_df<-tmp[,1:1145]
adonis_meta<-tmp[, 1146:1148]
# Check datasets match
all.equal(rownames(adonis_df), rownames(adonis_meta))

# Run PERMANOVA
m1<-adonis2(adonis_df ~ env, data = adonis_meta, permutations = 999)
m1

##### Plot heatmap. Needs work
# From counts data.

tmp<-df %>%
  left_join( prev_df, by = "system") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category))) %>%
  mutate( env = gsub("Host-associated:", "", env)) %>%
  mutate( env = gsub("Intertidal zone:", "", env)) %>%
  mutate( env = gsub("Urethra:", "", env)) %>%
  mutate( env = gsub("Water and sludge", "Wastewater", env))

ggplot(tmp, aes(assembly_accession, reorder(system, n)))+
  geom_tile( aes(fill = log10(count)))+
  facet_wrap( ~env, scale = "free_x", ncol = 16)+
  theme(axis.text.x=element_blank(), text = element_text(size = 15),
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text( size = 10))+
  xlab("Metagenomic Sample")+
  ylab("Defence System")+
  scale_fill_viridis_c( option = "mako")

##### Variation in defence system abundance across environments:

options(scipen=999)

df %>%
  group_by( assembly_accession, sub_sub_category) %>%
  summarise( tot = sum(count)) %>%
  ggplot(., aes(reorder(sub_sub_category, -tot), tot))+
  geom_jitter( width = 0.1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()

# Run with  just samples that had exactly 1 million reads:

r1<-read.table("PATH_TO_PROJECT/fastq_R1_wc_stats.txt", header = F, sep = " ")
r2<-read.table("PATH_TO_PROJECT//fastq_R2_wc_stats.txt", header = F, sep = " ")

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

write.csv(samples_to_keep, "PATH_TO_PROJECT/samples_to_keep_2024.csv",
          row.names = F, quote = F)

######### Figure with just samples that retained 1 million reads per sample:

df %>%
  filter( assembly_accession %in% samples_to_keep) %>%
  group_by( assembly_accession, sub_sub_category) %>%
  summarise( tot = sum(count)) %>%
  ggplot(., aes(reorder(sub_sub_category, -tot), tot))+
  geom_jitter( width = 0.1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()





counts %>%
  rename( assembly_accession = sample) %>%
  filter( assembly_accession %in% samples_to_keep) %>%
  select( -system) %>%
  unique() %>%
  group_by(assembly_accession) %>%
  summarise( tot = sum(read_count))


tmp_df<-counts %>%
  rename( assembly_accession = sample) %>%
  filter( assembly_accession %in% samples_to_keep) %>%
  select( -system) %>%
  unique() %>%
  group_by(assembly_accession) %>%
  summarise( tot = sum(read_count)) %>%
  left_join( meta, by = "assembly_accession") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category))) %>%
  group_by(env) %>%
  mutate( median = median(tot)) %>%
  ungroup()

# Plot results
ggplot(tmp_df, aes(reorder(env, -median), tot))+
  geom_jitter( width = 0.1)+
  geom_boxplot( alpha = 0.6)+
  scale_y_log10()+
  theme_classic()+
  xlab("Environmental Category")+
  ylab("Defence Abundance\n(Read Count)")+
  theme( axis.text.x = element_text(angle = 25, hjust = 1))+
  theme( text = element_text(size = 15))

# Include N50 (assembly fragmentation) info for models

n50<-read.csv("PATH_TO_PROJECT/N50_assembly_values.csv", header = F)

n50 %>%
  rename( assembly_accession = V1, n50 = V2) %>%
  mutate( assembly_accession = gsub("_FASTA", "", assembly_accession)) %>%
  left_join( ., tmp_df, by = "assembly_accession") %>%
  ggplot(., aes(tot, n50))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth( method = "lm", formula = y~poly(x, 2), linetype = "dashed", color = "black")+
  theme_classic()+
  xlab("Defence Abundance")+
  ylab("N50")+
  theme( text = element_text(size = 15))

tmp_df_3<-n50 %>%
  rename( assembly_accession = V1, n50 = V2) %>%
  mutate( assembly_accession = gsub("_FASTA", "", assembly_accession)) %>%
  left_join( ., tmp_df, by = "assembly_accession")

# Run GLM

hist(log10(tmp_df$tot)) 

m1<-glm(tot ~ n50 + env, data = tmp_df_3, family = "quasipoisson")
coefplot::coefplot(m1)
summary(m1)
m2<-update(m1,~.-env)
anova(m1, m2, test = "F")

# Pseudo R2
with(summary(m1), 1 - deviance/null.deviance)
summary(m1)

