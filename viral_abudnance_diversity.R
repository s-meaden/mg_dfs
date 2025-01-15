
library(tidyverse)

# Collect viral abundances from each metagenome sample.

# Mapping stats are derived from 1 million reads mapped to the assembly, then extracting
# viral (or chromosomal but carrying a prophage) read recruitment / coverage data.

# Update to use coverage from original BAM files rather than masked ones (samtools view doing weird things and inflating read count)
df<-read.csv("PATH_TO_PROJECT/master_dataset_viral_contig_read_count_v2.csv")

# Total read recruitment and total viruses
# per sample with 1+ read mapped.

df2<-df %>%
  group_by( sample) %>%
  summarise( tot = sum(read_count), richness = length(unique(seqid)))

# Remove samples with < 1 million reads and link up with metadata.

meta<-read.csv("PATH_TO_PROJECT/mgnify_analyses_5.0_assembly_list_sampled_01072024.csv")

link<-read.csv("PATH_TO_PROJECT/n100_link_file.csv")

meta<-meta %>%
  left_join( link, by = "sample_accession") %>%
  select( sub_sub_category, assembly_accession) %>%
  unique()


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

# List of samples to retain
sample_counts_keep<-sample_counts %>%
  mutate( keep = ifelse(read_count >= 1000000, "YES", "NO")) %>%
  filter( keep == "YES")

sample_counts_keep<-unique(sample_counts_keep$file)

samples_to_keep<-link %>%
  filter( sample_accession %in% sample_counts_keep) %>%
  select( assembly_accession) 

samples_to_keep<-samples_to_keep$assembly_accession

# Link with metadata and tidy up
df3<-df2 %>%
  filter( sample %in% samples_to_keep) %>%
  rename( assembly_accession = sample) %>%
  left_join( meta, by = "assembly_accession") %>%
  mutate( env = case_when(sub_sub_category == "gut" ~ "Gut",
                          sub_sub_category == "Engineered:Wastewater" ~ "Water and sludge",
                          sub_sub_category == "Oceanic" ~ "Marine",
                          sub_sub_category == "Oral:Saliva" ~ "Oral",
                          sub_sub_category == "Oral:Subgingival plaque" ~ "Oral",
                          .default = as.character(sub_sub_category))) %>%
  distinct()

####### Plots of viral count and richness per environment

# Link up with defence abundances per sample.

def_counts<-read.csv("PATH_TO_PROJECT/n100_total_defence_counts_v2.csv",
                     header = T)

full_df<-def_counts %>%
  select( -env) %>%
  rename( def_count = tot) %>%
  left_join( df3, by = "assembly_accession") %>% 
  rename( viral_count = tot)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

a<-full_df %>%
  na.omit() %>%
  ggplot(., aes( def_count, viral_count))+
  geom_point( aes(color = env))+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth( method = "lm", linetype = "dashed", color = "black")+
  xlab("Defence Abundance")+
  ylab("Viral Abundance")+
  theme_classic()+
  theme( text = element_text(size = 15))+
  labs(color = "Environment")+
  scale_color_brewer( palette = "Paired")+
  guides(color = "none")

a

# Include N50 data
n50<-read.csv("PATH_TO_PROJECT/N50_assembly_values.csv", header = F)

full_df2<-n50 %>%
  rename( assembly_accession = V1, n50 = V2) %>%
  mutate( assembly_accession = gsub("_FASTA", "", assembly_accession)) %>%
  left_join( ., full_df, by = "assembly_accession") %>%
  na.omit()

# Run GLM
hist(log10(full_df2$def_count))
# Include n50 in the model
m1<-glm(log10(def_count) ~ n50 + log10(viral_count), data = full_df2)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
m2<-update(m1,~.-log10(viral_count), data = full_df2)
anova(m1, m2, test = "F")


# Would expect similar trend with diversity. Although abundance and diversity are highly correlated.
full_df %>%
  na.omit() %>%
  ggplot(., aes( richness, def_count))+
  geom_point( aes(color = env))+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth( method = "lm", linetype = "dashed", color = "black")+
  xlab("Viral Richness")+
  ylab("Defence Abundance")+
  theme_classic()+
  theme( text = element_text(size = 15))+
  labs(color = "Environment")+
  scale_color_brewer( palette = "Paired")



# Abundance correlation on just caudovirales abundances.

# Read in Caudovirales abundances:

caudo<-read.csv("PATH_TO_PROJECT/n100_caudo_abundance.csv")

b<-full_df %>%
  rename( sample = assembly_accession) %>%
  left_join( caudo, by = "sample") %>%
  na.omit() %>%
  ggplot(., aes( def_count, tot))+
  geom_point( aes(color = env))+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth( method = "lm", linetype = "dashed", color = "black")+
  xlab("Defence Abundance")+
  ylab("Caudovirales Abundance")+
  theme_classic()+
  theme( text = element_text(size = 15))+
  labs(color = "Environment")+
  scale_color_brewer( palette = "Paired")

b

# Plot for paper:
cowplot::plot_grid(a, b, ncol = 2, align = "w")

# Update with egg
egg::ggarrange(a, b, ncol = 2, labels = c('A', 'B'))

grob_a <- ggplotGrob(a)
grob_b <- ggplotGrob(b)

legend_column <- 5
legend_width <- grob_b$widths[legend_column]
grob_combined <- gtable:::cbind_gtable(grob_a, grob_b) # use rbind_gtable if want single column
grid::grid.draw(grob_combined)

# Add panel labels manually.
grid::grid.text("A", x = 0.07, y = 0.975, gp = grid::gpar(fontsize = 16))
grid::grid.text("B", x = 0.49, y = 0.975, gp = grid::gpar(fontsize = 16))


# And GLM:

caudo_df<-full_df %>%
  rename( sample = assembly_accession) %>%
  left_join( caudo, by = "sample") %>%
  na.omit() %>%
  rename( assembly_accession = sample)

# Add N50 values
caudo_df2<-n50 %>%
  rename( assembly_accession = V1, n50 = V2) %>%
  mutate( assembly_accession = gsub("_FASTA", "", assembly_accession)) %>%
  left_join( ., caudo_df, by = "assembly_accession") %>%
  na.omit()

# Include n50 in the model
m1<-glm(log10(def_count) ~ n50 + log10(viral_count), data = caudo_df2)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
m2<-update(m1,~.-log10(viral_count), data = caudo_df2)
anova(m1, m2, test = "F")

