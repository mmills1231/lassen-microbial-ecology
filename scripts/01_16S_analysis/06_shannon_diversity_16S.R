# Create 16S Shannon diversity tables #####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S")

# Load filtered 16S phyloseq object.
nonrare_16S <- readRDS("nonrare_16S.rds")
rare_16S <- readRDS("rare_16S.rds")

# Calculate rarefied Shannon diversity and add metadata.
rare_rich_16S <- estimate_richness(rare_16S, measures = "Shannon") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_data(rare_16S) %>% data.frame() %>% rownames_to_column(var = "SampleID"), by = "SampleID")

# Calculate non-rarefied Shannon diversity and add metadata.
nonrare_rich_16S <- estimate_richness(nonrare_16S, measures = "Shannon") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_data(nonrare_16S) %>% data.frame() %>% rownames_to_column(var = "SampleID"), by = "SampleID")

# Save Shannon diversity tables.
write.csv(rare_rich_16S, "Shannon_metadata_16S_rarefied.csv", row.names = FALSE)
write.csv(nonrare_rich_16S, "Shannon_metadata_16S_nonrarefied.csv", row.names = FALSE)