# Create 18S Shannon diversity tables #####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Load filtered 18S phyloseq objects.
nonrare_18S <- readRDS("nonrare_18S.rds")
rare_18S <- readRDS("rare_18S.rds")

# Calculate rarefied Shannon diversity and add metadata.
rare_rich_18S <- estimate_richness(rare_18S, measures = "Shannon") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_data(rare_18S) %>% data.frame() %>% rownames_to_column(var = "SampleID"), by = "SampleID")

# Calculate non-rarefied Shannon diversity and add metadata.
nonrare_rich_18S <- estimate_richness(nonrare_18S, measures = "Shannon") %>%
  rownames_to_column(var = "SampleID") %>%
  left_join(sample_data(nonrare_18S) %>% data.frame() %>% rownames_to_column(var = "SampleID"), by = "SampleID")

# Save Shannon diversity tables.
write.csv(rare_rich_18S, "Shannon_metadata_18S_rarefied.csv", row.names = FALSE)
write.csv(nonrare_rich_18S, "Shannon_metadata_18S_nonrarefied.csv", row.names = FALSE)