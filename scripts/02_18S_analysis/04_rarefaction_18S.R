# Rarefaction for 18S dataset ####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Load filtered 16S phyloseq object.
physeq_18S_filtered <- readRDS("physeq_18S_filtered.rds")

# Rarefy 16S data for Shannon diversity analysis.
rare_18S <- rarefy_even_depth(physeq_18S_filtered, rngseed = 123, sample.size = 3088)

# Save phyloseq objects.
saveRDS(physeq_18S_filtered, "nonrare_18S.rds") # renaming for convenience
saveRDS(rare_16S, "rare_18S.rds")