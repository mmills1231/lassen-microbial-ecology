# Rarefaction for 16S dataset ####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S")

# Load filtered 16S phyloseq object.
physeq_16S_filtered <- readRDS("physeq_16S_filtered.rds")

# Rarefy 16S data for Shannon diversity analysis.
rare_16S <- rarefy_even_depth(physeq_16S_filtered, rngseed = 123, sample.size = 4400)

# Save phyloseq objects.
saveRDS(physeq_16S_filtered, "nonrare_16S.rds") # renaming for convenience
saveRDS(rare_16S, "rare_16S.rds")