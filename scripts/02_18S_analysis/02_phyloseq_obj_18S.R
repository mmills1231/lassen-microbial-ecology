# Create 18S phyloseq object for Lassen microbial ecology study ####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Read the 18S ASV count table.
asv_18S <- read.csv("Lassen_18S_ASV_table.csv", header=TRUE, row.names=1, check.names=FALSE)

# Read the 18S taxonomy table.
tax_18S <- read.csv("Lassen_18S_ASV_taxonomy.csv", header=TRUE, row.names=1, check.names=FALSE)

# Read the sample metadata.
meta_18S <- read.csv("lassen_metadata.csv", header=TRUE, row.names=1, check.names=TRUE, stringsAsFactors=FALSE)

# Keep samples from the locations and sample types.
meta_18S_filtered <- meta_18S %>%
  rownames_to_column(var="SampleID") %>%
  filter(Location %in% c("Bumpass Hell", "Devils Kitchen", "Pilot Pinnacle", "Sulphur Works Acid", "Sulphur Works Neutral"),
         Sampletype %in% c("spring fluid", "stream", "sediment", "crusted sediment", "biofilm", "precipitate")) %>%
  column_to_rownames(var="SampleID")

# Create phyloseq components.
ASV_18S <- otu_table(as.matrix(asv_18S), taxa_are_rows=TRUE)
TAX_18S <- tax_table(as.matrix(tax_18S))
META_18S <- sample_data(meta_18S_filtered)

# Create the unfiltered 18S phyloseq object.
physeq_18S <- phyloseq(ASV_18S, TAX_18S, META_18S)

# Keep eukaryotic ASVs, remove bacterial, metazoan, and plant-associated assignments, and remove samples with fewer than 1,000 reads.
physeq_18S_filtered <- physeq_18S %>%
  subset_taxa(Kingdom == "Eukaryota" & Kingdom != "Bacteria" & Division != "Metazoa" & Division != "Streptophyta") %>%
  prune_samples(sample_sums(.) >= 1000, .)

# Check the final filtered phyloseq object.
physeq_18S_filtered

# Save the filtered phyloseq objects.
saveRDS(physeq_18S_filtered, "physeq_18S_filtered.rds")
