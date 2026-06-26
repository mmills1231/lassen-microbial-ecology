# Create 16S phyloseq object for Lassen microbial ecology study ####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S")

# Read the 16S ASV count table.
asv_16S <- read.csv("Lassen_16S_ASV_table.csv", header=TRUE, row.names=1, check.names=FALSE)

# Read the 16S taxonomy table.
tax_16S <- read.csv("Lassen_16S_ASV_taxonomy.csv", header=TRUE, row.names=1, check.names=FALSE)

# Read the sample metadata.
meta_16S <- read.csv("lassen_metadata.csv", header=TRUE, row.names=1, check.names=TRUE, stringsAsFactors=FALSE)

# Keep samples from the locations and sample types.
meta_16S_filtered <- meta_16S %>%
  rownames_to_column(var="SampleID") %>%
  filter(Location %in% c("Bumpass Hell", "Devils Kitchen", "Pilot Pinnacle", "Sulphur Works Acid", "Sulphur Works Neutral"),
         Sampletype %in% c("spring fluid", "stream", "sediment", "crusted sediment", "biofilm", "precipitate")) %>%
  column_to_rownames(var="SampleID")

# Create phyloseq components.
ASV_16S <- otu_table(as.matrix(asv_16S), taxa_are_rows=TRUE)
TAX_16S <- tax_table(as.matrix(tax_16S))
META_16S <- sample_data(meta_16S_filtered)

# Create the unfiltered 16S phyloseq object.
physeq_16S <- phyloseq(ASV_16S, TAX_16S, META_16S)

# Keep bacterial and archaeal ASVs, remove chloroplasts and mitochondria, and remove samples with fewer than 1,000 reads.
physeq_16S_filtered <- physeq_16S %>%
  subset_taxa((Kingdom == "Bacteria" | Kingdom == "Archaea") & Order != "Chloroplast" & Family != "Mitochondria") %>%
  prune_samples(sample_sums(.) >= 1000, .)

# Check the final filtered phyloseq object.
physeq_16S_filtered

# Save the filtered phyloseq objects.
saveRDS(physeq_16S_filtered, "physeq_16S_filtered.rds")
