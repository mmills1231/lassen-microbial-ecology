# DADA2 processing of Lassen 16S rRNA amplicon sequencing data ####

# Load required packages.
library(dada2)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S")

# Define the folder containing demultiplexed FASTQ files.
path <- "~/OneDrive - Texas A&M University/Lassen2017/Lassen16S"

# List files in the input folder.
list.files(path)

# Import paired-end FASTQ filenames.
fnF <- sort(list.files(path, pattern="_R1_001.fastq", full.names=TRUE))
fnR <- sort(list.files(path, pattern="_R2_001.fastq", full.names=TRUE))

# Extract sample names from FASTQ filenames.
sample_names <- sapply(strsplit(basename(fnF), "_R"), `[`, 1)

# Plot forward and reverse read quality profiles.
plotQualityProfile(fnF[1:20])
plotQualityProfile(fnR[1:20])

# Create filtered read output filenames.
filtF <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtR <- file.path(path, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtF) <- sample_names
names(filtR) <- sample_names

# Filter and trim paired-end reads.
filter_out <- filterAndTrim(fnF, filtF, fnR, filtR, truncLen=c(200,190), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE)

# Inspect the filtering summary.
head(filter_out)

# Learn error rates.
errF <- learnErrors(filtF, multithread=FALSE)
errR <- learnErrors(filtR, multithread=FALSE)

# Plot estimated error rates.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate identical reads.
derepF <- derepFastq(filtF, verbose=TRUE)
derepR <- derepFastq(filtR, verbose=TRUE)

# Name dereplicated objects by sample.
names(derepF) <- sample_names
names(derepR) <- sample_names

# Infer ASVs.
dadaF <- dada(derepF, err=errF, multithread=FALSE)
dadaR <- dada(derepR, err=errR, multithread=FALSE)

# Inspect ASV inference for the first sample.
dadaF[[1]]

# Merge paired-end reads.
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE)

# Inspect merged reads for the first sample.
head(mergers[[1]])

# Create the ASV table.
seqtab_16S <- makeSequenceTable(mergers)

# Check ASV table dimensions.
dim(seqtab_16S)

# Inspect ASV sequence length distribution.
table(nchar(getSequences(seqtab_16S)))

# Keep ASVs within the expected 16S amplicon length range.
seqtab_16S_len <- seqtab_16S[, nchar(colnames(seqtab_16S)) %in% seq(370,378)]

# Remove chimeric ASVs.
seqtab_16S_nochim <- removeBimeraDenovo(seqtab_16S_len, method="consensus", multithread=FALSE, verbose=TRUE)

# Check ASV table dimensions after chimera removal.
dim(seqtab_16S_nochim)

# Calculate the proportion of reads retained after chimera removal.
sum(seqtab_16S_nochim)/sum(seqtab_16S_len)

# Track read retention.
getN <- function(x) sum(getUniques(x))
track_16S <- cbind(filter_out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab_16S_nochim))
colnames(track_16S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_16S) <- sample_names

# Inspect read retention table.
head(track_16S)

# Save read retention table.
write.csv(track_16S, "Lassen_16S_DADA2_read_tracking.csv")

# Save final 16S ASV count table.
write.csv(seqtab_16S_nochim, "Lassen_16S_ASV_table.csv")

# Assign 16S taxonomy using SILVA.
tax_16S <- assignTaxonomy(seqtab_16S_nochim, "data/references/silva_nr_v132_train_set.fa.gz", multithread=FALSE)

# Inspect taxonomy table.
tax_16S_print <- tax_16S
rownames(tax_16S_print) <- NULL
head(tax_16S_print)

# Save 16S taxonomy table.
write.csv(tax_16S, "Lassen_16S_ASV_taxonomy.csv")

