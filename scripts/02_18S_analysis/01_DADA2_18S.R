# DADA2 processing of Lassen 18S rRNA amplicon sequencing data####

# Load required packages.
library(dada2)
library(tidyverse)
library(vegan)

# Record R and package versions.
R.version.string
packageVersion("dada2")
packageVersion("tidyverse")
packageVersion("vegan")

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Define the folder containing demultiplexed FASTQ files.
path <- "~/OneDrive - Texas A&M University/Lassen2017/Lassen18S"

# List files in the input folder.
list.files(path)

# Import paired-end FASTQ filenames.
fnF <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names=TRUE))
fnR <- sort(list.files(path, pattern="_L001_R2_001.fastq", full.names=TRUE))

# Extract sample names from FASTQ filenames.
sample_names <- sapply(strsplit(basename(fnF), "_"), `[`, 1)

# Plot forward and reverse read quality profiles.
plotQualityProfile(fnF[1:20])
plotQualityProfile(fnR[1:20])

# Create filtered read output filenames.
filtF <- file.path(path, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtR <- file.path(path, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtF) <- sample_names
names(filtR) <- sample_names

# Filter and trim paired-end reads.
filter_out <- filterAndTrim(fnF, filtF, fnR, filtR, truncLen=c(200,190), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE)

# Inspect the filtering summary.
head(filter_out)

# Learn error rates.
errF <- learnErrors(filtF, multithread=TRUE)
errR <- learnErrors(filtR, multithread=TRUE)

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
dadaF <- dada(derepF, err=errF, multithread=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE)

# Inspect ASV inference for the first sample.
dadaF[[1]]

# Merge paired-end reads.
mergers <- mergePairs(dadaF, derepF, dadaR, derepR, verbose=TRUE, maxMismatch=0)

# Inspect merged reads for the first sample.
head(mergers[[1]])

# Create the ASV table.
seqtab_18S <- makeSequenceTable(mergers)

# Check ASV table dimensions.
dim(seqtab_18S)

# Inspect ASV sequence length distribution.
table(nchar(getSequences(seqtab_18S)))

# Remove chimeric ASVs.
seqtab_18S_nochim <- removeBimeraDenovo(seqtab_18S, method="consensus", multithread=TRUE, verbose=TRUE)

# Check ASV table dimensions after chimera removal.
dim(seqtab_18S_nochim)

# Calculate the proportion of reads retained after chimera removal.
sum(seqtab_18S_nochim)/sum(seqtab_18S)

# Track read retention.
getN <- function(x) sum(getUniques(x))
track_18S <- cbind(filter_out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(mergers, getN), rowSums(seqtab_18S_nochim))
colnames(track_18S) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_18S) <- sample_names

# Inspect read retention table.
head(track_18S)

# Save read retention table.
write.csv(track_18S, "Lassen_18S_DADA2_read_tracking.csv")

# Save final 18S ASV count table.
write.csv(seqtab_18S_nochim, "Lassen_18S_ASV_table.csv")

# Assign 18S taxonomy using PR2.
tax_18S <- assignTaxonomy(seqtab_18S_nochim, "data/references/pr2_version_5.1.1_SSU_dada2.fasta.gz", multithread=TRUE, minBoot=95, verbose=TRUE, taxLevels=c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"))

# Inspect taxonomy table.
tax_18S_print <- tax_18S
rownames(tax_18S_print) <- NULL
head(tax_18S_print)

# Save 18S taxonomy table.
write.csv(tax_18S, "Lassen_18S_ASV_taxonomy.csv")

