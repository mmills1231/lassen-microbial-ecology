# Construct 16S and 18S Pearson correlation network #####

# Load required packages.
library(phyloseq)
library(NetCoMi)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S18S")

# Load rarefied phyloseq objects.
rare_16S <- readRDS("rare_16S.rds")
rare_18S <- readRDS("rare_18S.rds")

# Merge rarefied 16S and 18S phyloseq objects.
merged_16S_18S <- merge_phyloseq(rare_16S, rare_18S)

# Extract ASV table with samples as rows and ASVs as columns.
network_data <- t(as.data.frame(otu_table(merged_16S_18S)))

# Construct Pearson correlation network.
net_pearson <- netConstruct(network_data,
                            measure = "pearson",
                            normMethod = "clr",
                            zeroMethod = "multRepl",
                            sparsMethod = "threshold",
                            thresh = 0.3,
                            verbose = 3,
                            seed = 10000)

# Analyze network properties.
props_pearson <- netAnalyze(net_pearson, clustMethod = "cluster_fast_greedy")


# Save network objects.
saveRDS(net_pearson, "net_pearson_16S_18S.rds")
saveRDS(props_pearson, "props_pearson_16S_18S.rds")

# Plot Pearson correlation network.
plot(props_pearson,
     nodeColor = "cluster",
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 0.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = "",
     showTitle = TRUE,
     cexTitle = 2.3)

# View network summary.
summary(props_pearson)