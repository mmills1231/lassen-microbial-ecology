# Export 16S and 18S network cluster assignments #####

# Load required packages.
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S18S")

# Load analyzed Pearson network object.
props_pearson <- readRDS("props_pearson_16S_18S.rds")

# Extract cluster assignments from NetCoMi object.
network_clusters_16S18S <- data.frame(ASV_ID = names(props_pearson$clustering), Cluster = props_pearson$clustering)

# Convert cluster numbers to cluster color names.
network_clusters_16S18S <- network_clusters_16S18S %>%
  mutate(Cluster_Color = case_when(
    Cluster == 1 ~ "Blue",
    Cluster == 2 ~ "Green",
    Cluster == 3 ~ "Purple",
    Cluster == 4 ~ "Red",
    Cluster == 5 ~ "Yellow",
    TRUE ~ paste0("Cluster_", Cluster)
  ))

# Save network cluster assignments.
write.csv(network_clusters_16S18S, "network_clusters_16S18S.csv", row.names = FALSE)