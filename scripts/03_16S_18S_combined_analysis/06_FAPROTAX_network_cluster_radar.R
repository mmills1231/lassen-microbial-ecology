# Match network clusters to taxonomy and summarize FAPROTAX functions #####

# Load required packages.
library(phyloseq)
library(tidyverse)
library(ggradar)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S18S")

# Load rarefied 16S phyloseq object.
rare_16S <- readRDS("rare_16S.rds")

# Extract 16S taxonomy table.
taxonomy_16S <- as.data.frame(tax_table(rare_16S))
taxonomy_16S$ASV_ID <- rownames(taxonomy_16S)
taxonomy_16S <- taxonomy_16S %>% mutate(ASV_ID = gsub("^ASV", "", ASV_ID))

# Load network cluster assignments.
network_clusters <- read.csv("network_clusters_16S18S.csv")
network_clusters <- network_clusters %>% mutate(ASV_ID = gsub("^ASV", "", ASV_ID))

# Join network clusters with 16S taxonomy.
network_clusters_taxonomy <- network_clusters %>% left_join(taxonomy_16S, by = "ASV_ID") %>% select(Cluster_Color, ASV_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species)

# Save network cluster taxonomy table.
write.csv(network_clusters_taxonomy, "network_clusters_16S18S_taxonomy.csv", row.names = FALSE)

# Load FAPROTAX groups-to-records output.
groups2records <- read.delim("groups2records_lassen.txt", header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

# Format FAPROTAX output.
colnames(groups2records) <- as.character(groups2records[1, ])
groups2records <- groups2records[-1, ]

# Convert FAPROTAX output to long format.
faprotax_long <- groups2records %>% pivot_longer(cols = -record, names_to = "Function", values_to = "Present") %>% mutate(Present = as.numeric(Present)) %>% filter(Present == 1) %>% select(record, Function)

# Create genus-level taxonomy strings from the 16S taxonomy table.
taxonomy_faprotax <- taxonomy_16S %>% mutate(Genus = ifelse(is.na(Genus), "NA", Genus), record = paste(Kingdom, Phylum, Class, Order, Family, Genus, sep = ";"))

# Trim FAPROTAX records to genus level.
faprotax_long <- faprotax_long %>% mutate(record_genus = sapply(strsplit(record, ";"), function(x) paste(x[1:6], collapse = ";")))

# Match ASVs to FAPROTAX functions.
asv_functions <- taxonomy_faprotax %>% select(ASV_ID, record) %>% inner_join(faprotax_long, by = c("record" = "record_genus"))

# Check number of annotated ASVs.
length(unique(asv_functions$ASV_ID))

# Join network clusters and FAPROTAX functions.
cluster_function_table <- network_clusters_taxonomy %>% select(ASV_ID, Cluster = Cluster_Color) %>% inner_join(asv_functions, by = "ASV_ID")

# Collapse FAPROTAX functions into selected functional groups.
cluster_function_table <- cluster_function_table %>%
  mutate(Functional_Group = case_when(
    Function %in% c("aerobic_ammonia_oxidation", "nitrification", "nitrate_respiration", "nitrate_reduction", "nitrogen_respiration") ~ "nitrogen cycling",
    Function %in% c("sulfur_respiration", "thiosulfate_respiration", "respiration_of_sulfur_compounds", "dark_sulfur_oxidation", "dark_thiosulfate_oxidation", "dark_oxidation_of_sulfur_compounds") ~ "sulfur cycling",
    Function %in% c("photosynthetic_cyanobacteria", "oxygenic_photoautotrophy", "photoautotrophy", "phototrophy") ~ "phototrophy",
    Function == "dark_hydrogen_oxidation" ~ "dark_hydrogen_oxidation",
    Function %in% c("aerobic_chemoheterotrophy", "chemoheterotrophy") ~ "chemoheterotrophy",
    Function == "fermentation" ~ "fermentation",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Functional_Group))

# Summarize functional groups by network cluster.
cluster_summary <- cluster_function_table %>% group_by(Cluster, Functional_Group) %>% summarize(n_ASVs = n_distinct(ASV_ID), .groups = "drop") %>% group_by(Cluster) %>% mutate(Percent = 100 * n_ASVs / sum(n_ASVs)) %>% arrange(Cluster, desc(Percent))

# Save cluster functional summary.
write.csv(cluster_summary, "network_clusters_16S18S_FAPROTAX_summary.csv", row.names = FALSE)

# Prepare radar plot data.
radar_df <- cluster_summary %>%
  select(Cluster, Functional_Group, Percent) %>%
  pivot_wider(names_from = Functional_Group, values_from = Percent, values_fill = 0) %>%
  rename(`nitrogen\ncycling` = `nitrogen cycling`,
         `sulfur\ncycling` = `sulfur cycling`,
         phototrophy = phototrophy,
         `dark hydrogen\noxidation` = dark_hydrogen_oxidation,
         chemoheterotrophy = chemoheterotrophy,
         fermentation = fermentation)

# Convert percentages to proportions.
radar_df[, -1] <- radar_df[, -1] / 100
radar_df[is.na(radar_df)] <- 0

# Define cluster colors.
cluster_colors <- c("Blue" = "#447ACA", "Green" = "#44E187", "Purple" = "#A55CB6", "Red" = "#D14242", "Yellow" = "#D1D049")

# Create radar plot.
cluster_radar <- ggradar(radar_df, group.line.width = 1.2, group.point.size = 2, grid.min = 0, grid.mid = 0.5, grid.max = 1, values.radar = c("0%", "50%", "100%"), legend.position = "bottom", axis.label.size = 4, grid.label.size = 3) +
  scale_color_manual(values = cluster_colors) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5), text = element_text(size = 12))

# View radar plot.
cluster_radar

# Save radar plot.
ggsave("network_clusters_16S18S_FAPROTAX_radar.pdf", plot = cluster_radar, width = 8, height = 6, bg = "transparent")