# Plot 16S and 18S Shannon diversity by location and sample type #####

# Load required packages.
library(tidyverse)
library(cowplot)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S18S")

# Load rarefied Shannon diversity tables.
richness_16S <- read.csv("Shannon_metadata_16S_rarefied.csv")
richness_18S <- read.csv("Shannon_metadata_18S_rarefied.csv")

# Add dataset labels.
richness_16S$Dataset <- "Prokaryotes"
richness_18S$Dataset <- "Eukaryotes"

# Combine 16S and 18S Shannon diversity tables.
combined_richness <- bind_rows(richness_16S, richness_18S)

# Order sample types along the hydrothermal gradient.
combined_richness$Sampletype <- factor(combined_richness$Sampletype, levels = c("spring fluid", "precipitate", "crusted sediment", "sediment", "biofilm", "stream"))

# Plot aesthetics #####

# Define dataset colors.
dataset_colors <- c("Prokaryotes" = "#5D47BD", "Eukaryotes" = "#DF9A00")
dataset_fills <- c("Prokaryotes" = "#785EF0", "Eukaryotes" = "#FFB000")

# Define shared alpha diversity theme.
alpha_theme <- theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank(),
        axis.text.y.left = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 10, angle = 45, colour = "black", vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 12, colour = "black"),
        title = element_text(size = 12))

# Plot Shannon diversity by location.
loc_plot <- ggplot(combined_richness, aes(x = Location, y = Shannon, fill = Dataset, color = Dataset)) +
  geom_boxplot(alpha = 0.1, lwd = 0.65) +
  geom_jitter(aes(color = Dataset), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1) +
  labs(x = "Location", y = "Shannon Diversity Index (H')") +
  scale_colour_manual(values = dataset_colors, name = "Dataset") +
  scale_fill_manual(values = dataset_fills, name = "Dataset") +
  guides(color = "none", fill = guide_legend(override.aes = list(alpha = 0.3))) +
  alpha_theme

# Plot Shannon diversity by sample type.
type_plot <- ggplot(combined_richness, aes(x = Sampletype, y = Shannon, fill = Dataset, color = Dataset)) +
  geom_boxplot(alpha = 0.1, lwd = 0.65) +
  geom_jitter(aes(color = Dataset), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1) +
  labs(x = "Sample type", y = "Shannon Diversity Index (H')") +
  scale_colour_manual(values = dataset_colors, name = "Dataset") +
  scale_fill_manual(values = dataset_fills, name = "Dataset") +
  guides(color = "none", fill = guide_legend(override.aes = list(alpha = 0.3))) +
  alpha_theme

# Combine plots with one shared legend.
shared_legend <- get_legend(type_plot + theme(legend.position = "bottom"))
type_noleg <- type_plot + theme(legend.position = "none")
loc_noleg <- loc_plot + theme(legend.position = "none")

combined_plot <- plot_grid(type_noleg, loc_noleg, labels = c("A", "B"), ncol = 1, align = "v")
final_plot <- plot_grid(combined_plot, shared_legend, ncol = 1, rel_heights = c(1, 0.08))

# View combined plot.
final_plot
ggsave("combined_shannon_location_sampletype.pdf", plot = final_plot, width = 8, height = 10)