# Plot Shannon diversity by pH regime for 16S and 18S #####

# Load required packages.
library(tidyverse)
library(car)
library(cowplot)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S18S")

# Load rarefied Shannon diversity tables.
rich_16S <- read.csv("Shannon_metadata_16S_rarefied.csv")
rich_18S <- read.csv("Shannon_metadata_18S_rarefied.csv")

# Add pH regime labels.
add_pH_regime <- function(df) {
  df$pH_regime <- NA
  df$pH_regime[df$pH <= 5.5] <- "Acidic"
  df$pH_regime[df$pH > 5.5] <- "Circumneutral"
  df$pH_regime[is.na(df$pH_regime) & df$Location == "Bumpass Hell"] <- "Acidic"
  df$pH_regime[is.na(df$pH_regime) & df$Location == "Sulphur Works Acid"] <- "Acidic"
  df$pH_regime[is.na(df$pH_regime) & df$Location == "Sulphur Works Neutral"] <- "Circumneutral"
  df$pH_regime <- factor(df$pH_regime, levels = c("Acidic", "Circumneutral"))
  df <- df[!is.na(df$pH_regime), ]
  df
}

# Assign pH regimes.
rich_16S_pH <- add_pH_regime(rich_16S)
rich_18S_pH <- add_pH_regime(rich_18S)

# Check pH regime assignments.
table(rich_16S_pH$Location, rich_16S_pH$pH_regime)
table(rich_18S_pH$Location, rich_18S_pH$pH_regime)

# Plot aesthetics #####

# Define location colors.
location_colors <- c("Bumpass Hell" = "#FF3131",
                     "Devils Kitchen" = "#0097B2",
                     "Sulphur Works Neutral" = "#00BF63",
                     "Pilot Pinnacle" = "#CB6CE6",
                     "Sulphur Works Acid" = "#FFBD59")

# Define shared pH plot theme.
pH_theme <- theme_classic() +
  theme(axis.text.x = element_text(size = 12),
        panel.background = element_rect(fill = "transparent", color = "black", linewidth = 1),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Create function for pH regime plots.
make_pH_plot <- function(df, y_label) {
  ggplot(df, aes(x = pH_regime, y = Shannon, color = Location)) +
    geom_boxplot(fill = "white", color = "black", alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.85) +
    scale_colour_manual(values = location_colors) +
    labs(x = "pH regime", y = y_label, color = "Location") +
    pH_theme
}

# Create 16S and 18S pH regime plots.
pH_plot_16S <- make_pH_plot(rich_16S_pH, "Shannon Diversity Index (H')")
pH_plot_18S <- make_pH_plot(rich_18S_pH, "Shannon Diversity Index (H')")

# Combine plots with one shared legend.
shared_legend <- get_legend(pH_plot_16S + theme(legend.position = "bottom"))
pH_16S_noleg <- pH_plot_16S + theme(legend.position = "none")
pH_18S_noleg <- pH_plot_18S + theme(legend.position = "none")

combined_pH_plot <- plot_grid(pH_16S_noleg, pH_18S_noleg, labels = c("A", "B"), ncol = 1, align = "v")
final_pH_plot <- plot_grid(combined_pH_plot, shared_legend, ncol = 1, rel_heights = c(1, 0.08))

# View combined plot.
final_pH_plot