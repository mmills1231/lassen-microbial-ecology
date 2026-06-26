# Create 18S relative abundance tables and bar plots #####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Load the filtered 18S phyloseq object.
physeq_18S_filtered <- readRDS("physeq_18S_filtered.rds")

# Create a function for relative abundance tables.
make_abundance_table <- function(ps_obj, tax_rank, output_file) {
  abundance_table <- ps_obj %>%
    tax_glom(taxrank = tax_rank) %>%
    transform_sample_counts(function(x) x / sum(x)) %>%
    psmelt() %>%
    arrange(.data[[tax_rank]])
  write.csv(abundance_table, output_file, row.names = FALSE)
  abundance_table
}

# Create relative abundance tables by taxonomic rank.
genus_abun_18S <- make_abundance_table(physeq_18S_filtered, "Genus", "GenusAbundance18S.csv")
family_abun_18S <- make_abundance_table(physeq_18S_filtered, "Family", "FamilyAbundance18S.csv")
order_abun_18S <- make_abundance_table(physeq_18S_filtered, "Order", "OrderAbundance18S.csv")
class_abun_18S <- make_abundance_table(physeq_18S_filtered, "Class", "ClassAbundance18S.csv")

# Group classes below 5% relative abundance as Others for individual samples.
class_sample_18S <- class_abun_18S %>%
  select(Division, Class, Sample, Location, Abundance) %>%
  group_by(Division, Class, Sample, Location) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Class = ifelse(avg_abundance < 0.05, "Others", as.character(Class)),
         Division = ifelse(Class == "Others", "Others", as.character(Division))) %>%
  group_by(Sample, Location, Division, Class) %>%
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop") %>%
  mutate(fill_group = paste(Division, Class, sep = "."))

# Group classes below 5% relative abundance as Others for sample types.
class_type_18S <- class_abun_18S %>%
  select(Division, Class, Sampletype, Abundance) %>%
  group_by(Division, Class, Sampletype) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(Class = ifelse(avg_abundance < 0.05, "Others", as.character(Class)),
         Division = ifelse(Class == "Others", "Others", as.character(Division))) %>%
  group_by(Sampletype, Division, Class) %>%
  summarize(avg_abundance = sum(avg_abundance), .groups = "drop") %>%
  mutate(fill_group = paste(Division, Class, sep = "."))

# Save filtered class abundance tables used for plotting.
write.csv(class_sample_18S, "ClassAbundance18S_samples_gt5percent.csv", row.names = FALSE)
write.csv(class_type_18S, "ClassAbundance18S_sampletypes_gt5percent.csv", row.names = FALSE)

# Plot aesthetics #####

# Define class colors for 18S bar plots.
class_colors_18S <- c(
  "Opisthokonta.Arthropoda"         = "#113b30",
  "Opisthokonta.Ascomycota"         = "#fb7156",
  "Stramenopiles.Bacillariophyceae" = "#b556b6",
  "Rhodophyta.Bangiophyceae"        = "#A3C400",
  "Opisthokonta.Basidiomycota"      = "#8687e9",
  "Stramenopiles.Bicoecea"          = "#D14285",
  "Discosea.Centramoebia"           = "#67c4f2",
  "Stramenopiles.Chrysophyceae"     = "#fab55c",
  "Alveolata.Colpodea"              = "#5E738F",
  "Opisthokonta.Cranatia"           = "#AD6F3B",
  "Discoba.Discoba_XX"              = "#599861",
  "Tubulinea.Elardia"               = "#5cfa66",
  "Discosea.Flabellinia"            = "#cdffea",
  "Opisthokonta.Gastrotricha"       = "#eacdff",
  "Discoba.Heterolobosea"           = "#0fcc90",
  "Alveolata.Litostomatea"          = "#00cae5",
  "Stramenopiles.Mediophyceae"      = "#535ab6",
  "Opisthokonta.Nematoda"           = "#fbc456",
  "Others.Others"                   = "grey80",
  "Alveolata.Oxyrrhea"              = "#C84248",
  "Centroplasthelida.Pterocystida"  = "#86b757",
  "Chlorophyta.Trebouxiophyceae"    = "#a0fa5c",
  "Evosea.Variosea"                 = "#CBD588",
  "Chlorophyta.Chlorophyceae"       = "#CD9BCD"
)

# Define a shared theme for 18S bar plots.
bar_theme_18S <- theme_bw() +
  theme(axis.text.y.left = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 1.5, hjust = 1.5, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        title = element_text(size = 12),
        strip.background = element_rect(linewidth = 0.5, color = "black", fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 0.5, color = "black"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))

# Bar plots #####

# Plot class-level relative abundance by sample.
barplot_samples_18S <- ggplot(class_sample_18S) +
  geom_col(aes(x = Sample, y = avg_abundance, fill = fill_group), position = "stack", show.legend = TRUE, color = "black", linewidth = 0.5) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = class_colors_18S, na.value = "grey80") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  facet_grid(. ~ Location, scales = "free", space = "free") +
  bar_theme_18S

# View the sample-level bar plot.
barplot_samples_18S

# Save the sample-level bar plot.
ggsave("Barplot_samples_18S.pdf", plot = barplot_samples_18S, bg = "transparent", width = 16, height = 8)

# Plot class-level relative abundance by sample type.
barplot_types_18S <- ggplot(class_type_18S) +
  geom_col(aes(x = Sampletype, y = avg_abundance, fill = fill_group), position = "fill", show.legend = TRUE, color = "black", linewidth = 0.5) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values = class_colors_18S, na.value = "grey80") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  bar_theme_18S

# View the sample-type bar plot.
barplot_types_18S

# Save the sample-type bar plot.
ggsave("Barplot_sampletypes_18S.pdf", plot = barplot_types_18S, bg = "transparent", width = 12, height = 8)