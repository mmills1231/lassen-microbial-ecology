# Create 16S relative abundance tables and bar plots #####

# Load required packages.
library(phyloseq)
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen16S")

# Load the filtered 16S phyloseq object.
physeq_16S_filtered <- readRDS("physeq_16S_filtered.rds")

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
genus_abun_16S <- make_abundance_table(physeq_16S_filtered, "Genus", "GenusAbundance16S.csv")
family_abun_16S <- make_abundance_table(physeq_16S_filtered, "Family", "FamilyAbundance16S.csv")
order_abun_16S <- make_abundance_table(physeq_16S_filtered, "Order", "OrderAbundance16S.csv")
class_abun_16S <- make_abundance_table(physeq_16S_filtered, "Class", "ClassAbundance16S.csv")

# Keep classes above 5% relative abundance within individual samples.
class_sample_16S <- class_abun_16S %>%
  select(Phylum, Class, Sample, Location, Abundance) %>%
  group_by(Phylum, Class, Sample, Location) %>%
  summarize(avg_abundance=mean(Abundance), .groups="drop") %>%
  filter(avg_abundance > 0.05) %>%
  mutate(Class=reorder(as.character(Class), as.character(Phylum), FUN=max))

# Keep classes above 5% relative abundance within sample types.
class_type_16S <- class_abun_16S %>%
  select(Phylum, Class, Sampletype, Abundance) %>%
  group_by(Phylum, Class, Sampletype) %>%
  summarize(avg_abundance=mean(Abundance), .groups="drop") %>%
  filter(avg_abundance > 0.05) %>%
  mutate(Class=reorder(as.character(Class), as.character(Phylum), FUN=max))

# Save filtered class abundance tables used for plotting.
write.csv(class_sample_16S, "ClassAbundance16S_samples_gt5percent.csv", row.names=FALSE)
write.csv(class_type_16S, "ClassAbundance16S_sampletypes_gt5percent.csv", row.names=FALSE)

#### Aesthetics ####

# Define class colors for 16S bar plots.
class_colors_16S <- c(
  "Acidobacteriota.Acidobacteriae"       = "#113b30",
  "Acidobacteriota.Blastocatellia"       = "#e4f3ff",
  "Actinobacteriota.Acidimicrobiia"      = "#262914",
  "Actinobacteriota.Actinobacteria"      = "#38571a",
  "Aquificota.Aquificae"                 = "#6a96e3",
  "Armatimonadota.Fimbriimonadia"        = "#d58400",
  "Bacteroidota.Bacteroidia"             = "#8ec8ff",
  "Bacteroidota.Ignavibacteria"          = "#b17242",
  "Campylobacterota.Desulfurellia"       = "#d9cafe",
  "Chloroflexi.Anaerolineae"             = "#646f1e",
  "Chloroflexi.Chloroflexia"             = "#4a128b",
  "Chloroflexi.Ktedonobacteria"          = "#fec700",
  "Cyanobacteria.Cyanobacteriia"         = "#8425df",
  "Desulfobacterota.Dissulfuribacteria"  = "#b17242",
  "Firmicutes.Bacilli"                   = "#6bb9c5",
  "Firmicutes.Clostridia"                = "#61177c",
  "Firmicutes.Desulfitobacteriia"        = "#7b219f",
  "Firmicutes.Sulfobacillia"             = "#C71585",
  "Firmicutes.Thermoanaerobacteria"      = "#ab0104",
  "Nitrospirota.Leptospirillia"          = "#ed6036",
  "Planctomycetota.Phycisphaerae"        = "#e4852b",
  "Planctomycetota.Planctomycetes"       = "#ffd9a8",
  "Proteobacteria.Alphaproteobacteria"   = "#669362",
  "Proteobacteria.Gammaproteobacteria"   = "#ecaf25",
  "Spirochaetota.Brachyspirae"           = "#4a128b",
  "Spirochaetota.Leptospirae"            = "#fff76b",
  "Spirochaetota.Spirochaetia"           = "#fff2d5",
  "Thermotogota.Thermotogae"             = "#cf5653",
  "Verrucomicrobiota.Chlamydiae"         = "#ffd9a8",
  "Crenarchaeota.Nitrososphaeria"        = "#FF6A00",
  "Crenarchaeota.Thermoprotei"           = "#B12104",
  "Thermoplasmatota.Thermoplasmata"      = "#ee131e",
  "Others.Others"                        = "grey80"
)

# Define a shared theme for 16S bar plots.
bar_theme_16S <- theme_bw() +
  theme(axis.text.y.left=element_text(size=12, color="black"),
        axis.text.x=element_text(size=10, angle=90, vjust=1.5, hjust=1.5, color="black"),
        axis.title.y=element_text(size=12, color="black"),
        title=element_text(size=12),
        strip.background=element_rect(linewidth=0.5, color="black", fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        panel.border=element_rect(linewidth=0.5, color="black"),
        plot.background=element_rect(fill="transparent", color=NA),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.background=element_rect(fill="transparent"))

#### Bar plots ####

# Plot class-level relative abundance by sample.
barplot_samples_16S <- ggplot(class_sample_16S) +
  geom_col(aes(x=Sample, y=avg_abundance, fill=interaction(Phylum, Class)), position="stack", show.legend=TRUE, color="black", size=0.5) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values=class_colors_16S) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  facet_grid(. ~ Location, scales="free", space="free") +
  bar_theme_16S

# View the sample-level bar plot.
barplot_samples_16S

# Save the sample-level bar plot.
ggsave("Barplot_samples_16S.pdf", plot=barplot_samples_16S, bg="transparent", width=16, height=8)

# Plot class-level relative abundance by sample type.
barplot_types_16S <- ggplot(class_type_16S) +
  geom_col(aes(x=Sampletype, y=avg_abundance, fill=interaction(Phylum, Class)), position="fill", show.legend=TRUE, color="black", size=0.5) +
  ylab("Relative Abundance") +
  xlab(NULL) +
  scale_fill_manual(values=class_colors_16S) +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0)) +
  bar_theme_16S

# View the sample-type bar plot.
barplot_types_16S

# Save the sample-type bar plot.
ggsave("Barplot_sampletypes_16S.pdf", plot=barplot_types_16S, bg="transparent", width=12, height=8)