# Shannon diversity and temperature regression for 18S #####

# Load required packages.
library(tidyverse)

# Set the working directory.
setwd("~/OneDrive - Texas A&M University/Lassen2017/Lassen18S")

# Load Shannon diversity tables.
rare_rich_18S <- read.csv("Shannon_metadata_18S_rarefied.csv")
nonrare_rich_18S <- read.csv("Shannon_metadata_18S_nonrarefied.csv")

# Plot rarefied Shannon diversity against temperature.
temp_line_18S <- ggplot(rare_rich_18S, aes(x = Temp, y = Shannon)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "orange", linetype = "dashed") +
  labs(x = "Temperature (°C)", y = "Shannon Diversity Index (H')") +
  theme_minimal(base_size = 12) +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent"))

# View plot.
temp_line_18S

# Save plot for cowplot.
saveRDS(temp_line_18S, "18S_shannon_temp.rds")

# Run rarefied temperature regression.
model_18S_temp <- lm(Shannon ~ Temp, data = rare_rich_18S)
summary(model_18S_temp)
nobs(model_18S_temp)

# Run non-rarefied temperature regression.
model_18S_temp_nonrare <- lm(Shannon ~ Temp, data = nonrare_rich_18S)
summary(model_18S_temp_nonrare)
nobs(model_18S_temp_nonrare)