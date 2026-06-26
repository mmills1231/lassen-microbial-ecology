# Lassen Microbial Ecology

This repository contains the R scripts and analysis workflow supporting the manuscript: **Prokaryotic and Microbial Eukaryotic Communities Across Acid-Sulfate and Chloride-Rich Hot Springs in Lassen Volcanic National Park** which is currently under review at **Scientific Reports**.

---

## Overview

This study assessed prokaryotic (16S rRNA gene) and microbial eukaryotic (18S rRNA gene) communities across acid-sulfate and chloride-rich geothermal environments within Lassen Volcanic National Park, CA, USA. Analyses examined microbial diversity, community composition, environmental relationships, predicted functional potential, and microbial co-occurrence networks.

---

## Data availability

Raw amplicon sequencing data are publicly available through the **NCBI Sequence Read Archive (SRA)** under **BioProject PRJNA1432923**.

This repository contains the R scripts and processed input data required to reproduce the analyses presented in the manuscript, including ASV count tables, taxonomy tables, and sample metadata. 
Scripts are numbered and intended to be run sequentially within each analysis folder.

---

## Software

### R (v4.6.0) using the following packages:

* dada2
* phyloseq
* tidyverse
* ggplot2
* vegan
* plotly
* car
* pairwiseAdonis
* NetCoMi
* ggradar
* cowplot

### Python

Predicted prokaryotic functional annotation was performed using **FAPROTAX v1.2.12** with **Python 3.14.6**.

The complete workflow is provided in:

```
scripts/01_16S_analysis/06_FAPROTAX_protocol.md
```
