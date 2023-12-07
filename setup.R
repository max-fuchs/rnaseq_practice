if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")

# Vector of all packages
packages <- c("GEOquery", "DESeq2", "ggplot2", "dplyr", "gplots", "pheatmap", 
              "EnhancedVolcano", "EDASeq", "biomaRt", "tibble", "writexl", 
              "clusterProfiler", "org.Hs.eg.db", "tidyr", "ReactomePA")

# Install packages
BiocManager::install(packages)

library(GEOquery)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(gplots)
library(pheatmap)
library(EnhancedVolcano)
library(EDASeq)
library(biomaRt)
library(tibble)
library(writexl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(ReactomePA)
