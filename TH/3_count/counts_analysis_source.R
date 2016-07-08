# counts_analysis_source.R

# user input
wd <- "~/Desktop/Metagenome_TS/3_count"
feature_type <- c("CDS", "rr", "rRNA","tRNA")

# load libraries
library(plyr)
library(lattice)
library(ggplot2)
library(reshape2)
library(cowplot)
library(grDevices)

# load required data
setwd(wd)
load("mapped_counts.RData") #
load("unmapped_counts.RData")
load("total_read_counts.RData")