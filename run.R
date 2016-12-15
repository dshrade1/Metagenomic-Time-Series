#user input 
ptm <- proc.time()
read_site <- function()
{
	site <- readline(prompt="Select sample site (TE, TH, or ME): ")
	return(as.character(site))
}
site <- read_site()

read_minid <- function()
{
	minid <- readline(prompt="Select minid (from 0.70 to 1.00): ")
	return(as.character(minid))
}
minid <- read_minid()

read_annotation <- function()
{
	annotation <- readline(prompt="Select preferred annotation scheme (COG, KO, EC, Pfam): ")
	return(as.character(annotation))
}
annotation <- read_annotation()

# load required libraries
library(plyr)
library(lattice)
library(ggplot2)
library(reshape2)
library(cowplot)
library(grDevices)

# read in counts
setwd("3_count")
source("count_features.R")
setwd("..")

# perform analysis of counts (optional)
#setwd("3_count")
#source("counts_analysis.R")
#setwd("..")

# normalize counts
setwd("4_normalize")
source("normalize_counts.R")
setwd("..")

# perform normalization analysis (optional)
#setwd("4_normalize")
#source("normalize_analysis.R")
#setwd("..")

# annotate normalized reads
setwd("5_annotate")
source("annotate_genes.R")

# to do yet:
# indicate which R packages need to be installed