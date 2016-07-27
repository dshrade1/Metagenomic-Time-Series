# Objective: 
# annotate the locus tags with COG annotations

######################################################################
# load required filepaths, libraries, and data
######################################################################

library(plyr)

# to do: have user select which annotation they would like to use #
# from ec, ko, cog, gene_product, pfam, etc.

# to do: have user select which annotation they would like to use #

normalize_dir <- "~/Desktop/Metagenome_TS/4_normalize"

annotate_dir <- "~/Desktop/Metagenome_TS/5_annotate"

annotation_dir <- "~/Desktop/Metagenome_TS/5_annotate/TH_coassembly_annotation"

analysis_dir <- "~/Desktop/Metagenome_TS/6_analysis"


setwd(annotation_dir)
cog <-read.delim("3300000553.a.cog.txt", sep="\t", header=F) 
names(cog) <- c("gene_id","ID","percent_identity", "align_length","query_start","query_end","subj_start","subj_end","evalue","bit_score")

ko <-read.delim("3300000553.a.ko.txt", sep="\t", header=F)
names(ko) <- c("gene_id", "img_ko_flag", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")

ec <-read.delim("3300000553.a.ec.txt", sep="\t", header=F)
names(ec) <- c("gene_id", "img_ko_flag", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")

gene_product <-read.delim("3300000553.a.gene_product.txt", sep="\t", header=F)
names(gene_product) <- c("gene_id", "ID", "source")

depth <-read.delim("3300000553.a.depth.txt", sep="\t", header=F)
names(depth) <- c("orig_id", "depth")

pfam <- read.delim("3300000553.a.pfam.txt", sep="\t", header=F)
names(pfam) <- c("gene_id", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")


# select your preferred annotator
ann <- cog

######################################################################
# analysis section: how many locus tags have multiple annotations?
######################################################################

# load the RPKM matrix data set
setwd(normalize_dir)
load("RPKM.RData")
setwd(annotate_dir)

# create a data frame out of the RPKM matrix - for the locus tags
RPKM_lt <- as.data.frame(RPKM) #dim(RPKM_lt) #[1] 690851     61

# create a column for the locus tag names
RPKM_lt$locus_tag <- as.factor(rownames(RPKM))

# merge in the annotation IDs
RPKM_lt_ann <- merge(RPKM_lt, subset(ann, select=c("gene_id","ID","evalue")), by.x="locus_tag", by.y="gene_id", all=T) #[1] 727660     63

# by this merge, 36809 locus tags are duplicated; this means that 5.3% of the genes have multiple annotations #dim(RPKM_lt_ann) 

# remove the rows for the genes with no annotation (is this necessary though)
RPKM_ann <- RPKM_lt_ann[-which(is.na(RPKM_lt_ann$ID)),]

length(which(duplicated(RPKM_ann$locus_tag)==T)) # 36809 # 9% of annotated genes have duplicate annotations (one locus tag has multiple annotations)

#identify which cogs have multiple annotations
ptm <- proc.time()
geneID_cog_counts <- ddply(cog, "gene_id", summarise, count=length(ID)) # takes ~90 seconds
proc.time() - ptm
# dim(dup_cog) # 366726      2

# some cogs have 11 annotations! 
max(geneID_cog_counts$count) #[1] 11
# This is because some locus_tags were encountered in multiple locations in the coassembly
# (if you look at their start and stop locations listed in the COG annotation file)

dup_cog <- subset(geneID_cog_counts, count>1)

# how many genes have multiple annotations?
# dim(dup_cog) # 29050     2

# what is the total number of annotations for these multiply-annotated genes?
# sum(dup_cog$count) # 65859


######################################################################
#### Here I select the best annotation for a locus tag
######################################################################

ann <- cog

# identify the cogs with the best evalue
ptm <- proc.time()
best_evalue <- ddply(ann, "gene_id", summarise, evalue=min(evalue)) #84s, dim(best_evalue) = 366726, 2
proc.time() - ptm

# retain only those with lowest evalue in the cog annotation file
ptm <- proc.time()
ann_best_evalue <- merge(ann, best_evalue, all.y=T)
proc.time() - ptm

# of those, identify the cogs with the highest bit score
ptm <- proc.time()
best_bit <- ddply(ann_best_evalue, "gene_id", summarise, bit_score=max(bit_score)) # 85s, dim(best_evalue) = 366726, 2
proc.time() - ptm

# retain only those with highest bit score in the annotation file
ptm <- proc.time()
ann_best_bit <- merge(ann_best_evalue, best_bit, all.y=T) # 3s, dim(ann_best_bit) = 366777     11
proc.time() - ptm

# of those, identify those with the highest %ID*align_length
ann_best_bit <- transform(ann_best_bit, align_score=percent_identity*align_length) 
ptm <- proc.time()
best_score <- ddply(ann_best_bit, "gene_id", summarise, align_score=max(align_score)) #92s, dim(best_score) = 366726, 2
proc.time() - ptm

# retain only those with highest align_score in the annotation file
ptm <- proc.time()
ann_best_score <- merge(ann_best_bit, best_score, all.y=T) # 
proc.time() - ptm

# One gene with multiple annotations remain, but they are for the same cog. Choose the first one.
ann_best_score <- ann_best_score[-duplicated(ann_best_score$gene_id),] # dim(ann_best_score) = 366726,    11


setwd(annotate_dir)
save(ann_best_score, file="best_annotation.RData")


######################################################################
# annotate each annotatable gene
######################################################################

# load the RPKM matrix data set
setwd(normalize_dir)
load("RPKM.RData")
setwd(annotate_dir)
load("best_annotation.RData") #ann_best_score

# create a data frame out of the RPKM matrix - for the locus tags
RPKM_lt <- as.data.frame(RPKM) #dim(RPKM_lt) #[1] 690851     61

# create a column for the locus tag names
RPKM_lt$gene_id <- as.factor(rownames(RPKM))

# merge in the annotation IDs
RPKM_lt_ann <- merge(RPKM_lt, subset(ann_best_score, select=c("gene_id","ID")), all.x=T) 

# remove the duplicate row
RPKM_lt_ann <- RPKM_lt_ann[!duplicated(RPKM_lt_ann),]

# remove the rows for the genes with no annotation (is this necessary though?)
RPKM_ann <- RPKM_lt_ann[-which(is.na(RPKM_lt_ann$ID)),]

######################################################################
# quantify total RPKM value per (COG) per date
######################################################################

# collapse annotated matrix to sum 1 RPKM value per (COG)
library(plyr)
RPKM_ann_collapsed <- ddply(RPKM_ann, .(ID), numcolwise(sum)) # dim(RPKM_ann_collapsed) = 4627   61 (same number of rows as cogs)
# row.names(RPKM_ann_collapsed) <- RPKM_ann_collapsed$ID
# RPKM_ann_collapsed$ID <- NULL
# which version of COG to use?


# transpose (COG) table to make it amenable to merging in dates
row.names(RPKM_ann_collapsed) <- RPKM_ann_collapsed$ID
RPKM_ann_collapsed$ID <- NULL
RPKM_ann_t <- data.frame(t(RPKM_ann_collapsed))
RPKM_ann_t$Library.ID <- substr(row.names(RPKM_ann_t), 1,5)

# assign dates to metagenomes
setwd(analysis_dir)
dates <- read.csv("TH_Metagenome_Dates.csv")
dates$date <- as.Date(dates$Date, format="%m/%d/%y")
dates$Date <- NULL

# merge in the dates; 
RPKM_dates <- merge(dates, RPKM_ann_t, all.y=T)

# save for analysis of COG values for duplicate dates
setwd(annotate_dir)
save(RPKM_dates, file="RPKM_all_metagenomes.RData")

# collapse (COG) RPKM values to generate 1 RPKM value per date
# i.e., average RPKM values for COGs on the same date
RPKM_ts <- ddply(subset(RPKM_dates,select=-Library.ID), .(date), numcolwise(mean)) #28s, dim(RPKM_ts) = 53, 4628

# merge back in the Library.ID values
RPKM_ts <- merge(dates,RPKM_ts, all.y=T)

# save time series as .RData file and CSV for downstream analysis 
setwd(annotate_dir) # or save in analysis directory?
save(RPKM_ts, file="RPKM_ts.RData")
write.csv(RPKM_ts,"RPKM_time_series.csv")
