# Objective: annotate the locus tags

######################################################################
# load required filepaths, libraries, and data
######################################################################

#library(plyr)
count_dir <- "../3_count"

normalize_dir <- "../4_normalize"

annotate_dir <- "../5_annotate"

analysis_dir <- "../6_analysis"

data_dir <- "../raw_data"

output_dir <- paste(site,minid,"output",sep="_")

dir.create(output_dir)

setwd(data_dir)

cog_file <- list.files(site)[grep(".cog.txt",list.files(site))]
ko_file <- list.files(site)[grep(".ko.txt",list.files(site))]
ec_file <- list.files(site)[grep(".ec.txt",list.files(site))]
gene_product_file <- list.files(site)[grep(".gene_product.txt",list.files(site))]
pfam_file <- list.files(site)[grep(".pfam.txt",list.files(site))]

setwd(site)
#cog <-read.delim("3300000553.a.cog.txt", sep="\t", header=F) 

scheme <- list()

scheme[["cog"]] <-read.delim(cog_file, sep="\t", header=F)
names(scheme[["cog"]]) <- c("gene_id","ID","percent_identity", "align_length","query_start","query_end","subj_start","subj_end","evalue","bit_score")

scheme[["ko"]] <-read.delim(ko_file, sep="\t", header=F)
names(scheme[["ko"]]) <- c("gene_id", "img_ko_flag", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")

scheme[["ec"]] <-read.delim(ec_file, sep="\t", header=F)
names(scheme[["ec"]]) <- c("gene_id", "img_ko_flag", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")

scheme[["gene_product"]] <-read.delim(gene_product_file, sep="\t", header=F)
names(scheme[["gene_product"]]) <- c("gene_id", "ID", "source")

# depth <-read.delim("3300000439.a.depth.txt", sep="\t", header=F)
# names(depth) <- c("orig_id", "depth")

scheme[["pfam"]] <- read.delim(pfam_file, sep="\t", header=F)
names(scheme[["pfam"]]) <- c("gene_id", "ID", "percent_identity", "query_start", "query_end", "subj_start", "subj_end", "evalue", "bit_score", "align_length")
setwd("..")


# select preferred annotator
ann <- scheme[[annotation]]

######################################################################
# analysis section: how many locus tags have multiple annotations?
######################################################################

# load the RPKM matrix data set
setwd(normalize_dir)
load(paste(output_dir, "RPKM.RData", sep="/"))
setwd(annotate_dir)

# create a data frame out of the RPKM matrix - for the locus tags
RPKM_lt <- as.data.frame(RPKM) #dim(RPKM_lt) [1] 352330     47

# create a column for the locus tag names
RPKM_lt$locus_tag <- as.factor(rownames(RPKM))

# merge in the annotation IDs
RPKM_lt_ann <- merge(RPKM_lt, subset(ann, select=c("gene_id","ID","evalue")), by.x="locus_tag", by.y="gene_id", all=T) #[1] 365417     49

# by this merge, 13087 locus tags are duplicated; this means that 3.7% of the genes have multiple annotations #dim(RPKM_lt_ann) 

# remove the rows for the genes with no annotation (is this necessary though)
RPKM_ann <- RPKM_lt_ann[-which(is.na(RPKM_lt_ann$ID)),]  # dim(RPKM_ann) [1] 191404     49
#RPKM_ann <- RPKM_lt_ann

# length(which(duplicated(RPKM_ann$locus_tag)==T))

#identify which cogs have multiple annotations
# geneID_cog_counts <- ddply(cog, "gene_id", summarise, count=length(ID)) # takes ~90 seconds
# dim(dup_cog) # 366726      2

# some cogs have 10 annotations! 
# max(geneID_cog_counts$count) #[1] 10
# This is because some locus_tags were encountered in multiple locations in the coassembly
# (if you look at their start and stop locations listed in the COG annotation file)

# dup_cog <- subset(geneID_cog_counts, count>1)

# how many genes have multiple annotations?
# dim(dup_cog) # 11032     2

# what is the total number of annotations for these multiply-annotated genes?
# sum(dup_cog$count) # 24119


######################################################################
#### Here I select the best annotation for a locus tag
######################################################################

# identify the cogs with the best evalue
ptm <- proc.time()
best_evalue <- ddply(ann, "gene_id", summarise, evalue=min(evalue)) #42s, dim(best_evalue) = 178317      2
proc.time() - ptm

# retain only those with lowest evalue in the cog annotation file
ptm <- proc.time()
ann_best_evalue <- merge(ann, best_evalue, all.y=T)
proc.time() - ptm

# of those, identify the cogs with the highest bit score
ptm <- proc.time()
best_bit <- ddply(ann_best_evalue, "gene_id", summarise, bit_score=max(bit_score)) # 43s
proc.time() - ptm

# retain only those with highest bit score in the annotation file
ptm <- proc.time()
ann_best_bit <- merge(ann_best_evalue, best_bit, all.y=T) # 3s, dim(ann_best_bit) = 178335     10
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
ann_best_score <- ann_best_score[-duplicated(ann_best_score$gene_id),]


setwd(annotate_dir)
save(ann_best_score, file=paste(output_dir,"best_annotation.RData", sep="/"))


######################################################################
# annotate each annotatable gene
######################################################################

# load the RPKM matrix data set
setwd(normalize_dir)
load(paste(output_dir,"RPKM.RData", sep="/"))
setwd(annotate_dir)
load(paste(output_dir, "best_annotation.RData", sep="/")) #ann_best_score

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

# transpose (COG) table to make it amenable to merging in dates
row.names(RPKM_ann_collapsed) <- RPKM_ann_collapsed$ID
RPKM_ann_collapsed$ID <- NULL
RPKM_ann_t <- data.frame(t(RPKM_ann_collapsed))
RPKM_ann_t$Library.ID <- substr(row.names(RPKM_ann_t), 1,5)

# assign dates to metagenomes
setwd(count_dir)
load(paste(output_dir,"dates.RData", sep="/"))

# merge in the dates; 
RPKM_dates <- merge(dates, RPKM_ann_t, all.y=T)

# save for analysis of COG values for duplicate dates
setwd(annotate_dir)
save(RPKM_dates, file=paste(output_dir,"RPKM_all_metagenomes.RData", sep="/"))

# collapse (COG) RPKM values to generate 1 RPKM value per date
# i.e., average RPKM values for COGs on the same date
RPKM_ts <- ddply(subset(RPKM_dates,select=-Library.ID), .(date), numcolwise(mean))

# merge back in the Library.ID values
RPKM_ts <- merge(dates,RPKM_ts, all.y=T)

# save time series as .RData file and CSV for downstream analysis 
setwd(annotate_dir)
save(RPKM_ts, file=paste(output_dir,"RPKM_ts.RData", sep="/"))
write.csv(RPKM_ts,paste(output_dir,"RPKM_time_series.csv", sep="/"))
