# Goal of this script: quantify differences in the percent that were read and counted

# Step 1: READ in htseq-count feature count files; store as a list of 4 lists of data.frames
# Step 2: Convert list of counts to 4-member list ("df") of 1 data frame per feature type
# Step 3: Separate counts of mapped reads and unmapped reads into separate tables
# Step 4: Create table summarizing #reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
# Step 5: Create table summarizing %reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
# Step 6: Sanity checks

##########################################################################################################
# Step 1: READ in htseq-count feature count files; store as a list of 4 lists of data.frames
##########################################################################################################

# user input:
wd <- "~/Desktop/Metagenome_TS/3_count/"
count_dir <- "all_counts"

#load libraries
library(plyr)

# create a list of lists of count_files
feature_type <- c("CDS", "rr", "rRNA","tRNA")

# set wd
setwd(wd)

list.count.files <- function(x) list.files(count_dir)[grep(x,list.files(count_dir))]
count_files_list <-  sapply(feature_type, list.count.files, simplify = FALSE)
length(count_files_list) #4, perfect
length(count_files_list[[1]]) #60, perfect
names(count_files_list) <- feature_type 
head(count_files_list[["CDS"]]) #name lines up with content of each list of count files

setwd(count_dir)

# read in all files in the list of files for each feature (4 features, 60 files to read per feature)
read.apply <- function (ft) counts[[ft]] <- sapply(count_files_list[[ft]], read.table, simplify=F)
counts <- sapply(feature_type, function(x) NULL)
counts <- sapply(feature_type, read.apply, simplify=F)
setwd(wd)
#save(counts,file="counts.RData")

##########################################################################################################
# Step 2: Convert list of counts to 4-member list ("df") of 1 data frame per feature type
##########################################################################################################

#setwd(wd)
#load("counts.RData")

# write function to rename a data.frame's rows and remove the now-redundant first column

rename.rows <- function(x) { # where x is a file
	data.frame(x)
	row.names(x) <- x$V1
	x$V1 <- NULL
	return(x)
	}

make.df <- function(ft) {
	counts[[ft]] <- sapply(counts[[ft]], rename.rows, simplify=F)
	for (i in 1:length(counts[[ft]])) {
	names(counts[[ft]][[i]]) <- substr(names(counts[[ft]][i]), nchar(ft)+2,nchar(names(counts[[ft]][i]))-6)
	}
	df[[ft]] <- as.data.frame(do.call(cbind,counts[[ft]]))
}
df <- list()
df <- sapply(feature_type, make.df, simplify=F)

# save(df,file="counts_df.RData")

##########################################################################################################
# Step 3: Separate counts of mapped reads and unmapped reads into separate tables
##########################################################################################################
# setwd(wd)
# load("counts_df.RData")

# total sums of reads for each metagenome:
df_totals <- sapply(feature_type, function(ft) as.data.frame(apply(df[[ft]],2,sum)), simplify=F)
# notice that the "total" values for CDS, rr, rRNA, and tRNA are the same; this represents the TOTAL number of reads input into mapping.
total_reads <- df_totals[[1]]
names(total_reads) <- "total_reads"

# CREATE A DATA FRAME OF ALL THE READS NOT MAPPED PER METAGENOME
select.unmapped <- function(ft) df[[ft]][((nrow(df[[ft]])-4):nrow(df[[ft]])),]
unmapped <- sapply(feature_type, select.unmapped, simplify=F)
unmapped <- sapply(feature_type, function(ft) t(unmapped[[ft]]), simplify=F)

# CREATE A DATA FRAME OF ALL THE READS MAPPED PER METAGENOME
remove.unmapped <- function(ft) df[[ft]][-((nrow(df[[ft]])-4):nrow(df[[ft]])),]
mapped <- sapply(feature_type, remove.unmapped, simplify=F)

# save(mapped, file="mapped_counts.RData")
# save(unmapped, file="unmapped_counts.RData")

#######################################################################################################
# Step 4: Create table summarizing #reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
#######################################################################################################
# setwd(wd)
# load("mapped_counts.RData")
# load("unmapped_counts.RData")

# CREATE A DATA FRAME OF the SUM of READS MAPPED PER METAGENOME
mapped_sums <- sapply(feature_type, function(ft) as.data.frame(apply(mapped[[ft]],2,sum)), simplify=F)

create.summary.table <- function(ft) {
	x <- cbind(unmapped[[ft]], 
	sum_mapped=mapped_sums[[ft]],
	#total_reads=sapply(feature_type, function(ft) as.data.frame(apply(df[[ft]],2,sum)), simplify=F)[,1] )
	total_reads=as.data.frame(apply(df[[ft]],2,sum)) )
	names(x) <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique", "sum_mapped", "total_reads")
	return(x)
}

summary_table <- sapply(feature_type,create.summary.table, simplify=F)

summary_numbers <- data.frame(sapply(feature_type, function(ft) apply(summary_table[[ft]],2,sum) ))

#######################################################################################################
# Step 5: Create table summarizing %reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
#######################################################################################################

summary_percents <- cbind(
pct_CDS=(summary_numbers$CDS[which(rownames(summary_numbers)=="sum_mapped")]/summary_numbers$CDS[which(rownames(summary_numbers)=="total_reads")])*100,
pct_ambiguous_CDS=(summary_numbers$CDS[which(rownames(summary_numbers)=="ambiguous")]/summary_numbers$CDS[which(rownames(summary_numbers)=="total_reads")])*100,
pct_rr=(summary_numbers$rr[which(rownames(summary_numbers)=="sum_mapped")]/summary_numbers$rr[which(rownames(summary_numbers)=="total_reads")])*100,
pct_rRNA=(summary_numbers$rRNA[which(rownames(summary_numbers)=="sum_mapped")]/summary_numbers$rRNA[which(rownames(summary_numbers)=="total_reads")])*100,
pct_tRNA=(summary_numbers$tRNA[which(rownames(summary_numbers)=="sum_mapped")]/summary_numbers$tRNA[which(rownames(summary_numbers)=="total_reads")])*100,
pct_not_aligned=(summary_numbers$CDS[which(rownames(summary_numbers)=="not_aligned")]/summary_numbers$CDS[which(rownames(summary_numbers)=="total_reads")])*100,

pct_aligned_noncoding=
((summary_numbers$CDS[which(rownames(summary_numbers)=="total_reads")]-
(summary_numbers$CDS[which(rownames(summary_numbers)=="sum_mapped")] +
 summary_numbers$CDS[which(rownames(summary_numbers)=="ambiguous")] + 
 summary_numbers$rr[which(rownames(summary_numbers)=="sum_mapped")] + 
 summary_numbers$rRNA[which(rownames(summary_numbers)=="sum_mapped")] + 
 summary_numbers$tRNA[which(rownames(summary_numbers)=="sum_mapped")] + 
 summary_numbers$CDS[which(rownames(summary_numbers)=="not_aligned")]))/summary_numbers$CDS[which(rownames(summary_numbers)=="total_reads")])*100) #417413506

print(summary_percents, digits=3)

#######################################################################################################
# Step 6: Sanity checks
#######################################################################################################

# check: does rowSums of the first 7 columns equal the total_reads column (total reads input into the mapping)? It should.
#sum(which(rowSums(summary_table[["CDS"]][,-ncol(summary_table[["CDS"]])])==summary_table[["CDS"]][,ncol(summary_table[["CDS"]])])==F) # 0
#sum(which(rowSums(summary_table[["rr"]][,-ncol(summary_table[["rr"]])])==summary_table[["rr"]][,ncol(summary_table[["rr"]])])==F) # 0
#sum(which(rowSums(summary_table[["rRNA"]][,-ncol(summary_table[["rRNA"]])])==summary_table[["rRNA"]][,ncol(summary_table[["rRNA"]])])==F) # 0
#sum(which(rowSums(summary_table[["tRNA"]][,-ncol(summary_table[["tRNA"]])])==summary_table[["tRNA"]][,ncol(summary_table[["tRNA"]])])==F) # 0
# And it does. Good. 

# Note: not_aligned column is the same for all 4 feature types.
# The sum of the following does not sum to equal the total number of reads.
# this tells us that there are reads that align but not to regions categorized as CDS, rr, tRNA, rRNA - maybe junk
# summary_table[["CDS"]][,"ambiguous"] + 
# summary_table[["CDS"]][,"not_aligned"] + 
# summary_table[["CDS"]][,"sum_mapped"] + 
# summary_table[["rr"]][,"sum_mapped"] + 
# summary_table[["tRNA"]][,"sum_mapped"] + 
# summary_table[["rRNA"]][,"sum_mapped"]
