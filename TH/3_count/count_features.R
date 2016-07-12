# Goal of this script: quantify differences in the percent of reads of each feature type

# Step 1: READ in htseq-count feature count files; store as a list of 4 lists of data.frames
# Step 2: Convert list of counts to 4-member list ("df") of 1 data frame per feature type
# Step 3: Separate counts of mapped reads and unmapped reads into separate tables

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
# Note: these are the only 4 features that appear in our GFF file. 
# i.e., this list is the output of levels(as.factor(gff$feature)) !! <- encode this?

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
setwd(wd)
# load("counts_df.RData")

# total sums of reads for each metagenome:
df_totals <- sapply(feature_type, function(ft) as.data.frame(apply(df[[ft]],2,sum)), simplify=F)
# notice that the "total" values for CDS, rr, rRNA, and tRNA are the same; this represents the TOTAL number of reads input into mapping.
total_reads <- df_totals[[1]]
names(total_reads) <- "total_reads"
save(total_reads, file="total_read_counts.RData")

# CREATE A DATA FRAME OF ALL THE READS NOT MAPPED PER METAGENOME
select.unmapped <- function(ft) df[[ft]][((nrow(df[[ft]])-4):nrow(df[[ft]])),]
unmapped <- sapply(feature_type, select.unmapped, simplify=F)
unmapped <- sapply(feature_type, function(ft) t(unmapped[[ft]]), simplify=F)

# CREATE A DATA FRAME OF ALL THE READS MAPPED PER METAGENOME
remove.unmapped <- function(ft) df[[ft]][-((nrow(df[[ft]])-4):nrow(df[[ft]])),]
mapped <- sapply(feature_type, remove.unmapped, simplify=F)

save(mapped, file="mapped_counts.RData")
save(unmapped, file="unmapped_counts.RData")
