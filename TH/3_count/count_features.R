# Goal of this script: quantify differences in the percent that were read and counted

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
# THIS WORKED!! and is super fast! (well, relatively)
setwd(wd)
save(counts,file="counts.RData")

##########################################################################################################
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
# everything works perfectly!

# total sums reads of each metagenome:
df_totals <- sapply(feature_type, function(ft) apply(df[[ft]],2,sum), simplify=F)

# now, section off the last 5 rows of each metagenome's count table
### write a function "remove.unmapped" to section off the last 5 rows of 1 metagenome's count table

unmapped <- list()
remove.unmapped <- function(x) { #e.g., x= df[["CDS"]]
	unmapped[[ft]] <- df[[ft]][((nrow(df[[ft]])-4):nrow(df[[ft]])),]
	df[[ft]]       <- df[[ft]][-((nrow(df[[ft]])-4):nrow(df[[ft]])),]
}


# create a table of no_feature, ambiguous, too_low_aQual, not_aligned, and alignment_not_unique
select.unmapped <- function(ft) df[[ft]][((nrow(df[[ft]])-4):nrow(df[[ft]])),]
unmapped <- sapply(feature_type, select.unmapped, simplify=F)

# create a table of only mapped features
remove.unmapped <- function(ft) df[[ft]][-((nrow(df[[ft]])-4):nrow(df[[ft]])),]
mapped <- sapply(feature_type, remove.unmapped, simplify=F)

# now there is 1 table of mapped and one table of unmapped reads per feature type.

##########################################################################################################
##########################################################################################################

# Next, combine these data into a descriptive table
