################ Objective for this script ###########################

# normalize read counts by RPKM

# RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000)
# where
# numReads = number of reads that mapped as a particular gene
# geneLength = length of that particular gene
# totalNumReads = total #reads in a metagenome that aligned to coassembly

######################################################################
# Step 1: load and set filepaths
######################################################################

count_dir <- "../3_count"
normalize_dir <- "../4_normalize"
data_dir <- paste("../raw_data", site, sep="/")
gff_file <-list.files(data_dir)[grep(".gff",list.files(data_dir))]
output_dir <- paste(site,minid,"output",sep="_")

setwd(normalize_dir)
dir.create(output_dir)

######################################################################
# Step 2: Create RPKM dividend (numerator)
######################################################################
setwd(count_dir)
load(paste(output_dir,"mapped_counts.RData", sep="/")) #mapped

# create the RPKM dividend by row-binding the matrices of all features
RPKMdividend <- as.matrix(do.call(rbind,unname(mapped)))
setwd(normalize_dir)
save(RPKMdividend, file=paste(output_dir,"RPKM_dividend.RData", sep="/"))

######################################################################
# Step 3: Create RPKM divisor (denominator)
######################################################################

# create a function to read a GFF file into R as a data.frame
gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

# create function to extract an GFF file attribute field as a column in an R data.frame
# the following function is from https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html
getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gff <- gffRead(paste(data_dir, gff_file, sep="/"))
gff$locus_tag <- getAttributeField(gff$attributes, "locus_tag") # slow
gff$length <- gff$end - gff$start

# 3.a: create a data.frame that describes the locus tag length
tag_lengths <- gff[c(which(names(gff)=="locus_tag"), which(names(gff)=="length"))]
row.names(tag_lengths) <- tag_lengths$locus_tag
tag_lengths$locus_tag <- NULL
tag_lengths <- as.matrix(tag_lengths)

# divide locus_tag length by 1000
tag_lengths <- tag_lengths/1000

# 3.b: create a table that represents the number of reads mapped per metagenome.
setwd(count_dir)
load(paste(output_dir,"unmapped_counts.RData", sep="/")) # unmapped
load(paste(output_dir,"total_read_counts.RData", sep="/")) # total_reads
setwd(normalize_dir)
not_aligned <- data.frame(not_aligned=unmapped[[1]][,"__not_aligned"][match(rownames(unmapped[[1]]), rownames(total_reads))])
total_reads$not_aligned <- not_aligned$not_aligned

# subtract the number of reads that aligned from the total number of reads 
total_reads <- transform(total_reads, aligned=total_reads-not_aligned)

# create a vector of the numbers of reads aligned
total_aligned <- t(as.matrix(subset(total_reads, select="aligned")))

# divide the numbers of reads aligned by one million
total_aligned <- total_aligned/1000000

RPKMdivisor <- tag_lengths %*% total_aligned

save(RPKMdivisor, file=paste(output_dir,"RPKM_divisor.RData", sep="/"))

######################################################################
# Step 4: Create RPKM table
######################################################################

RPKM <- RPKMdividend/RPKMdivisor
save(RPKM, file=paste(output_dir,"RPKM.RData", sep="/"))
