# counts_analysis.R
# analysis of the feature count output from htseq-count

# Step 1:  Create tables summarizing #reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
# Step 2:  Create table summarizing %reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
# Step 3:  Determine total % of reads mapped in each metagenome
# Step 4:  Which dates have the highest % of total reads mapping, and which have the lowest? Plot.
# Step 5:  How many reads mapped as each feature type on each date? Plot as stacked barplot.
# Step 6:  Are the ratios of [CDS + aligned noncoding]:total reads consistent?
# Step 7:  For samples taken on the same date, how similar are the percents of reads mapped? total number of reads input into mapping (metagenome size)?
# Step 8:  Sanity checks

#######################################################################################################
# Step 1:  Create tables summarizing #reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")

load("mapped_counts.RData")
load("unmapped_counts.RData")
load("total_read_counts.RData")

# CREATE A DATA FRAME OF the SUM of READS MAPPED PER METAGENOME
mapped_sums <- sapply(feature_type, function(ft) as.data.frame(apply(mapped[[ft]],2,sum)), simplify=F)

create.summary.table <- function(ft) {
	x <- cbind(unmapped[[ft]], 
	sum_mapped=mapped_sums[[ft]],
	#total_reads=sapply(feature_type, function(ft) as.data.frame(apply(df[[ft]],2,sum)), simplify=F)[,1] )
	#total_reads=as.data.frame(apply(df[[ft]],2,sum)) )
	total_reads=total_reads)
	names(x) <- c("no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique", "sum_mapped", "total_reads")
	return(x)
}

summary_table <- sapply(feature_type,create.summary.table, simplify=F)

summary_numbers <- data.frame(sapply(feature_type, function(ft) apply(summary_table[[ft]],2,sum) ))

save(summary_numbers, file="summary_numbers.RData")
save(summary_table, file="summary_table.RData")
summary_numbers

#######################################################################################################
# Step 2:  Create table summarizing %reads/metagenome that mapped as CDS, repeat_region, rRNA, tRNA, etc.
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")


load("summary_numbers.RData") # loads summary_numbers

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
quartz(width=13)
pie(summary_percents[1,], labels= c(
"% CDS", # pct_CDS
"", # pct_ambiguous_CDS
"", # pct_rr
"", # pct_rRNA
"% rRNA, tRNA, repeat_region, or ambiguous", #pct_tRNA
"% not aligned", 
"% aligned noncoding" ), col=cm.colors(7))

quartz.save(file="Pie Chart Overall Mapping.jpg", type="jpg", dpi=200)

#######################################################################################################
# Step 3:  Determine total % of reads mapped in each metagenome
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")

load("summary_table.RData") # loads summary_table

# create a vector the total sum of reads mapped (CDS + rr + rRNA + tRNA)
mapped_sums_all <- apply(data.frame(lapply(summary_table,`[`, "sum_mapped")), 1, sum)
# divide this by the total number of reads 
reads_mapped <- cbind(mapped=mapped_sums_all, total_reads=total_reads)
reads_mapped <- transform(reads_mapped,percent_mapped=mapped/total_reads)
quartz(width=13)
hist(reads_mapped$percent_mapped, breaks=15, col="gray", main="Percent of Reads/Metagenome that Mapped to Coassembly as Feature", xlab="% of Total Reads Mapped")
quartz.save("Distribution of Percent Reads Mapped.jpg", type="jpg", dpi=200)
range(reads_mapped$percent_mapped) #0.3037097 0.5284729

setwd(wd)
save(reads_mapped, file="reads_mapped.RData")

#######################################################################################################
# Step 4:  Which dates have the highest % of total reads mapping, and which have the lowest? Plot.
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")


load("reads_mapped.RData")

dates <- read.csv("~/Desktop/Metagenome_TS/6_analysis/TH_Metagenome_Dates.csv")
dates$date <- as.Date(dates$Date, format="%m/%d/%y")
dates$Date <- NULL
save(dates, file="dates.RData")

reads_mapped_dates <- reads_mapped

reads_mapped_dates$Library.ID <- substring(row.names(reads_mapped_dates), 1,5)
# in this merge, I want to get rid of non-required Library IDs but keep the duplicates

reads_mapped_dates <- merge(reads_mapped_dates, dates)

quartz(width=13)
ggplot() + 
  geom_point(data = reads_mapped_dates, aes(x = date, y = percent_mapped, color="red")) + # , color = "Trout Bog Hypo metagenomes"))  +
   theme_gray() +
   theme(legend.position="none") +
   xlab("Date") +
   ylab("Percent of Reads that Mapped")
   ggtitle("Percent of Reads Mapped as CDS, rRNA, tRNA, or repeat region")
quartz.save("Percent Reads Mapped by Date.jpg", type="jpg", dpi=200)
save(reads_mapped_dates, file="reads_mapped_dates.RData")

# Finding: each year's % mapped as feature is different

#######################################################################################################
# Step 5:  How many reads mapped as each feature type on each date? Plot as stacked barplot.
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")

load("summary_table.RData")
load("dates.RData")
load("reads_mapped_dates.RData")

# pull different columns from the 4-member summary_table list to create the table that contains the components I want in the stacked bar chart.

# stacks in the stacked bar chart: 

summary_df <- data.frame(summary_table)
summary_df <- subset(summary_df,select=c(
grep("CDS.not_aligned", names(summary_df)),
grep("sum_mapped", names(summary_df)),
grep("CDS.ambiguous", names(summary_df)),
grep("CDS.total_reads", names(summary_df))
))
names(summary_df)[1] <- "not_aligned"
names(summary_df)[ncol(summary_df)] <- "total_reads"
summary_df <- transform(summary_df, aligned_noncoding=total_reads-(not_aligned+CDS.sum_mapped+rr.sum_mapped+rRNA.sum_mapped+tRNA.sum_mapped+CDS.ambiguous))
summary_df$total_reads <- NULL
summary_df$Library.ID <- substring(row.names(summary_df), 1,5)
summary_df <- merge(summary_df, dates)

# collapse by date
# http://stackoverflow.com/questions/2558191/how-can-i-collapse-a-dataframe-by-some-variables-taking-mean-across-others
summary_df <- aggregate(summary_df,list(summary_df$Library.ID), mean)
# warnings ok here.
summary_df$Library.ID <- NULL
# sort summary_df by date; reorder columns so that it's got the following order: 
# top : tRNA, rRNA, rr, CDS, CDS ambiguous, aligned_noncoding, not_aligned : bottom
summary_df <- summary_df[order(summary_df$date),c(9,1,6,5,4,3,7,8,2)]
summary_df$Group.1 <- as.factor(summary_df$Group.1)
row.names(summary_df) <- 1:nrow(summary_df)
save(summary_df, file="summary_df.RData")

# reorder it # reshape data from wide to long for plotting
# http://stackoverflow.com/questions/25664007/reorder-bars-in-geom-bar-ggplot2
summary_long <- melt(summary_df, id.vars=c("date","Group.1"))

total_input <- data.frame(apply(summary_df[,3:ncol(summary_df)], 1, sum))
total_input <- cbind(date=summary_df$date , total_input =data.frame(total_input))
names(total_input)[2] <- "total_input"
save(total_input, file="total_input.RData")

# create a data matrix of percents that fall into each category per date
summary_prop <- prop.table(data.matrix(summary_df[,3:ncol(summary_df)]), 1)
summary_prop  <- cbind(Group.1=summary_df$Group.1, date=summary_df$date , data.frame(summary_prop))
summary_prop_long <- melt(summary_prop, id.vars=c("date", "Group.1"))

# plot stacked barplot of different mapping types per Library.ID
stacked <- ggplot(summary_long, aes(x = reorder(Group.1, date), y = value,fill=variable)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("metagenome") + ylab("reads in metagenome") + ggtitle("Mapping Efficiency by Date")
# plot barplot of all reads input into mapping per Library.ID
mapped <- ggplot(reads_mapped_dates, aes(x = reorder(Library.ID, date), y = mapped)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("metagenome") + ylab("reads")
# the following plot gives what we want: stacked barplot of different mapping types per date
mapped_with_date <- ggplot(summary_long, aes(x = reorder(as.factor(date), date), y = value,fill=variable)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("reads") + ggtitle("Mapping Efficiency by Date")
    
input_with_date <- ggplot(total_input, aes(x = reorder(as.factor(date), date), y = total_input)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("reads") + ggtitle("Metagenome Size (in # of Reads) by Date")    
    
stacked_prop <- ggplot(summary_prop_long, aes(x = reorder(as.factor(date), date), y = value,fill=variable)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("percent of reads in metagenome") + ggtitle("Mapping Efficiency by Date")

quartz(width=13)
input_with_date
quartz(width=13)
mapped
quartz.save("Reads Input into Mapping by Date.jpg", type="jpg", dpi=200)
quartz(width=13)
mapped_with_date
quartz.save("Stacked Mapping Efficiency - reads.jpg", type="jpg", dpi=200)
quartz(width=13)
stacked_prop
quartz.save("Stacked Mapping Efficiency by Date - proportions.jpg", type="jpg", dpi=200)

#######################################################################################################
# Step 6:  Are the ratios of [CDS + aligned noncoding]:total reads consistent?
#######################################################################################################

# the big question here is: is there some number that is constant?
# is the % of ALL the reads that mapped proportionally low compared to the amount of reads that mapped as 
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")

load("summary_df.RData")
load("total_input.RData")

# what is the total % aligned over time (not just aligned as feature)?
# hypothesis: not consistent/dependent on what got into the assembly
summary_aligned <- mutate(summary_df, aligned=tRNA.sum_mapped+rRNA.sum_mapped+rr.sum_mapped+CDS.sum_mapped+CDS.ambiguous+aligned_noncoding)
summary_aligned <- summary_aligned[,c(1,2,ncol(summary_aligned),ncol(summary_aligned)-1)]

aligned_prop <- prop.table(data.matrix(summary_aligned[,3:ncol(summary_aligned)]), 1)
aligned_prop  <- cbind(Group.1=summary_aligned$Group.1, date=summary_aligned$date , data.frame(aligned_prop))
aligned_prop_long <- melt(aligned_prop, id.vars=c("date","Group.1"))

aligned <- ggplot(aligned_prop_long, aes(x = reorder(as.factor(date), date), y = value,fill=variable)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("reads") + ggtitle("Proportions of Metagenome Reads that Aligned to Coassembly")
aligned
# Hypothesis is supported: inconsistent number of reads


# what is the % of the # of reads aligned over time that are feature?
# hypothesis: the % of reads that aligned which align as feature is consistent over time.

summary_feature <- mutate(summary_df, feature=tRNA.sum_mapped+rRNA.sum_mapped+rr.sum_mapped+CDS.sum_mapped+CDS.ambiguous)
summary_feature <- summary_feature[,c(1,2,ncol(summary_feature),ncol(summary_feature)-2)]
feature_prop <- prop.table(data.matrix(summary_feature[,3:ncol(summary_feature)]), 1)
feature_prop  <- cbind(Group.1=summary_feature$Group.1, date=summary_feature$date , data.frame(feature_prop))
feature_prop_long <- melt(feature_prop, id.vars=c("date","Group.1"))

prop_feature <- ggplot(feature_prop_long, aes(x = reorder(as.factor(date), date), y = value,fill=variable)) +
    geom_bar(stat='identity') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("reads") + ggtitle("Proportions of Metagenome Reads that Aligned to Coassembly which are Features")
    
quartz(width=13)
prop_feature
quartz.save("Proportion of aligned reads that are features.jpg", type="jpg", dpi=200)
range(feature_prop$feature)[2] - range(feature_prop$feature)[1] #0.06941593
range(feature_prop$aligned_noncoding)[2] - range(feature_prop$aligned_noncoding)[1] #0.06941593

# my hypothesis is supported: the proportion of reads that map as feature is fairly steady.
# the 7% could be natural variation in the % that aligns, as opposed to the percent
# So maybe we should use the number that aligned as the denominator.


# we can say:
# of all the reads input into the mapping, the proportion that mapped to the combined assembly was inconsistent.
# this is due to how the combined assembly favored certain genomes and not others.
# but if we take as our baseline, as our denominator, the number of reads that did align,
# then the number of reads annotated as feature is much much more consistent
# and it's not based on how well certain metagenomes were represented in the combined assembly.


#######################################################################################################
# Step 7:  For samples taken on the same date, how similar are the percents of reads mapped? total number of reads input into mapping (metagenome size)?
#######################################################################################################
rm(list = ls())
setwd("~/Desktop/Metagenome_TS/3_count")
source("counts_analysis_source.R")

load("reads_mapped_dates.RData")

duplicate_dates <- reads_mapped_dates[c(which(duplicated(reads_mapped_dates$Library.ID)),which(duplicated(reads_mapped_dates$Library.ID))-1),]
date_comp <- duplicate_dates[order(duplicate_dates$date),]
date_comp$version <- rep(c(1,2))

library(lattice)
barchart(total_reads~as.factor(date),data=date_comp,groups=version,scales=list(x=list(rot=90,cex=0.8)), xlab="sample date", ylab="total reads", main="Metagenome sizes for duplicate samples taken on same date")
quartz.save("Duplicates_reads.jpg",type="jpg", dpi=200)

barchart(percent_mapped~as.factor(date),data=date_comp,groups=version,scales=list(x=list(rot=90,cex=0.8)), xlab="sample date", ylab="Percent mapped as feature", main="% mapped as feature for duplicate samples taken on same date")
quartz.save("Duplicates_percents.jpg",type="jpg", dpi=200)

# Finding: They are Very tightly linked.

#######################################################################################################
# Step 8:  Sanity checks
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

