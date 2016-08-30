library(ggplot2)
library(plyr)

# analysis of RPKM output
normalize_dir <- "../4_normalize"
counts_dir <- "../3_count"
output_dir <- paste(site,minid,"output",sep="_")
data_dir <- paste("../raw_data", site, sep="/")


setwd(normalize_dir)
load(paste(output_dir,"RPKM.RData", sep="/"))

#######################################################################################################
# what are the sums of the RPKM values per metagenome? Are they consistent?
#######################################################################################################

RPKM_sums <- data.frame(rowSums(data.frame(t(RPKM)))) # slow
names(RPKM_sums) <- "sum"
RPKM_sums$Library.ID <- substr(row.names(RPKM_sums), 1, 5)

RPKM_values <- ggplot(RPKM_sums, aes(x = Library.ID, y = sum)) + geom_bar(stat='identity')+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("metagenome") + ylab("RPKM")

# sort these bars by date
setwd(counts_dir)
load(paste(output_dir,"dates.RData", sep="/"))
setwd(normalize_dir)
RPKM_sums_dates <- merge(RPKM_sums, dates)
RPKM_sums_dates <- aggregate(RPKM_sums_dates,list(RPKM_sums_dates$Library.ID), mean) #warnings ok here
RPKM_sums_dates$Library.ID <- NULL
RPKM_sums_dates$Group.1 <- as.factor(RPKM_sums_dates$Group.1)
RPKM_sums_dates$Group.1 <- substr(RPKM_sums_dates$Group.1, 1,4)
RPKM_sums_dates <- ddply(RPKM_sums_dates, .(date,Group.1), numcolwise(mean))
save(RPKM_sums_dates,file=paste(output_dir,"RPKM_sums_dates.RData", sep="/"))

RPKM_values_dates <- ggplot(RPKM_sums_dates, aes(x = reorder(as.factor(date), date), y = sum)) + geom_bar(stat='identity')+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("metagenome") + ylab("RPKM")

jpeg(filename=paste(output_dir,"Sum of RPKM by Date.jpg",sep="/"), width=1300, height= 750, quality=100)
print(RPKM_values_dates)
dev.off()
setwd(normalize_dir)
#quartz.save("Sum of RPKM by Date.jpg", type="jpg", dpi=200)

# Finding: this generally mirrors the plot of proportions of reads that aligned to coassembly
# (found in count_dir, "Proportion of aligned reads that are features.jpg")
# In other words, normalization appears to have been successful. 
# Because of the relative consistency, we should be able to compare relative abundance of RPKM.
# Perhaps the small variations are due to the normalization by gene length?
# Is there a way that we can quantify the semblance?
