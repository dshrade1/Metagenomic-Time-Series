wd <- "~/Desktop/Metagenome_TS/0_minid"
percent_reads_dir <- "percent_reads"
setwd(wd)

# read in file that lists the percent of reads that mapped for each minid value
reads_mapped_files <- list.files(percent_reads_dir)

setwd(percent_reads_dir)

minid_seq <- as.character(format(seq(0.70,1.0,0.01), nsmall=2))

#create a list of 31 tables, where each table represents a single minid value
read.mapped.tables <- function(x) read.table(reads_mapped_files[[grep(x,reads_mapped_files)]], col.names=c("mapped","pct.reads", "num.reads", "pct.bases", "num.bases"))
reads_mapped_tables <- sapply(minid_seq, read.mapped.tables, simplify = FALSE) #almost instantaneous
# using simplify=T could be handy in some situations, too

# create a table that represents percent reads mapped per minid value
# create function that extracts first column
select_pct.reads <- function(x) subset(reads_mapped_tables[[x]], select="pct.reads")
pct.reads_table <- data.frame(sapply(minid_seq, select_pct.reads))

# convert percentage values to numbers
pct2numeric <- function(x) as.numeric(substr(as.character(x), 1,nchar(as.character(x))-1))
pct.reads_table <- data.frame(apply(pct.reads_table, c(1,2), pct2numeric))
# calculate mean values per minid value
mean_pct.reads.mapped <- apply(pct.reads_table,2, mean)

min_id <- seq(0.70,1.0,0.01)




# create a table that represents percent bases mapped per minid value
# create function that extracts first column
select_pct.bases <- function(x) subset(reads_mapped_tables[[x]], select="pct.bases")
pct.bases_table <- data.frame(sapply(minid_seq, select_pct.bases))
pct.bases_table <- data.frame(apply(pct.bases_table, c(1,2), pct2numeric))
# calculate mean values per minid value
mean_pct.bases.mapped <- apply(pct.bases_table,2, mean)


# plot in base R
plot(min_id, mean_pct.reads.mapped, main="Mean Percent of Reads Mapped", xlab="min_id value", ylab="mean percent of reads mapped")
points(min_id, mean_pct.bases.mapped, main="Mean Percent of Bases Mapped", xlab="min_id value", ylab="mean percent of bases mapped", col="red")
legend("topright", c("Reads Mapped", "Bases Mapped"),col=c("black", "red"), pch=1)

# plot in ggplot
library(ggplot2)
ggplot() +
geom_point(aes(x = min_id, y = mean_pct.reads.mapped, color = "Reads Mapped"))  +
  geom_point(aes(x = min_id, y = mean_pct.bases.mapped, color = "Bases Mapped"))  +
  xlab('BBMap min_id value') +
  ylab('Percent Mapped') +
  labs(color="Value Calculated") + 
  ggtitle("Percent Mapped of a Random Subset of TH Reads") +
  scale_x_continuous(breaks = seq(0.7, 1.0, by=0.01))
  #theme(legend.key.height=unit(5.35,"line")) +
  #scale_fill_discrete("") +
  #guides(fill = guide_legend(keyheight = 10)) +
  #xlim(as.Date("04-15-2008", format = "%m-%d-%Y"),as.Date("11-15-2008", format = "%m-%d-%Y")) +
  #theme(axis.ticks = element_blank(), axis.text.y = element_blank())
  


# find a curve to fit this plot:

x <- min_id
y <- mean_pct.reads.mapped
x1 <- min_id[which(min_id>0.92)]
y1 <- mean_pct.reads.mapped[which(min_id>0.92)]
fit4 <- lm(y1~poly(x1,4,raw=T))
xx <- seq(0.93,1.0,0.01)
plot(x,y,pch=19, main = "Fitting a curve to percent mapped to find optimal min_id value", xlab="% reads mapped", ylab="min_id values tested")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple")

x2 <- min_id[which(min_id<=0.92)]
y2 <- mean_pct.reads.mapped[which(min_id<=0.92)]
fit2 <- lm(y2~poly(x2, 2, raw=T))
zz <- seq(0.70,0.92,0.01)
lines(zz, predict(fit2, data.frame(x=zz)), col="green")

