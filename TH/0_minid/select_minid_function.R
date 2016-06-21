####### This function extracts data on % mapped, % mapped perfectly, and % mapped with error from '.err' outputs from BBMap #############
####### Then the plotting scheme plots #############

# load libraries
library(ggplot2)

# user inputs
wd <- "~/Desktop/Metagenome_TS/0_minid"
percent_reads_dir <- "percent_reads"
percent_perfect_dir <- "percent_perfect"
percent_error_dir <- "percent_error"
minid_seq <- as.character(format(seq(0.70,1.0,0.01), nsmall=2))
min_id <- seq(0.70,1.0,0.01)

# function to quantify the various percents of reads mapped
quantify_mapping <- function(mapping_output_dir) { # where 'mapping_output_dir' is either percent_reads_dir or percent_perfect_dir
	setwd(wd)
	reads_mapped_files <- list.files(mapping_output_dir)
	setwd(mapping_output_dir)
	print(mapping_output_dir)
	#create a list of 31 tables, where each table represents a single minid value
	read.mapped.tables <- function(x) read.table(reads_mapped_files[[grep(x,reads_mapped_files)]] ) 
	reads_mapped_tables <- sapply(minid_seq, read.mapped.tables, simplify = FALSE) #almost instantaneous
	select.relevant.columns <- function(x) subset(x,select=((ncol(x)-3):(ncol(x))))
	reads_mapped_tables <- sapply(reads_mapped_tables, select.relevant.columns, simplify = FALSE)
	# change names of all tables in the list reads_mapped_tables
	rename.columns <- function(x) {
		x <- data.frame(x)
		names(x) <- c("pct.reads", "num.reads", "pct.bases", "num.bases")
		return(x)
		}
	reads_mapped_tables <- lapply(reads_mapped_tables, rename.columns)
	print(head(reads_mapped_tables[[1]]))
	# create a table that represents percent reads mapped per minid value
	select_pct.reads <- function(x) subset(reads_mapped_tables[[x]], select="pct.reads")
	pct.reads_table <- data.frame(sapply(minid_seq, select_pct.reads))
	#print(head(pct.reads_table))
	# convert percentage values to numbers
	pct2numeric <- function(x) as.numeric(substr(as.character(x), 1,nchar(as.character(x))-1))
	pct.reads_table <- data.frame(apply(pct.reads_table, c(1,2), pct2numeric))
	print(head(pct.reads_table))
	# calculate mean values per minid value
	mean_pct.reads.mapped <- apply(pct.reads_table,2, mean)
	# ^ looks good
	# create a table that represents percent bases mapped per minid value
	select_pct.bases <- function(x) subset(reads_mapped_tables[[x]], select="pct.bases")
	pct.bases_table <- data.frame(sapply(minid_seq, select_pct.bases))
	pct.bases_table <- data.frame(apply(pct.bases_table, c(1,2), pct2numeric))
	print(head(pct.bases_table))
	# calculate mean values per minid value
	mean_pct.bases.mapped <- apply(pct.bases_table,2, mean)
	# make a quick base-R plot to ensure that this is the outcome we want, but don't save it
	plot(min_id, mean_pct.bases.mapped, main="Mean Percent of Reads or Bases Mapped", xlab="min_id value", ylab="mean percent mapped", ylim=c(0,100))
	points(min_id, mean_pct.reads.mapped, col="red")
	legend("topright", c("Reads Mapped", "Bases Mapped"),col=c("black", "red"), pch=1)
	return(list(mean_pct.reads.mapped, mean_pct.bases.mapped))
	rm(mapping_output_dir)
	setwd(wd)
}


# plot all plots on the same plot

setwd(wd)
p <- ggplot() +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_reads_dir)[[1]], color = "Percent of All Reads that Mapped"))  +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_reads_dir)[[2]], color = "Percent of All Bases that Mapped"))  +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_perfect_dir)[[1]], color = "Percent of Reads that Mapped Perfectly"))  +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_perfect_dir)[[2]], color = "Percent of Bases that Mapped Perfectly"))  +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_error_dir)[[1]], color = "Percent of Reads that Mapped with Error"))  +
	geom_point(aes(x = min_id, y = quantify_mapping(percent_error_dir)[[2]], color = "Percent of Bases that Mapped with Error"))  +
	xlab('BBMap min_id value') +
	ylab('Percent') +
	labs(color="Value Calculated") + 
	ggtitle(paste("Mapping Efficiencies of a Random Subset of TH Reads")) +
	scale_x_continuous(breaks = seq(0.7, 1.0, by=0.01)) + 
	ylim(0,100)
	p
	setwd(wd)
	

savePlot <- function(myPlot) {
        #pdf(paste("Plot_", mapping_output_dir ,".pdf"))
        pdf("Plot.pdf")
        print(myPlot)
        dev.off()
}

	savePlot(p)
	
# optional find a curve to fit this plot:

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
