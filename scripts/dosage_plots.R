# Script computes normalied read depth based on a single control sample,
# then generates a plot similar to the ones we routinely use for illustrating
# chromosome dosage in JMP

# Only a single sample is processed at a time. This enables parallelization over a cluster, if desired
# Note, this approach requires comparison between case-control pairs.
# If you want to use the population-wide average as a control, use bin by sam

# USAGE: Rscript dosage_plots.R <controlfile> <testfile> <output plot>

library(ggplot2)
library(dplyr)
# library(grid)
# library(gridExtra)

args = commandArgs(trailingOnly = TRUE)

# define header for BED format table output from bedtools coverage
header.temp <- c("chrom","start","end","readcount","bases_covered", "binsize", "breadth")

# read in control
fh.con <- read.csv(args[1],
                   header=F,
                   sep = '\t',
                   comment.char = "#",
                   col.names = header.temp)

# read in a single bedtools coverage output table
fh <- read.csv(args[2],
               header=F,
               sep='\t',
               comment.char="#",
               col.names=header.temp)

# Compute normalized coverage
fh$normcov <- 2*(fh$readcount/sum(fh$readcount)) / (fh.con$readcount / sum(fh.con$readcount))

# Filter out non-pseudomolecule scaffolds
chrs <- c("chr01","chr02","chr03","chr04","chr05","chr06",
          "chr07", "chr08", "chr09", "chr10", "chr11","chr12")

fh.filtered <- subset(fh, chrom %in% chrs)
fh.filtered$chrom <- factor(fh.filtered$chrom)

# Define functions for plot aesthetics

mid.dosage <- function(x,y) { # x, y is title of plot
  n=floor(nrow(x)/2)
  return(x$bin[n])
}

plot.dosage <- function(x,y) {
  numblanks.dosage <- 15
  stuf.d <- c(rep(NA, numblanks.dosage))
  dosage.stuffer <- data.frame("chrom"=stuf.d, "start"=stuf.d, "end"=stuf.d,
                               "readcount"=stuf.d, "bases_covered"=stuf.d,
                               'binsize'=stuf.d,"breadth"=stuf.d,"normcov"=stuf.d,
                               "bin"=stuf.d)
  x$bin <- seq(1,nrow(x))
  chr.list.dosage <- split(x, f=x$chrom)
  chr.list.dosage.stuffed <- lapply(chr.list.dosage[1:11], function(x) rbind(x,dosage.stuffer))
  chr.list.dosage.stuffed <- dplyr::bind_rows(chr.list.dosage.stuffed, chr.list.dosage[12])
  chr.list.dosage.stuffed$bin2 <- seq(1:nrow(chr.list.dosage.stuffed))
  chr.list.dosage.stuffed$normcov <- as.numeric(chr.list.dosage.stuffed$normcov) # force as numeric if not
  
  midpoints.dosage <- sapply(chr.list.dosage[1:12],mid.dosage)
  
  plt.dosage <- ggplot(chr.list.dosage.stuffed, aes(x=bin2,y=normcov)) +
    labs(x="",y="Copy Number") +
    geom_line(color="#008080",fill="white",size=1) + # color currently hardcoded
    ggtitle(y) +
    geom_point(aes(x=bin2,y=normcov), size=1.0, color="black") +
    guides(fill=F,color=F) +
    scale_x_continuous(breaks=which(chr.list.dosage.stuffed$bin %in% midpoints.dosage),
                       labels=names(midpoints.dosage)) +
    scale_y_continuous(limits=c(0.5,5), breaks=seq(1,5), labels=seq(1,5)) +
    theme(panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_line(color="black",linetype="dashed"),
          panel.background=element_rect(fill="white",color="black"),
          axis.text.x=element_text(size=12,color="black"),
          axis.text.y=element_text(size=12,color="black"),
          axis.title.y=element_text(size=12,angle=90,vjust=0.5),
          axis.ticks=element_blank(),
          plot.title=element_text(size=12,face="bold",hjust=0))
  write.table(chr.list.dosage.stuffed, gsub(".pdf", ".tsv", args[3]), 
              quote=F,eol='\n',row.names=F)
  ggsave(args[3], width=10,height=2,units="in",plot=plt.dosage,device="pdf")
}

# modify plot title
title <- gsub(".+/[0-9]x_([A-Za-z0-9]+)", "\\1", args[2])
plot.dosage(fh.filtered,title)