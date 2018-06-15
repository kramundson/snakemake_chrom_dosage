# Script computes normalized read depth based on a single control sample.
# File names of individuals to be analyzed are specified in bedfiles.txt
# Output file overlays chromosome dosage for all individuals in a population
# File naming convention is handled by Snakemake

# TODO: add back in ability to have control be population mean. Will be an optional command line argument that defaults to FALSE
# Then, write function that computes per-bin depth relative to population mean

# USAGE: Rscript dosage_plots_population.R <controlfile> bedfiles.fofn <output>

library(ggplot2)
library(dplyr)
args=commandArgs(trailingOnly = T)

# define header for BED format output file from bedtools coverage
header.temp <- c("chrom","start","end","readcount","bases_covered", "binsize", "breadth")

# read in control file
fh.con <- read.csv(args[1],
                   header=F,
                   sep='\t',
                   comment.char="#",
                   col.names=(header.temp))

# Read in all test files using bedfiles.fofn as file of file names
filenames <- read.csv(args[2], header=F) # will be used in snakemake
input <- apply(filenames, 1, function(x) read.csv(x, sep='\t', comment.char="#", header=F, col.names=header.temp))

# Biological sample names are too long. Chop off leading and trailing extraneous text
names(input) <- gsub("-windowcov.bed", "", filenames$V1)
names(input) <- gsub("^data/bedtools_coverage/", "", names(input))

# Compute normalized per-bin coverage and add to each read-in sample dataframe
for (i in seq(1, length(input))) {
  input[[i]] <- cbind(input[[i]], "sample"=names(input)[i], "chrbin"=floor(input[[i]]$start/input[[i]]$end[1]),
                      "normcov"=2*(input[[i]]$readcount/sum(input[[i]]$readcount)) /(fh.con$readcount / sum(fh.con$readcount)))
}

# Concatenate all sample data frames together
df <- dplyr::bind_rows(input)

# Filter out non-pseudomolecule scaffolds. Needed for mapping but do not want to plot.
chrs <- c("chr01","chr02","chr03","chr04","chr05","chr06",
          "chr07", "chr08", "chr09", "chr10", "chr11","chr12")
df.filtered <- subset(df, chrom %in% chrs)

# Reassign df.filtered$chrom as factor to remove chrom factor levels that were filtered out
df.filtered$chrom <- factor(df.filtered$chrom)

# Plot and save output to disk
plt <- ggplot(df.filtered, aes(x=chrbin,y=normcov)) +
  geom_line(aes(fill=sample), color="black",alpha=0.4) +
  facet_wrap(~ chrom, scales="free_x") +
  guides(color = F) +
  scale_y_continuous(limits=c(0.5,4))
ggsave(args[3], width=16, height=9, units="in",plot=plt,device="pdf")

write.table(df.filtered, file="pop_dosag.tsv", quote = F, sep = '\t', eol = '\n', na = "NA", row.names = FALSE, col.names = TRUE)