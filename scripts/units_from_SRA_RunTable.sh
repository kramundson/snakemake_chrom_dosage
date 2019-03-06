#!/bin/bash

# table for paired end read data, not used for karyotype analysis but nice to have on hand
tail -n +2 Tang_etal_2018_GenomBiol_SRA.txt | \
	awk -v OFS="\t" 'BEGIN {print "sample","unit","fq1","fq2","parhap"} {print $2,$6,$6"_1.fastq.gz",$6"_2.fastq.gz","haploid"}' > units_PE.tsv

# table for forward mates only, used for karyotype analysis
tail -n +2 Tang_etal_2018_GenomBiol_SRA.txt | \
	awk -v OFS="\t" 'BEGIN {print "sample","unit","fq1","fq2","parhap"} {print $2,$6,$6"_1.fastq.gz","NaN","haploid"}' > units.tsv
