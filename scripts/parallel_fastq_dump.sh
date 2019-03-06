#!/bin/bash
source activate dosage
parallel -j 6 fastq-dump --gzip --split-3 -B --readids -O data/reads ::: sra/*.sra
